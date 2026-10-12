# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
# Copyright (C) 2026 Jack S. Hale
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Domains.

A domain is a set over which something can be integrated.

Some domains are coordinate domains: their points are tuples of numbers,
so a tuple names a point. Others are parametrized domains, presented as
the image of a map out of a coordinate domain; their points are not
tuples, and naming one means naming a tuple in the source domain and
composing.

That map is a parametrization. How it is described is not fixed: a finite
element basis summed against a mesh's coordinate dofs is one description
and a closed form expression is another, and two of them compose. A mesh
is a parametrized domain whose parametrization is a finite element one,
and what UFLx holds is the shape of that parametrization rather than the
map itself, since the coordinate dofs stay external.

In finite element terms, the two kinds of domain are the reference cell
and the mesh. For a surface mesh of triangles in three dimensions::

    mesh = parametrized_domain(P1_vector)       # tdim 2, gdim 3
    X = Point([a, b], entity_domain(triangle))  # tdim == gdim == 2
    x = ImagePoint(X, mesh).expand_geometry()  # a point of RD(3)

The pair (a, b) names no point of the mesh: it has many triangles, and a
point of the surface needs three ambient coordinates. It names a point of
the triangle's coordinate domain, and the push forward carries it through
one cell's parametrization to get ambient coordinates.

A parametrized domain makes no assumption that it only contains cells of
a single type: one could contain (eg) a mixture of triangles and
quadrilaterals, or even a mixture of (eq) tetrahedra and intervals.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

from uflx.entities import AbstractEntity
from uflx.tensors import Identity, Vector

if TYPE_CHECKING:
    from uflx.expressions import AbstractExpression
    from uflx.variables import AbstractVariable


class AbstractDomain(ABC):
    """Base class for a domain."""

    @property
    @abstractmethod
    def geometric_dimension(self) -> int:
        """The number of coordinates needed to name a point of this domain."""

    @property
    @abstractmethod
    def topological_dimension(self) -> int | None:
        """The topological dimension of the domain.

        This returns None iff the domain contains entities of a mixture
        of topological dimensions.
        """

    @property
    def coordinate_domain(self) -> AbstractCoordinateDomain:
        """The coordinate domain this domain's points are named in.

        How many coordinates name a point is the geometric dimension;
        which coordinates those are is this. They are the ambient ones
        unless a domain says otherwise, and a facet of a reference cell
        is where the difference shows: its points are named in the cell's
        coordinates and not in R^d.
        """
        return RD(self.geometric_dimension)


class AbstractCoordinateDomain(AbstractDomain):
    """Base class for a domain whose points are coordinate tuples.

    A tuple of numbers names a point of such a domain outright, so tuple
    equality is point equality and arithmetic on components is
    meaningful. Its geometry is the identity, so its topological and
    geometric dimensions agree and a subclass need only give one of them.
    """

    @property
    def topological_dimension(self) -> int:
        """The topological dimension of the domain.

        The geometry is the identity, so this is the geometric dimension.
        """
        return self.geometric_dimension

    @property
    def coordinate_domain(self) -> AbstractCoordinateDomain:
        """A coordinate domain's points are named in its own coordinates."""
        return self


class AbstractParametrizedDomain(AbstractDomain):
    """Base class for a domain presented through maps out of coordinate domains.

    A parametrization is a map out of a parameter region: a set whose
    points are coordinate tuples, carried by the map into this domain. The
    region is what identifies the map, since naming a point of this domain
    means naming a point of some region and composing.

    A reference cell is one kind of parameter region and not the only
    kind, so being parametrized does not imply being made of cells. A
    region of R^d parametrized by itself is a domain of this kind with no
    cells in it, and its boundary is another. What a geometric quantity needs of a
    domain is a map out of the region its point lies in, which is why it
    asks for one here rather than for a cell's.
    """

    @property
    def coordinate_domain(self) -> AbstractCoordinateDomain:
        """The coordinate domain this domain's points are named in.

        Where its parametrizations land. They must agree, a domain's points being
        named one way.

        Raises:
            ValueError: If its parametrizations land in different coordinates
        """
        targets = {parametrization.target for parametrization in self.parametrizations}
        if len(targets) != 1:
            raise ValueError(
                f"The parametrizations of {self!r} land in {len(targets)} coordinates."
            )
        (target,) = targets
        return target

    @property
    @abstractmethod
    def parametrizations(self) -> tuple[AbstractParametrization, ...]:
        """The maps this domain is parametrized by, each knowing its own region.

        Not keyed by that region. Several parametrizations may start from one of
        them -- every facet of a triangle is parametrized from the reference
        interval -- so the region does not name the parametrization, and a domain
        is not asked to choose between them. Where exactly one parametrization is
        needed, as when a geometric quantity is evaluated at a point of
        some region, that is a question about the quantity rather than
        about the domain, and `uflx.geometry` asks it.
        """


class RD(AbstractCoordinateDomain, AbstractParametrizedDomain):
    """R^d, the ambient coordinate domain.

    Parametrized by itself, through the identity: its points are coordinate
    tuples already, so nothing has to carry them anywhere. The measure
    that follows is the Lebesgue one, whose density is one.
    """

    def __init__(self, dim: int):
        """Initialise."""
        self._dim = dim

    @property
    def geometric_dimension(self) -> int:
        """The number of coordinates needed to name a point of this domain."""
        return self._dim

    @property
    def parametrizations(self) -> tuple[AbstractParametrization, ...]:
        """R^d is parametrized by itself, through the identity."""
        return (IdentityParametrization(self),)

    def __repr__(self) -> str:
        """Representation."""
        return f"RD({self._dim})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, RD) and self._dim == other._dim

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.RD", self._dim))


class AbstractParametrization(ABC):
    """Base class for a map out of a coordinate domain.

    A parametrization takes a point of its source, whose points are
    coordinate tuples, and gives the coordinates it lands on. How the map
    is described is the implementer's business, so this class asks only
    that it can be evaluated and differentiated at a point.

    Two parametrizations compose when the first lands where the second
    starts, which is why the source is a domain and not a dimension.

    Implementations must be comparable by value and hashable: a
    parametrization is an argument of the geometric quantities built on
    it, so rewriting an expression and simplifying a product of Jacobians
    both compare them.
    """

    @property
    @abstractmethod
    def source(self) -> AbstractCoordinateDomain:
        """The coordinate domain this map starts from."""

    @property
    @abstractmethod
    def target(self) -> AbstractCoordinateDomain:
        """The coordinate domain this map lands in.

        Where a map lands, and not only how many coordinates it lands in:
        the inclusion of a facet into a cell lands in the cell's
        coordinates, where a mesh's own map lands in the ambient ones.
        That is what tells a domain built from a map which coordinates
        name its points.
        """

    @property
    def target_dimension(self) -> int:
        """The number of coordinates this map lands in."""
        return self.target.geometric_dimension

    @abstractmethod
    def value(self, point: AbstractVariable) -> AbstractExpression:
        """Evaluate this map at a point of its source.

        Args:
            point: A point of this map's source, or a variable standing for one

        Returns:
            The coordinates it maps to, of shape ``(target_dimension,)``
        """

    @abstractmethod
    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """Differentiate this map at a point of its source.

        Args:
            point: A point of this map's source, or a variable standing for one

        Returns:
            The derivative of each target coordinate with respect to each
            source coordinate, of shape ``(target_dimension,
            source.geometric_dimension)``
        """

    @property
    def is_affine(self) -> bool:
        """Whether this map is affine.

        Conservative by default (False): a subclass overrides this when
        it knows its map is affine.
        """
        return False

    @property
    def is_identity(self) -> bool:
        """Whether this map leaves the coordinates it is given untouched.

        Conservative by default (False): a subclass overrides this when
        it knows its map is the identity.
        """
        return False

    @abstractmethod
    def __eq__(self, other) -> bool:
        """Check for equality."""

    @abstractmethod
    def __hash__(self) -> int:
        """Hash."""

    @abstractmethod
    def __repr__(self) -> str:
        """Representation."""


@runtime_checkable
class ReadOnADomain(Protocol):
    """An expression that says which domain it is read on.

    Which domain an expression is on is read off the parts of it that
    know: a function knows from its space, a geometric quantity from the
    domain it differentiates, and a restriction from the domain it
    restricts to. :func:`uflx.functions.extract_domain` asks, and a node
    that does not implement this is asked of its successors instead --
    arithmetic and tensors have no domain of their own.

    A node that has already been positioned claims nothing. A function
    pulled back onto a cell no longer says which domain it came from, and
    a quantity told where it is read was put there by whatever knew the
    point. Both answer None, so the question is about an expression as it
    is written rather than one mid-transformation.
    """

    @property
    def reads_on(self) -> AbstractDomain | None:
        """The domain this is read on, or None where it does not say."""


@runtime_checkable
class IncludesASubEntity(Protocol):
    """A map that says which sub-entity of a cell it includes.

    The inclusion an entity supplies is what identifies a kernel along
    with the cell type: one kernel per sub-entity, since each has its own
    map. A consumer that generates code needs the index rather than the
    matrix, so an inclusion that knows which sub-entity it is says so,
    and the index survives into a pulled back integral in the quantities
    that name the domain it parametrizations.

    Not required of an inclusion, so a map that does not know is still a
    parametrization like any other. :class:`uflx.parametrizations.SubEntityInclusion`
    is the affine implementation.
    """

    @property
    def sub_entity(self) -> tuple[int, int]:
        """The dimension of the sub-entity included, and which one it is."""


@runtime_checkable
class ComposesADomain(Protocol):
    """A domain presented as another domain with a map applied to it.

    What a composite is a composite *of*, which is how a restriction asks
    a facet of a mesh which sub-entity of which cell it is, and how the
    sub-entity a parametrization includes is read back off it.
    :class:`uflx.parametrizations.ComposedDomain` is the implementation.

    A label is not one of these. A labelled domain is the domain it
    labels in everything geometric, so whatever asks this question asks
    it of :func:`unlabelled` first rather than of the label.
    """

    @property
    def inner_domain(self) -> AbstractCellularDomain:
        """The domain being mapped on."""

    @property
    def outer_map(self) -> AbstractParametrization:
        """The map applied to that domain's own coordinates."""


@runtime_checkable
class MadeOfPieces(Protocol):
    """A domain integrated over as the sum over its pieces.

    A union is one, each piece a parametrized domain in its own right,
    and so is an interface of unions, each piece a pair of them. An
    integral over such a domain is split into one integral per piece
    before anything is pulled back, since a piece is what names one
    parametrization.
    """

    @property
    def pieces(self) -> tuple[AbstractParametrizedDomain, ...]:
        """The domains this one is integrated over as a sum of."""


@runtime_checkable
class EmbedsAnInterface(Protocol):
    """A domain that is one side of an interface.

    An interface is a domain two others meet along, and a side is that
    domain presented from one of them: the same set, parametrized by the
    one domain's own map. A field on that domain is read on the
    interface through the side, and a quantity of the side -- the normal
    pointing out of it -- is read at the point of the side's own
    parameter region that the interface's point corresponds to, which
    ``transition`` gives. :class:`uflx.interfaces.Embedding` is the
    implementation; this is what the geometry and the integrals ask of
    it without naming it.
    """

    @property
    def interface(self) -> AbstractParametrizedDomain:
        """The interface this is a side of."""

    @property
    def image(self) -> AbstractParametrizedDomain:
        """The interface as this side's domain presents it."""

    def transition(self, variable: AbstractVariable) -> AbstractVariable:
        """Where in this side's parameter region a point of the interface is.

        Args:
            variable: A point of a parameter region of the interface

        Returns:
            The point of this side's parameter region that lands on it,
            or the variable itself where this side's region is the
            interface's own
        """


@runtime_checkable
class StatesAnOutwardNormal(Protocol):
    """A domain that says which way is out of the domain it bounds.

    A face of a box knows from its own description: the held coordinate
    increases out of the box at its upper end. Nothing is computed and no
    orientation is left open, so :class:`uflx.geometry.OutwardNormal`
    takes such a domain at its word rather than transporting a reference
    normal it does not have.
    """

    @property
    def outward_normal(self) -> tuple[float, ...]:
        """The direction out of the bounded domain, as numbers."""


@runtime_checkable
class CarriesADomain(Protocol):
    """A map that carries a whole domain somewhere, rather than a region of coordinates.

    A map out of a coordinate domain is what a parametrization is, and
    most say nothing about where their coordinates came from. One that
    is a domain's own parametrization followed by something, as
    :class:`uflx.interfaces.CoordinateRepresentation` is, starts from that
    domain in all but name, and says so, which is how an interface knows
    a side is presented through its own domain and needs no chart
    transition.
    """

    @property
    def carries(self) -> AbstractParametrizedDomain:
        """The domain this map carries."""


class IdentityParametrization(AbstractParametrization):
    """The map that leaves a coordinate domain's points where they are."""

    def __init__(self, domain: AbstractCoordinateDomain):
        """Initialise.

        Args:
            domain: The coordinate domain this map starts from and lands in
        """
        self._domain = domain

    @property
    def source(self) -> AbstractCoordinateDomain:
        """The coordinate domain this map starts from."""
        return self._domain

    @property
    def target(self) -> AbstractCoordinateDomain:
        """This map lands where it starts."""
        return self._domain

    def value(self, point: AbstractVariable) -> AbstractExpression:
        """Give back the coordinates of the point."""
        return Vector([point.component(i) for i in range(self.target_dimension)])

    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """Differentiating the identity gives the identity."""
        return Identity(self.target_dimension)

    @property
    def is_affine(self) -> bool:
        """The identity is affine."""
        return True

    @property
    def is_identity(self) -> bool:
        """The identity is the identity."""
        return True

    def __repr__(self) -> str:
        """Representation."""
        return f"IdentityParametrization({self._domain!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, IdentityParametrization) and self._domain == other._domain

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.IdentityParametrization", self._domain))


class UnionDomain(AbstractParametrizedDomain):
    """A domain presented as several parametrized pieces.

    Every other parametrized domain here has one parametrization per cell type, keyed by
    the parameter region it starts from. The facets of a triangle are not
    like that: all three are parametrized from the reference interval, so they
    cannot be told apart by their region. A union holds them as pieces
    instead, each a parametrized domain in its own right, and integration over
    the union is the sum of the integrals over the pieces.

    Not an atlas. An atlas's parametrizations cover one set and are related by
    transition maps; these are separate pieces with no claim of a relation
    between them. Whether they meet, and which of them a mesh actually
    has, is external, as the cells of a mesh are.

    Which selection of them is meant is carried by the pieces, each a
    :class:`LabelledDomain` where a consumer has marked one, and not by
    the union: a union of differently labelled pieces is already a
    different domain, and a label on the union alone would be dropped the
    moment an integral over it were split onto one piece.
    """

    def __init__(self, pieces: tuple[AbstractParametrizedDomain, ...]):
        """Initialise.

        Args:
            pieces: The parametrized domains being united, in order

        Raises:
            ValueError: If there are no pieces, or they are not in the same
                coordinates
        """
        if len(pieces) == 0:
            raise ValueError("A union of no pieces is not a domain.")
        dimensions = {piece.geometric_dimension for piece in pieces}
        if len(dimensions) != 1:
            raise ValueError(
                f"A union's pieces must be in the same coordinates, not {sorted(dimensions)}."
            )
        self._pieces = pieces

    @property
    def pieces(self) -> tuple[AbstractParametrizedDomain, ...]:
        """The parametrized domains being united."""
        return self._pieces

    @property
    def geometric_dimension(self) -> int:
        """The coordinates every piece is named in."""
        return self._pieces[0].geometric_dimension

    @property
    def topological_dimension(self) -> int | None:
        """The pieces' common topological dimension, or None when they differ."""
        dimensions = {piece.topological_dimension for piece in self._pieces}
        if len(dimensions) == 1:
            (dimension,) = dimensions
            return dimension
        return None

    @property
    def parametrizations(self) -> tuple[AbstractParametrization, ...]:
        """Every piece's parametrizations, in order.

        Several of them may start from one region, which is the whole
        reason a union exists rather than a domain with a parametrization per cell
        type.
        """
        return tuple(
            parametrization for piece in self._pieces for parametrization in piece.parametrizations
        )

    def __repr__(self) -> str:
        """Representation."""
        pieces = ", ".join(repr(piece) for piece in self._pieces)
        return f"UnionDomain(({pieces}))"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, UnionDomain) and self._pieces == other._pieces

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.UnionDomain", self._pieces))


class AbstractCellularDomain(AbstractParametrizedDomain):
    """Base class for a domain decomposed into cells.

    Being made of cells is what lets a finite element be attached per
    cell, so these are the domains a finite element function space can
    live on. The decomposition itself stays external, like the mesh: such
    a domain says which cell types occur in it, not how many cells there
    are or where they sit.

    Each cell type has a parametrization: a map out of that cell's
    coordinate domain into the coordinates this domain lives in. Being
    made of cells and being the image of a map per cell are one thing
    here, since the cells are what the maps are out of, and a cell's
    coordinate domain is the parameter region that parametrization starts from.
    """

    @property
    @abstractmethod
    def cell_types(self) -> tuple[AbstractEntity, ...]:
        """Get the cell types that occur in this domain."""

    @abstractmethod
    def on_cell_type(self, cell: AbstractEntity) -> AbstractCellularDomain:
        """Get the part of this domain made of cells of one type.

        Integrating over a domain of several cell types means integrating
        over each of these and adding the results, since each cell type
        has its own coordinate domain and its own map out of it.

        Args:
            cell: A cell type of this domain

        Returns:
            A domain whose only cell type is the given one, parametrized
            by the same map for the cell type that is kept

        Raises:
            ValueError: If the cell is not a cell type of this domain
        """

    @property
    def parametrizations(self) -> tuple[AbstractParametrization, ...]:
        """One parametrization per cell type, being that cell type's parametrization."""
        return tuple(self.parametrization(cell) for cell in self.cell_types)

    @abstractmethod
    def parametrization(self, cell: AbstractEntity) -> AbstractParametrization:
        """Get the map out of the given cell type's coordinate domain.

        Args:
            cell: A cell type of this domain

        Returns:
            That cell type's parametrization
        """

    @property
    def has_affine_parametrization(self) -> bool:
        """Whether every cell's parametrization is affine."""
        return all(self.parametrization(c).is_affine for c in self.cell_types)


class EntityDomain(AbstractCoordinateDomain, AbstractCellularDomain):
    """The coordinate realization of a single topological entity.

    An entity domain's points are already coordinate tuples, so its
    parametrization is the identity. The coordinates of those points are
    fixed by whoever defines the elements on the entity, not by UFLx.
    """

    def __init__(self, entity: AbstractEntity):
        """Initialise.

        Args:
            entity: The entity whose coordinates this domain carries
        """
        self._entity = entity

    @property
    def geometric_dimension(self) -> int:
        """The number of coordinates needed to name a point of this domain."""
        return self._entity.topological_dimension

    @property
    def cell_types(self) -> tuple[AbstractEntity, ...]:
        """Get the cell types in this domain, which is the entity itself."""
        return (self._entity,)

    def parametrization(self, cell: AbstractEntity) -> AbstractParametrization:
        """Get the identity map on this domain."""
        if cell != self._entity:
            raise ValueError(f"{cell} is not the entity of this domain.")
        return IdentityParametrization(self)

    def on_cell_type(self, cell: AbstractEntity) -> EntityDomain:
        """An entity domain has one cell type already, which is its entity."""
        if cell != self._entity:
            raise ValueError(f"{cell} is not the entity of this domain.")
        return self

    def __repr__(self) -> str:
        """Representation."""
        return f"EntityDomain({self._entity!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, EntityDomain) and self._entity == other._entity

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.EntityDomain", self._entity))


class LabelledParametrization(AbstractParametrization):
    """A map and the label it is told apart by.

    To a map what :class:`LabelledDomain` is to a domain: everything
    geometric is the map's, and what changes is identity. Two sides of
    one mesh are parametrized by one map, structurally, so a preimage
    through one side's map says which side only if the map does, and
    this is how it does.
    """

    def __init__(self, parametrization: AbstractParametrization, label: Any):
        """Initialise.

        Args:
            parametrization: The map being named
            label: Anything hashable it is told apart by
        """
        self._parametrization = parametrization
        self._label = label

    @property
    def parametrization(self) -> AbstractParametrization:
        """The map being named."""
        return self._parametrization

    @property
    def label(self) -> Any:
        """What it is told apart by."""
        return self._label

    @property
    def source(self) -> AbstractCoordinateDomain:
        """The map's source."""
        return self._parametrization.source

    @property
    def target(self) -> AbstractCoordinateDomain:
        """The map's target."""
        return self._parametrization.target

    def value(self, point: AbstractVariable) -> AbstractExpression:
        """The map's value."""
        return self._parametrization.value(point)

    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """The map's Jacobian."""
        return self._parametrization.jacobian(point)

    @property
    def is_affine(self) -> bool:
        """Whether the map is affine."""
        return self._parametrization.is_affine

    @property
    def is_identity(self) -> bool:
        """Whether the map is the identity."""
        return self._parametrization.is_identity

    def __repr__(self) -> str:
        """Representation."""
        return f"LabelledParametrization({self._parametrization!r}, {self._label!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, LabelledParametrization)
            and self._parametrization == other._parametrization
            and self._label == other._label
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.LabelledParametrization", self._parametrization, self._label))


class LabelledDomain(AbstractCellularDomain):
    """A domain and the label a consumer selects it by.

    Everything geometric is the domain's: the same cells, the same parametrizations
    and the same coordinates, so no quantity, no density and no kernel
    changes. What changes is identity, which is the point: a marked part
    of a domain is a domain, as INTEGRATION.md note 3 says, and two
    selections are two domains rather than one measure carrying a
    subdomain id.

    The label is opaque and hashable. Which entities a consumer has
    marked is external, as a mesh's cells are; this records only that a
    selection was meant and which one.

    It labels a piece rather than a union of them, so that splitting an
    integral keeps it: the pull back leaves quantities naming the piece in
    the integrand, so a kernel says which selection it is for.
    """

    def __init__(self, domain: AbstractCellularDomain, label: Any):
        """Initialise.

        Args:
            domain: The domain being selected
            label: Anything hashable the consumer selects it by
        """
        self._domain = domain
        self._label = label

    @property
    def domain(self) -> AbstractCellularDomain:
        """The domain being selected."""
        return self._domain

    @property
    def label(self) -> Any:
        """What the consumer selects it by."""
        return self._label

    @property
    def geometric_dimension(self) -> int:
        """A label changes no geometry."""
        return self._domain.geometric_dimension

    @property
    def topological_dimension(self) -> int | None:
        """A label changes no topology."""
        return self._domain.topological_dimension

    @property
    def coordinate_domain(self) -> AbstractCoordinateDomain:
        """A label changes no coordinates."""
        return self._domain.coordinate_domain

    @property
    def cell_types(self) -> tuple[AbstractEntity, ...]:
        """The cells are the domain's."""
        return self._domain.cell_types

    def parametrization(self, cell: AbstractEntity) -> AbstractParametrization:
        """The parametrization is the domain's."""
        return self._domain.parametrization(cell)

    def on_cell_type(self, cell: AbstractEntity) -> LabelledDomain:
        """Restricting keeps the label, a cell type being no part of it."""
        return LabelledDomain(self._domain.on_cell_type(cell), self._label)

    def __repr__(self) -> str:
        """Representation."""
        return f"LabelledDomain({self._domain!r}, {self._label!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, LabelledDomain)
            and self._domain == other._domain
            and self._label == other._label
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.LabelledDomain", self._domain, self._label))


def entity_domain(entity: AbstractEntity) -> EntityDomain:
    """Create the coordinate domain of an entity.

    Args:
        entity: The entity whose coordinates the domain carries

    Returns:
        The entity's coordinate domain
    """
    return EntityDomain(entity)


def union_domain(pieces: Sequence[AbstractParametrizedDomain]) -> UnionDomain:
    """Create the domain several parametrized pieces make up.

    Args:
        pieces: The parametrized domains being united, in order

    Returns:
        The union of those pieces
    """
    return UnionDomain(tuple(pieces))


def labelled_domain(domain: AbstractCellularDomain, label: Any) -> LabelledDomain:
    """Create the domain a consumer has selected under a label.

    Args:
        domain: The domain being selected
        label: Anything hashable the consumer selects it by

    Returns:
        The same domain, told apart by that label
    """
    return LabelledDomain(domain, label)


def unlabelled(domain: AbstractDomain) -> AbstractDomain:
    """The domain a selection selects from, or the domain itself.

    A label says which entities a consumer means and nothing about the
    fields on them, so a function on a domain is a function on any
    selection of it. This is how to ask what it is a selection of.

    Args:
        domain: A domain, labelled or not

    Returns:
        The domain with any labels taken off
    """
    while isinstance(domain, LabelledDomain):
        domain = domain.domain
    return domain


def pieces_of(domain: AbstractParametrizedDomain) -> tuple[AbstractParametrizedDomain, ...]:
    """The pieces a domain is integrated over as a sum of, or the domain as its own one.

    A domain not made of pieces is its own one piece, and so is one that
    says so by answering itself, as an interface of one piece a side does.

    Args:
        domain: Any parametrized domain

    Returns:
        Its pieces, or ``(domain,)``
    """
    if isinstance(domain, MadeOfPieces):
        pieces = domain.pieces
        if pieces != (domain,):
            return pieces
    return (domain,)


def parameter_regions(domain: AbstractParametrizedDomain) -> set[AbstractCoordinateDomain]:
    """The regions a domain is parametrized from.

    Args:
        domain: Any parametrized domain

    Returns:
        The sources of its parametrizations
    """
    return {parametrization.source for parametrization in domain.parametrizations}


def parametrization_out_of(
    domain: AbstractParametrizedDomain, region: AbstractCoordinateDomain
) -> AbstractParametrization:
    """The one parametrization of a domain that starts from a region.

    A domain holds its parametrizations rather than keying them by their
    region, since several may share one. Needing exactly one is what
    evaluating a quantity at a point asks for, so it is asked here.

    Args:
        domain: The domain
        region: The region a point lies in

    Returns:
        The parametrization out of it

    Raises:
        ValueError: If the domain has no parametrization out of that
            region, or several, so that a point there names no one map
    """
    found = [p for p in domain.parametrizations if p.source == region]
    if len(found) == 0:
        raise ValueError(f"{domain!r} has no parametrization out of {region!r}.")
    if len(found) > 1:
        raise ValueError(
            f"{domain!r} has {len(found)} parametrizations out of {region!r}, so a point "
            f"there names no one map. Split the integral by piece first."
        )
    return found[0]
