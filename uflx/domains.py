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
from typing import TYPE_CHECKING, Any

from uflx.entities import AbstractEntity
from uflx.tensors import Identity, Vector

if TYPE_CHECKING:
    from uflx.expressions import AbstractExpression
    from uflx.functions import AbstractVariable


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
    @abstractmethod
    def parameter_regions(self) -> tuple[AbstractCoordinateDomain, ...]:
        """The parameter regions this domain is parametrized out of."""

    @abstractmethod
    def parametrization_out_of(self, source: AbstractCoordinateDomain) -> AbstractParametrization:
        """Get the map out of the given parameter region.

        Args:
            source: One of this domain's parameter regions

        Returns:
            The map out of that region

        Raises:
            ValueError: If this domain has no map out of that region
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
    def parameter_regions(self) -> tuple[AbstractCoordinateDomain, ...]:
        """R^d is parametrized by itself."""
        return (self,)

    def parametrization_out_of(self, source: AbstractCoordinateDomain) -> AbstractParametrization:
        """R^d's own parametrization is the identity.

        Raises:
            ValueError: If the region is not this R^d
        """
        if source != self:
            raise ValueError(f"{self!r} is parametrized by itself, not by {source!r}.")
        return IdentityParametrization(self)

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
    def target_dimension(self) -> int:
        """The number of coordinates this map lands in."""

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
    def target_dimension(self) -> int:
        """This map lands where it starts."""
        return self._domain.geometric_dimension

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

    A label rides along, opaque: two unions of the same pieces and
    different labels are different domains. It is what lets a consumer
    name a part of a boundary without the measure acquiring a subdomain
    id, since a subdomain is a domain.
    """

    def __init__(self, pieces: tuple[AbstractParametrizedDomain, ...], label: Any = None):
        """Initialise.

        Args:
            pieces: The parametrized domains being united, in order
            label: Anything hashable the consumer uses to tell this union
                from another of the same pieces, compared by value and
                never interpreted

        Raises:
            ValueError: If there are no pieces, or they are not in the same
                coordinates
        """
        if len(pieces) == 0:
            raise ValueError("A union of no pieces is not a domain.")
        dimensions = {piece.geometric_dimension for piece in pieces}
        if len(dimensions) != 1:
            raise ValueError(
                f"Every piece of a union must be named in the same coordinates, but "
                f"these are in {sorted(dimensions)}."
            )
        self._pieces = pieces
        self._label = label

    @property
    def pieces(self) -> tuple[AbstractParametrizedDomain, ...]:
        """The parametrized domains being united."""
        return self._pieces

    @property
    def label(self) -> Any:
        """What tells this union from another of the same pieces."""
        return self._label

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
    def parameter_regions(self) -> tuple[AbstractCoordinateDomain, ...]:
        """The parameter regions the pieces are parametrized by, each once."""
        return tuple(
            dict.fromkeys(source for piece in self._pieces for source in piece.parameter_regions)
        )

    def parametrization_out_of(self, source: AbstractCoordinateDomain) -> AbstractParametrization:
        """Get the map out of the given parameter region.

        Only where one piece is parametrized by it. Several pieces of one
        region is the whole reason a union exists, and choosing between
        them is not this domain's to do: an integral over it is split into
        one per piece first.

        Args:
            source: One of this domain's parametrization sources

        Returns:
            The map that piece is parametrized by

        Raises:
            ValueError: If no piece is parametrized by that region, or more than
                one is
        """
        parametrized = [piece for piece in self._pieces if source in piece.parameter_regions]
        if len(parametrized) == 0:
            raise ValueError(f"No piece of {self!r} is parametrized by {source!r}.")
        if len(parametrized) > 1:
            raise ValueError(
                f"{len(parametrized)} pieces of {self!r} are parametrized by {source!r}, so there "
                f"is no one map out of it. Split an integral over this domain into one "
                f"per piece first."
            )
        (piece,) = parametrized
        return piece.parametrization_out_of(source)

    def __repr__(self) -> str:
        """Representation."""
        pieces = ", ".join(repr(piece) for piece in self._pieces)
        if self._label is None:
            return f"UnionDomain(({pieces}))"
        return f"UnionDomain(({pieces}), {self._label!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, UnionDomain)
            and self._pieces == other._pieces
            and self._label == other._label
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.UnionDomain", self._pieces, self._label))


class AbstractCellularDomain(AbstractDomain):
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
    def parameter_regions(self) -> tuple[AbstractCoordinateDomain, ...]:
        """The coordinate domains of this domain's cell types."""
        return tuple(EntityDomain(cell) for cell in self.cell_types)

    def parametrization_out_of(self, source: AbstractCoordinateDomain) -> AbstractParametrization:
        """Get the map out of the given cell's coordinate domain.

        Args:
            source: The coordinate domain of one of this domain's cell types

        Returns:
            That cell type's parametrization

        Raises:
            ValueError: If the region is not a cell of this domain
        """
        if not isinstance(source, EntityDomain):
            raise ValueError(
                f"{self!r} is parametrized by its cells' coordinate domains, not by {source!r}."
            )
        (cell,) = source.cell_types
        return self.parametrization(cell)

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


class SubEntityDomain(AbstractCellularDomain):
    """One sub-entity of an entity, as a domain in the entity's coordinates.

    Facet 0 of a triangle is a domain of topological dimension one whose
    points are named by two coordinates, those of the triangle, and whose
    parametrization is the inclusion the entity supplies. It is of codimension one
    in the cell without being a boundary: which facets of a mesh lie on
    one is a question about a cell complex, and this is only the geometry
    of a sub-entity sitting in its cell.

    The dimension is any, not only a facet's. The machinery is the same
    either way, an edge of a tetrahedron being a domain of this kind too,
    and what does not generalise is a claim rather than a construction: a
    normal is a single direction only in codimension one.
    """

    def __init__(self, cell: AbstractEntity, dim: int, index: int):
        """Initialise.

        Args:
            cell: The entity the sub-entity sits in
            dim: The dimension of the sub-entity
            index: Which sub-entity of that dimension

        Raises:
            ValueError: If the cell has no such sub-entity
        """
        if dim < 0 or dim > cell.topological_dimension:
            raise ValueError(f"{cell} has no sub-entities of dimension {dim}.")
        if index < 0 or index >= cell.sub_entity_count(dim):
            raise ValueError(f"{cell} has no sub-entity {index} of dimension {dim}.")
        self._cell = cell
        self._dim = dim
        self._index = index

    @property
    def cell(self) -> AbstractEntity:
        """The entity the sub-entity sits in."""
        return self._cell

    @property
    def dim(self) -> int:
        """The dimension of the sub-entity."""
        return self._dim

    @property
    def index(self) -> int:
        """Which sub-entity of that dimension."""
        return self._index

    @property
    def sub_entity(self) -> AbstractEntity:
        """The sub-entity itself."""
        return self._cell.sub_entities(self._dim)[self._index]

    @property
    def geometric_dimension(self) -> int:
        """A sub-entity is named in the coordinates of the cell it sits in."""
        return self._cell.topological_dimension

    @property
    def topological_dimension(self) -> int:
        """The dimension of the sub-entity."""
        return self._dim

    @property
    def coordinate_domain(self) -> AbstractCoordinateDomain:
        """Those coordinates are the cell's."""
        return EntityDomain(self._cell)

    @property
    def cell_types(self) -> tuple[AbstractEntity, ...]:
        """The one cell type is the sub-entity."""
        return (self.sub_entity,)

    def parametrization(self, cell: AbstractEntity) -> AbstractParametrization:
        """Get the inclusion of the sub-entity into its cell.

        Raises:
            ValueError: If the cell is not this domain's sub-entity
        """
        if cell != self.sub_entity:
            raise ValueError(f"{cell} is not the sub-entity of this domain.")
        return self._cell.sub_entity_inclusion(self._dim, self._index)

    def on_cell_type(self, cell: AbstractEntity) -> SubEntityDomain:
        """A sub-entity domain has one cell type already, which is the sub-entity."""
        if cell != self.sub_entity:
            raise ValueError(f"{cell} is not the sub-entity of this domain.")
        return self

    def __repr__(self) -> str:
        """Representation."""
        return f"SubEntityDomain({self._cell!r}, {self._dim}, {self._index})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, SubEntityDomain)
            and self._cell == other._cell
            and self._dim == other._dim
            and self._index == other._index
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.SubEntityDomain", self._cell, self._dim, self._index))


def entity_domain(entity: AbstractEntity) -> EntityDomain:
    """Create the coordinate domain of an entity.

    Args:
        entity: The entity whose coordinates the domain carries

    Returns:
        The entity's coordinate domain
    """
    return EntityDomain(entity)


def sub_entity_domain(cell: AbstractEntity, dim: int, index: int) -> SubEntityDomain:
    """Create the domain one sub-entity of an entity is.

    Args:
        cell: The entity the sub-entity sits in
        dim: The dimension of the sub-entity
        index: Which sub-entity of that dimension

    Returns:
        That sub-entity, as a domain in the cell's coordinates
    """
    return SubEntityDomain(cell, dim, index)


def union_domain(pieces: Sequence[AbstractParametrizedDomain], label: Any = None) -> UnionDomain:
    """Create the domain several parametrized pieces make up.

    Args:
        pieces: The parametrized domains being united, in order
        label: Anything hashable telling this union from another of the
            same pieces

    Returns:
        The union of those pieces
    """
    return UnionDomain(tuple(pieces), label)
