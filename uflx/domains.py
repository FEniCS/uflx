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
    x = PushedForwardPoint(X, mesh).expand_geometry()  # a point of RD(3)

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
from typing import TYPE_CHECKING

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


class RD(AbstractCoordinateDomain):
    """R^d, the ambient coordinate domain."""

    def __init__(self, dim: int):
        """Initialise."""
        self._dim = dim

    @property
    def geometric_dimension(self) -> int:
        """The number of coordinates needed to name a point of this domain."""
        return self._dim

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


class AbstractCellularDomain(AbstractDomain):
    """Base class for a domain decomposed into cells.

    Being made of cells is what lets a finite element be attached per
    cell, so these are the domains a finite element function space can
    live on. The decomposition itself stays external, like the mesh: such
    a domain says which cell types occur in it, not how many cells there
    are or where they sit.
    """

    @property
    @abstractmethod
    def cell_types(self) -> tuple[AbstractEntity, ...]:
        """Get the cell types that occur in this domain."""


class AbstractParametrizedDomain(AbstractCellularDomain):
    """Base class for a domain presented as the image of a map.

    Each cell type has a parametrization: a map out of that cell's
    coordinate domain into the coordinates this domain lives in.
    """

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


class EntityDomain(AbstractCoordinateDomain, AbstractParametrizedDomain):
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

    def __repr__(self) -> str:
        """Representation."""
        return f"EntityDomain({self._entity!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, EntityDomain) and self._entity == other._entity

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.EntityDomain", self._entity))


def entity_domain(entity: AbstractEntity) -> EntityDomain:
    """Create the coordinate domain of an entity.

    Args:
        entity: The entity whose coordinates the domain carries

    Returns:
        The entity's coordinate domain
    """
    return EntityDomain(entity)
