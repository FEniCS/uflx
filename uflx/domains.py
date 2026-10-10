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
composing. A mesh is a parametrized domain, and the mesh itself stays
external to UFLx.

In finite element terms, the two are the reference cell and the mesh. For
a surface mesh of triangles in three dimensions::

    mesh = parametrized_domain(P1_vector)              # tdim 2, gdim 3
    X = Point([a, b], entity_domain(triangle))         # tdim == gdim == 2
    x = PushedForwardPoint(X, mesh).expand_geometry()  # a point in RD(3)

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

if TYPE_CHECKING:
    from uflx.expressions import AbstractExpression
    from uflx.points import AbstractPoint


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


class EntityDomain(AbstractCoordinateDomain, AbstractCellularDomain):
    """The coordinate realization of a single topological entity.

    The geometry of an entity domain is the identity, so unlike a
    parametrized domain it has no parametrization. The coordinates of its
    points are fixed by whoever defines the elements on the entity, not
    by UFLx.
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

    def __repr__(self) -> str:
        """Representation."""
        return f"EntityDomain({self._entity!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, EntityDomain) and self._entity == other._entity

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.EntityDomain", self._entity))


class AbstractParametrizedDomain(AbstractCellularDomain):
    """Base class for a domain presented as the image of a map.

    Each cell type has a parametrization: a map out of that cell's
    coordinate domain into the ambient coordinates. How the map is
    described is the implementer's business, so this class asks only that
    it can be evaluated and differentiated at a point. A finite element
    basis summed against a mesh's coordinate dofs is one such
    description, and a closed form expression is another.
    """

    @abstractmethod
    def parametrization_component(
        self,
        cell: AbstractEntity,
        point: AbstractPoint,
        component: int,
        derivative: tuple[int, ...],
    ) -> AbstractExpression:
        """Evaluate one ambient component of a cell's parametrization at a point.

        Args:
            cell: The cell type whose parametrization to evaluate, one of
                this domain's cell types
            point: A point of that cell's coordinate domain
            component: Which ambient coordinate to return, in
                ``range(geometric_dimension)``
            derivative: How many times to differentiate in each of the
                cell's coordinate directions, so a tuple as long as the
                cell's topological dimension

        Returns:
            The component, as an expression
        """

    @property
    def has_affine_parametrization(self) -> bool:
        """Whether every cell's parametrization is affine.

        Conservative by default (False): a subclass overrides this when
        it knows its maps are affine.
        """
        return False


def entity_domain(entity: AbstractEntity) -> EntityDomain:
    """Create the coordinate domain of an entity.

    Args:
        entity: The entity whose coordinates the domain carries

    Returns:
        The entity's coordinate domain
    """
    return EntityDomain(entity)
