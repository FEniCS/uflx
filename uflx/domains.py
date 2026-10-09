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

For example, take a surface mesh of triangles in three dimensions, whose
geometry is described by a degree 1 Lagrange element with three
components::

    domain = parametrized_domain(P1_vector)   # tdim 2, gdim 3

That domain is parametrized. The tuple (0.25, 0.25) does not name one of
its points: the mesh has many triangles, so there is no one point those
two numbers refer to, and a point of the surface needs three ambient
coordinates anyway. What (0.25, 0.25) does name is a point of the
triangle's own coordinate domain::

    a = b = RealScalar(0.25)
    X = point([a, b], entity_domain(triangle))   # tdim == gdim == 2

A point of the mesh is that tuple carried through the parametrization of
one cell, which is what PushedForwardPoint does: given the coordinate
dofs of a cell, it evaluates the geometry element's basis at X and sums,
giving three explicit ambient coordinates in R^3. So the two kinds of
domain are what a finite element method already distinguishes as the
reference cell and the mesh, and the distinction here is about which of
them a tuple of numbers can name.

The dimensions follow from this. A coordinate domain's geometry is the
identity, so its topological and geometric dimensions always agree --
the triangle's coordinate domain has both equal to 2. A parametrized
domain is where they come apart, with tdim 2 and gdim 3 for the surface
mesh above.

There is no assumption that a domain only contains cells of a single type:
one could contain (eg) a mixture of triangles and quadrilaterals, or even
a mixture of (eq) tetrahedra and intervals.
"""

from abc import ABC, abstractmethod
from collections.abc import Sequence

from uflx.entities import AbstractEntity
from uflx.finite_elements import AbstractMappedFiniteElement


class AbstractDomain(ABC):
    """Base class for a domain."""

    @property
    @abstractmethod
    def geometric_dimension(self) -> int:
        """The dimension of the space this domain is embedded in."""

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
    meaningful. Its geometry is the identity, which means its geometric
    and topological dimensions always agree.
    """


class AbstractFiniteElementDomain(AbstractDomain):
    """Base class for a domain of a finite element function."""

    @property
    @abstractmethod
    def cells(self) -> tuple[AbstractEntity, ...]:
        """Get the cell types in the finite element mesh."""


class AbstractParametrizedDomain(AbstractFiniteElementDomain):
    """Base class for a coordinate element.

    In a coordinate element, the geometry of the domain is defined using a
    finite element.
    """

    @abstractmethod
    def parametrization(self, cell: AbstractEntity) -> AbstractMappedFiniteElement:
        """Get the element describing the geometry of the given cell type."""

    @property
    def is_affine_map(self) -> bool:
        """Is the parametrization of this domain affine?"""
        return all(
            c.is_simplex and self.parametrization(c).lagrange_superdegree == 1 for c in self.cells
        )


class ParametrizedDomain(AbstractParametrizedDomain):
    """A coordinate element."""

    def __init__(self, elements: tuple[AbstractMappedFiniteElement, ...]):
        """Initialise."""
        self._elements = {e.cell: e for e in elements}
        (self._gdim,) = elements[0].entity_value_shape
        for e in elements[1:]:
            assert e.entity_value_shape == (self._gdim,)

    @property
    def geometric_dimension(self) -> int:
        """Dimension of the space this domain is embedded in."""
        return self._gdim

    @property
    def cells(self) -> tuple[AbstractEntity, ...]:
        """Get the cells in the domain."""
        return tuple(self._elements.keys())

    def parametrization(self, cell: AbstractEntity) -> AbstractMappedFiniteElement:
        """Get the element describing the geometry of the given cell type."""
        return self._elements[cell]

    @property
    def topological_dimension(self) -> int | None:
        """The topological dimension of the domain.

        This returns None iff the domain contains entities of a mixture
        of topological dimensions.
        """
        dims = {c.topological_dimension for c in self.cells}
        if len(dims) == 1:
            (dim,) = dims
            return dim
        else:
            return None

    def __repr__(self) -> str:
        """Representation."""
        elements = ", ".join(repr(e) for e in self._elements.values())
        return f"ParametrizedDomain({elements})"

    def __eq__(self, other) -> bool:
        """Check for equality.

        Two parametrized domains are equal when they have the same
        geometric description. This says nothing about the meshes a
        consumer may attach to them, which UFLx never sees.
        """
        return isinstance(other, ParametrizedDomain) and self._elements == other._elements

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.ParametrizedDomain", *sorted(self._elements.items(), key=repr)))


class EntityDomain(AbstractCoordinateDomain, AbstractFiniteElementDomain):
    """The coordinate realization of a single topological entity.

    The geometry of an entity domain is the identity, so unlike a
    :class:`ParametrizedDomain` it has no parametrization. The
    coordinates of its points are fixed by whoever defines the elements
    on the entity, not by UFLx.
    """

    def __init__(self, entity: AbstractEntity):
        """Initialise.

        Args:
            entity: The entity whose coordinates this domain carries
        """
        self._entity = entity

    @property
    def entity(self) -> AbstractEntity:
        """The entity whose coordinates this domain carries."""
        return self._entity

    @property
    def geometric_dimension(self) -> int:
        """The dimension of the space this domain is embedded in."""
        return self._entity.topological_dimension

    @property
    def topological_dimension(self) -> int | None:
        """The topological dimension of the domain."""
        return self._entity.topological_dimension

    @property
    def cells(self) -> tuple[AbstractEntity, ...]:
        """Get the cell types in the finite element mesh."""
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


class RD(AbstractCoordinateDomain):
    """R^d, the ambient coordinate domain."""

    def __init__(self, dim: int):
        """Initialise."""
        self._dim = dim

    @property
    def geometric_dimension(self) -> int:
        """The dimension of the space this domain is embedded in."""
        return self._dim

    @property
    def topological_dimension(self) -> int | None:
        """The topological dimension of the domain.

        This returns None iff the domain contains entities of a mixture
        of topological dimensions.
        """
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


def entity_domain(entity: AbstractEntity) -> EntityDomain:
    """Create the coordinate domain of an entity.

    Args:
        entity: The entity whose coordinates the domain carries

    Returns:
        The entity's coordinate domain
    """
    return EntityDomain(entity)


def parametrized_domain(
    elements: Sequence[AbstractMappedFiniteElement] | AbstractMappedFiniteElement,
):
    """Create a domain.

    Args:
        elements: The finite element(s) used to define the geometry of the cells in this domain
    """
    if isinstance(elements, AbstractMappedFiniteElement):
        elements = (elements,)
    assert len(elements[0].entity_value_shape) == 1
    (gdim,) = elements[0].entity_value_shape
    for e in elements:
        assert e.entity_value_shape == (gdim,)

    return ParametrizedDomain(tuple(elements))
