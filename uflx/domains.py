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

from abc import ABC, abstractmethod
from collections.abc import Sequence

from uflx.entities import AbstractEntity
from uflx.finite_elements import AbstractMappedFiniteElement


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

    The map out of a cell's coordinate domain into the ambient
    coordinates is its parametrization, and a finite element per cell
    type gives that map's basis. The element alone is not the map:
    evaluating it on a cell also needs that cell's coordinate dofs, which
    come from the external mesh.
    """

    @abstractmethod
    def parametrization_element(self, cell: AbstractEntity) -> AbstractMappedFiniteElement:
        """Get the element giving the basis of the given cell type's parametrization.

        The parametrization itself cannot be returned: it is this basis
        summed against a particular cell's coordinate dofs, and those
        come from the external mesh.
        """

    @property
    def has_affine_parametrization(self) -> bool:
        """Is the parametrization of this domain affine?"""
        return all(self.parametrization_element(c).describes_affine_map for c in self.cell_types)


class ParametrizedDomain(AbstractParametrizedDomain):
    """A domain whose geometry is described by a finite element per cell."""

    def __init__(self, elements: tuple[AbstractMappedFiniteElement, ...]):
        """Initialise."""
        self._elements = {e.cell: e for e in elements}
        for e in elements:
            if not e.value_map.is_identity:
                raise ValueError(
                    "A parametrization's values are the ambient coordinates, "
                    "so its element must be identity mapped."
                )
        # Hence a parametrization's value shape is (gdim,), and reading it
        # in the entity's coordinates is the same as reading it in the
        # ambient ones -- which is just as well, since the ambient shape
        # would need the gdim being computed here.
        shapes = {e.entity_value_shape for e in elements}
        if len(shapes) != 1:
            raise ValueError(
                f"Every parametrization of a domain must have the same value shape, got {shapes}."
            )
        (shape,) = shapes
        if len(shape) != 1:
            raise ValueError(
                f"A parametrization's values are a point of R^gdim, so it must be vector "
                f"valued, but its value shape is {shape}."
            )
        (self._gdim,) = shape

    @property
    def geometric_dimension(self) -> int:
        """The number of coordinates needed to name a point of this domain."""
        return self._gdim

    @property
    def topological_dimension(self) -> int | None:
        """The topological dimension of the domain.

        This returns None iff the domain contains entities of a mixture
        of topological dimensions.
        """
        dims = {c.topological_dimension for c in self.cell_types}
        if len(dims) == 1:
            (dim,) = dims
            return dim
        else:
            return None

    @property
    def cell_types(self) -> tuple[AbstractEntity, ...]:
        """Get the cell types in this domain."""
        return tuple(self._elements.keys())

    def parametrization_element(self, cell: AbstractEntity) -> AbstractMappedFiniteElement:
        """Get the element giving the basis of the given cell type's parametrization."""
        return self._elements[cell]

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
        return hash(("uflx.ParametrizedDomain", frozenset(self._elements.items())))


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
) -> ParametrizedDomain:
    """Create a parametrized domain.

    Args:
        elements: The finite element(s) used to define the geometry of each cell type

    Returns:
        A domain whose geometry is described by the given element(s)
    """
    if isinstance(elements, AbstractMappedFiniteElement):
        elements = (elements,)
    return ParametrizedDomain(tuple(elements))
