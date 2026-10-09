# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Domains.

A domain is a subset of R^d over which something can be integrated.
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


class RD(AbstractDomain):
    """R^d."""

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
