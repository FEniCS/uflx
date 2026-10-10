# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Function spaces.

A function space is a space containing functions defined on a domain.
In most if not all cases, these will be finite dimensional.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence

from uflx.domains import AbstractCellularDomain, AbstractDomain
from uflx.entities import AbstractEntity
from uflx.finite_elements import AbstractMappedFiniteElement


class AbstractFunctionSpace(ABC):
    """Abstract base class for a function space."""

    @property
    @abstractmethod
    def domain(self) -> AbstractDomain:
        """Domain of the function space."""

    @property
    @abstractmethod
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the function space."""

    @property
    @abstractmethod
    def real_valued(self) -> bool:
        """Check if this space is real-valued."""

    @abstractmethod
    def __eq__(self, other) -> bool:
        """Check if this space is equal to another space."""

    @abstractmethod
    def __hash__(self) -> int:
        """Hash."""


class AbstractMappedFunctionSpace(AbstractFunctionSpace):
    """Abstract base class for a function space whose functions are mapped from an entity."""

    @property
    @abstractmethod
    def domain(self) -> AbstractCellularDomain:
        """Domain of the function space.

        An element is attached per cell, so the domain is decomposed into
        cells.
        """

    @property
    @abstractmethod
    def elements(self) -> tuple[AbstractMappedFiniteElement, ...]:
        """Elements in the function space."""

    def restricted_to(self, cell: AbstractEntity) -> AbstractMappedFunctionSpace:
        """Get this space on the part of its domain made of one cell type.

        Every element on that cell is kept, which is what makes this right
        both for a space holding several elements on one cell and for a
        space holding one element per cell type of a mixed domain.

        Args:
            cell: A cell type of this space's domain

        Returns:
            The space with that cell's elements on the restricted domain

        Raises:
            ValueError: If this space has no element on that cell
        """
        elements = tuple(e for e in self.elements if e.cell == cell)
        if len(elements) == 0:
            raise ValueError(f"This space has no element on {cell}.")
        return function_space(self.domain.restricted_to(cell), elements)


class FunctionSpace(AbstractMappedFunctionSpace):
    """Function space."""

    def __init__(
        self,
        domain: AbstractCellularDomain,
        elements: tuple[AbstractMappedFiniteElement, ...],
    ):
        """Initialise."""
        self._domain = domain
        self._elements = elements
        for element in elements:
            if element.cell not in domain.cell_types:
                raise ValueError(
                    f"Element on cell {element.cell} is not defined on a cell of its domain."
                )
        shape = elements[0].ambient_value_shape(domain)
        for element in elements[1:]:
            if element.ambient_value_shape(domain) != shape:
                raise ValueError(
                    "Elements in a functions space must have the same ambient value shape."
                )

    @property
    def domain(self) -> AbstractCellularDomain:
        """Domain of the function space."""
        return self._domain

    @property
    def real_valued(self) -> bool:
        """Check if this space is real-valued."""
        return all(e.real_valued for e in self._elements)

    @property
    def elements(self) -> tuple[AbstractMappedFiniteElement, ...]:
        """Elements in the function space."""
        return self._elements

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the function space."""
        return self.elements[0].ambient_value_shape(self.domain)

    def __repr__(self) -> str:
        """Representation."""
        elements = ", ".join(repr(e) for e in self._elements)
        return f"FunctionSpace({self._domain!r}, {elements})"

    def __eq__(self, other) -> bool:
        """Check if this space is equal to another space.

        Two spaces are equal when they are the same elements on the same
        domain. The element order is part of that, since elements[0] is
        what a pull back reaches for.
        """
        return (
            isinstance(other, FunctionSpace)
            and self._domain == other._domain
            and self._elements == other._elements
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.FunctionSpace", self._domain, self._elements))


def function_space(
    domain: AbstractCellularDomain,
    elements: Sequence[AbstractMappedFiniteElement] | AbstractMappedFiniteElement,
) -> FunctionSpace:
    """Create a function space.

    Args:
        domain: The domain on which the function space is defined.
        elements: The elements in the function space.
    """
    if isinstance(elements, AbstractMappedFiniteElement):
        elements = (elements,)
    return FunctionSpace(domain, tuple(elements))
