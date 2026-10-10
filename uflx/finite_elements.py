# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Finite element.

A finite element is an object that is used to define basis functions on a single mesh entity.
The entity on which the element is defined is called the cell.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from math import prod
from typing import Protocol, runtime_checkable

from uflx.entities import AbstractEntity
from uflx.maps import AbstractValueMap


@runtime_checkable
class HasAGeometricDimension(Protocol):
    """Anything that says how many coordinates name a point of it.

    All a value map needs of a domain: the shape a value takes once
    mapped depends on how many coordinates it is mapped into, and on
    nothing else about where it is mapped to. Asking for the capability
    rather than importing the domain leaves this module where the
    layering puts it, below `domains`, and still leaves no room to pass a
    topological dimension, which is what the number used to be before it
    was asked for this way.
    """

    @property
    def geometric_dimension(self) -> int:
        """The number of coordinates needed to name a point."""


class AbstractFiniteElement(ABC):
    """Abstract base class for a finite element.

    To make your element library compatible with UFL, you should make a
    subclass of AbstractFiniteElement and provide implementations of all
    the abstract methods and properties. All methods and properties that
    are not marked as abstract are implemented here and should not need
    to be overwritten in your subclass.
    """

    @abstractmethod
    def __eq__(self, other) -> bool:
        """Check if this element is equal to another element."""

    @property
    @abstractmethod
    def cell(self) -> AbstractEntity:
        """Return the cell that this element is defined on."""

    @property
    @abstractmethod
    def real_valued(self) -> bool:
        """Check if this element is real-valued."""

    @abstractmethod
    def mapped_value_shape(self, domain: HasAGeometricDimension) -> tuple[int, ...]:
        """Return the shape of the value space in a domain's own coordinates.

        The shape a value has once the value map has been applied, which
        is why it is named for the map: those coordinates are the ambient
        ones only when the domain says so, as `coordinate_domain` has it.

        Args:
            domain: The domain whose coordinates the values are expressed
                in. Asking it for its dimension rather than being handed
                one leaves no room to pass a topological dimension.
        """

    def mapped_value_size(self, domain: HasAGeometricDimension) -> int:
        """Return the value size of the value space in a domain's own coordinates."""
        return prod(self.mapped_value_shape(domain))

    @property
    @abstractmethod
    def lagrange_superdegree(self) -> int | None:
        """Degree of the minimum degree Lagrange space that spans this element.

        This returns the degree of the lowest degree Lagrange space such
        that the polynomial space of the Lagrange space is a superspace
        of this element's polynomial space. If this element contains
        basis functions that are not in any Lagrange space, this
        function should return None.

        Note that on a simplex cells, the polynomial space of Lagrange
        space is a complete polynomial space, but on other cells this is
        not true. For example, on quadrilateral cells, the degree 1
        Lagrange space includes the degree 2 polynomial xy.
        """

    @property
    @abstractmethod
    def dim(self) -> int:
        """The dimension of the finite element, ie the number of basis functions."""

    @abstractmethod
    def __hash__(self):
        """Hash."""

    @abstractmethod
    def __repr__(self) -> str:
        """Representation."""


class AbstractMappedFiniteElement(AbstractFiniteElement):
    """Abstract base class for a finite element whose values are mapped from an entity.

    To make your element library compatible with UFL, you should make a
    subclass of AbstractFiniteElement and provide implementations of all
    the abstract methods and properties. All methods and properties that
    are not marked as abstract are implemented here and should not need
    to be overwritten in your subclass.
    """

    @property
    @abstractmethod
    def entity_value_shape(self) -> tuple[int, ...]:
        """Return the shape of the value space in the entity's coordinates."""

    @property
    def entity_value_size(self) -> int:
        """Return the value size of the value space in the entity's coordinates."""
        return prod(self.entity_value_shape)

    @property
    @abstractmethod
    def value_map(self) -> AbstractValueMap:
        """Get the push forward and pull back map."""

    @property
    def describes_affine_map(self) -> bool:
        """Whether the map this element describes is affine.

        A degree 1 element on a simplex gives an affine map. On a
        tensor-product cell even a degree 1 element's map is multilinear,
        not affine. A library that knows more about its own elements may
        override this.
        """
        return self.cell.is_simplex and self.lagrange_superdegree == 1

    def mapped_value_shape(self, domain: HasAGeometricDimension) -> tuple[int, ...]:
        """Return the shape of the value space in a domain's own coordinates."""
        return self.value_map.mapped_value_shape(
            self.entity_value_shape, domain.geometric_dimension
        )
