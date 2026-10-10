# Copyright (C) 2026 Jack S. Hale
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Parametrizations described by a finite element.

A parametrized domain needs a map out of each cell's coordinate domain,
and one way to describe that map is with a finite element: the element
gives the basis, and a particular cell's coordinate dofs say which
combination of it that cell uses. Evaluating the map is then an
interpolation sum, which is an expression, so this lives above the
function space and basis function layers rather than in `uflx.domains`.
"""

from collections.abc import Sequence

from uflx.basis_functions import EvaluatedBasisFunction
from uflx.domains import AbstractParametrizedDomain
from uflx.entities import AbstractEntity
from uflx.expressions import AbstractExpression, expression_sum
from uflx.finite_elements import AbstractMappedFiniteElement
from uflx.function_spaces import function_space
from uflx.points import AbstractPoint
from uflx.tensors import FlattenedTensorMap


class ParametrizedDomain(AbstractParametrizedDomain):
    """A domain whose parametrization is described by a finite element per cell."""

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

    def parametrization_component(
        self,
        cell: AbstractEntity,
        point: AbstractPoint,
        component: int,
        derivative: tuple[int, ...],
    ) -> AbstractExpression:
        """Evaluate one ambient component of a cell's parametrization at a point.

        The element gives only the basis, so the map is this basis summed
        against the cell's coordinate dofs, which come from the external
        mesh.
        """
        element = self._elements[cell]
        return expression_sum(
            FlattenedTensorMap((i // self._gdim, i % self._gdim), (self._gdim,))
            * EvaluatedBasisFunction(
                function_space(self, element),
                i,
                point,
                derivative=derivative,
                component=component,
            )
            for i in range(element.dim)
        )

    @property
    def has_affine_parametrization(self) -> bool:
        """Whether every cell's parametrization is affine."""
        return all(self.parametrization_element(c).describes_affine_map for c in self.cell_types)

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
