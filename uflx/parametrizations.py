# Copyright (C) 2026 Jack S. Hale
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Parametrizations built from a finite element, and compositions of maps.

A parametrized domain needs a map out of each cell's coordinate domain.
One way to describe that map is with a finite element: the element gives
the basis, and a particular cell's coordinate dofs say which combination
of it that cell uses. Evaluating such a map is an interpolation sum, which
is an expression, so this lives above the function space and basis
function layers rather than in `uflx.domains`.

Maps compose, so a mesh can be post-composed with a chart that is not a
finite element map at all. A mesh of intervals in R^1 composed with
``y -> (y, y^2)`` is a parabola represented exactly, rather than
interpolated through the nodes of a higher degree coordinate element.
"""

from collections.abc import Sequence

from uflx.basis_functions import EvaluatedBasisFunction
from uflx.domains import (
    AbstractCoordinateDomain,
    AbstractParametrization,
    AbstractParametrizedDomain,
    EntityDomain,
)
from uflx.entities import AbstractEntity
from uflx.expressions import AbstractExpression, MatrixProduct, expression_sum
from uflx.finite_elements import AbstractMappedFiniteElement
from uflx.function_spaces import function_space
from uflx.functions import AbstractVariable
from uflx.points import Point
from uflx.tensors import FlattenedTensorMap, Matrix, Vector


class FiniteElementParametrization(AbstractParametrization):
    """A map described by a finite element basis and a cell's coordinate dofs.

    The element gives only the basis, so the map is that basis summed
    against the coordinate dofs of a particular cell. Those dofs belong to
    the external mesh, which is why this holds the domain they are indexed
    against as well as the element.
    """

    def __init__(self, domain: AbstractParametrizedDomain, element: AbstractMappedFiniteElement):
        """Initialise.

        Args:
            domain: The domain whose coordinate dofs the basis is summed against
            element: The element giving the basis of the map
        """
        self._domain = domain
        self._element = element
        (self._target_dimension,) = element.entity_value_shape

    @property
    def source(self) -> AbstractCoordinateDomain:
        """The coordinate domain of the element's cell."""
        return EntityDomain(self._element.cell)

    @property
    def target_dimension(self) -> int:
        """The map's values are a point of the ambient coordinates."""
        return self._target_dimension

    @property
    def element(self) -> AbstractMappedFiniteElement:
        """The element giving the basis of this map."""
        return self._element

    def _dof_sum(
        self, point: AbstractVariable, component: int, derivative: tuple[int, ...]
    ) -> AbstractExpression:
        """Sum the basis against the coordinate dofs, for one component."""
        dim = self._target_dimension
        return expression_sum(
            FlattenedTensorMap((i // dim, i % dim), (dim,))
            * EvaluatedBasisFunction(
                function_space(self._domain, self._element),
                i,
                point,
                derivative=derivative,
                component=component,
            )
            for i in range(self._element.dim)
        )

    def value(self, point: AbstractVariable) -> AbstractExpression:
        """Interpolate the coordinate dofs at the point."""
        no_derivative = (0,) * self.source.geometric_dimension
        return Vector(
            [self._dof_sum(point, j, no_derivative) for j in range(self._target_dimension)]
        )

    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """Interpolate the coordinate dofs against the basis's derivatives."""
        tdim = self.source.geometric_dimension
        return Matrix(
            [
                [
                    self._dof_sum(point, row, tuple(1 if d == col else 0 for d in range(tdim)))
                    for col in range(tdim)
                ]
                for row in range(self._target_dimension)
            ]
        )

    @property
    def is_affine(self) -> bool:
        """The element knows whether the map it describes is affine."""
        return self._element.describes_affine_map

    def __repr__(self) -> str:
        """Representation."""
        return f"FiniteElementParametrization({self._domain!r}, {self._element!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, FiniteElementParametrization)
            and self._domain == other._domain
            and self._element == other._element
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.FiniteElementParametrization", self._domain, self._element))


class ComposedParametrization(AbstractParametrization):
    """Two maps applied in turn, the second starting where the first lands."""

    def __init__(self, inner: AbstractParametrization, outer: AbstractParametrization):
        """Initialise.

        Args:
            inner: The map applied first
            outer: The map applied to the result of the first
        """
        if outer.source.geometric_dimension != inner.target_dimension:
            raise ValueError(
                f"Cannot compose a map landing in {inner.target_dimension} coordinates "
                f"with one starting from {outer.source!r}."
            )
        self._inner = inner
        self._outer = outer

    @property
    def source(self) -> AbstractCoordinateDomain:
        """Composition starts where the first map starts."""
        return self._inner.source

    @property
    def target_dimension(self) -> int:
        """Composition lands where the second map lands."""
        return self._outer.target_dimension

    def _intermediate_point(self, point: AbstractVariable) -> Point:
        """The point the first map gives the second."""
        value = self._inner.value(point)
        return Point(
            [value.component(k) for k in range(self._inner.target_dimension)],
            self._outer.source,
        )

    def value(self, point: AbstractVariable) -> AbstractExpression:
        """Apply the second map to the first map's value."""
        return self._outer.value(self._intermediate_point(point))

    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """The chain rule."""
        return MatrixProduct(
            [
                self._outer.jacobian(self._intermediate_point(point)),
                self._inner.jacobian(point),
            ]
        )

    @property
    def is_affine(self) -> bool:
        """A composition of affine maps is affine."""
        return self._inner.is_affine and self._outer.is_affine

    @property
    def is_identity(self) -> bool:
        """A composition of identities is the identity."""
        return self._inner.is_identity and self._outer.is_identity

    def __repr__(self) -> str:
        """Representation."""
        return f"ComposedParametrization({self._inner!r}, {self._outer!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, ComposedParametrization)
            and self._inner == other._inner
            and self._outer == other._outer
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.ComposedParametrization", self._inner, self._outer))


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

    def parametrization(self, cell: AbstractEntity) -> FiniteElementParametrization:
        """Get the map out of the given cell type's coordinate domain."""
        return FiniteElementParametrization(self, self._elements[cell])

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


class ComposedDomain(AbstractParametrizedDomain):
    """A parametrized domain post-composed with a further map.

    The cells and the topology are the inner domain's; the coordinates
    landed in are the outer map's. The outer map need not be described by
    a finite element, so this is how a mesh is carried onto a surface
    given in closed form.
    """

    def __init__(self, domain: AbstractParametrizedDomain, map: AbstractParametrization):
        """Initialise.

        Args:
            domain: The domain whose cells are being mapped on
            map: The map applied to the domain's own coordinates
        """
        if map.source.geometric_dimension != domain.geometric_dimension:
            raise ValueError(
                f"Cannot map a domain in {domain.geometric_dimension} coordinates "
                f"with a map starting from {map.source!r}."
            )
        self._domain = domain
        self._map = map

    @property
    def geometric_dimension(self) -> int:
        """The composite lands wherever the outer map lands."""
        return self._map.target_dimension

    @property
    def topological_dimension(self) -> int | None:
        """Composing with a map does not change the topology."""
        return self._domain.topological_dimension

    @property
    def cell_types(self) -> tuple[AbstractEntity, ...]:
        """The cells are those of the domain being mapped on."""
        return self._domain.cell_types

    def parametrization(self, cell: AbstractEntity) -> ComposedParametrization:
        """Get the inner map for this cell, followed by the outer map."""
        return ComposedParametrization(self._domain.parametrization(cell), self._map)

    def __repr__(self) -> str:
        """Representation."""
        return f"ComposedDomain({self._domain!r}, {self._map!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, ComposedDomain)
            and self._domain == other._domain
            and self._map == other._map
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.ComposedDomain", self._domain, self._map))


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


def composed_domain(
    domain: AbstractParametrizedDomain, map: AbstractParametrization
) -> ComposedDomain:
    """Create a domain by mapping another domain onward.

    Args:
        domain: The domain whose cells are being mapped on
        map: The map applied to the domain's own coordinates

    Returns:
        A domain whose geometry is the domain's followed by the map
    """
    return ComposedDomain(domain, map)
