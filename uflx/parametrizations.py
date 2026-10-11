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

Maps compose, so a mesh can be post-composed with a map that is not a
finite element one at all. A mesh of intervals in R^1 composed with
``y -> (y, y^2)`` is a parabola represented exactly, rather than
interpolated through the nodes of a higher degree coordinate element.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from uflx.basis_functions import EvaluatedBasisFunction
from uflx.domains import (
    RD,
    AbstractCoordinateDomain,
    AbstractParametrization,
    EntityDomain,
    UnionDomain,
)
from uflx.entities import AbstractEntity
from uflx.expressions import AbstractExpression, MatrixProduct, RealScalar, expression_sum
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

    def __init__(self, domain: AbstractCellularDomain, element: AbstractMappedFiniteElement):
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
    def target(self) -> AbstractCoordinateDomain:
        """The map's values are a point of the ambient coordinates."""
        return RD(self._target_dimension)

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
    def target(self) -> AbstractCoordinateDomain:
        """Composition lands where the second map lands."""
        return self._outer.target

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


class AffineParametrization(AbstractParametrization):
    """A map given in closed form by its matrix and its offset, ``X -> A X + b``.

    Enough to describe the inclusion of a reference sub-entity into its
    cell, or a linear change of the coordinates a domain lives in, with no
    finite element involved. The numbers are the caller's: an inclusion's
    matrix and offset come from the reference coordinates of the cell,
    which belong to whoever defines the elements on it.
    """

    def __init__(
        self,
        source: AbstractCoordinateDomain,
        target: AbstractCoordinateDomain,
        matrix: Sequence[Sequence[float]],
        offset: Sequence[float],
    ):
        """Initialise.

        Args:
            source: The coordinate domain this map starts from
            target: The coordinate domain it lands in
            matrix: The rows of the matrix of the map
            offset: Where the source's origin lands

        Raises:
            ValueError: If the matrix and the offset do not describe a map
                from the source to the target
        """
        self._source = source
        self._target = target
        self._matrix = tuple(tuple(float(a) for a in row) for row in matrix)
        self._offset = tuple(float(b) for b in offset)
        if len({len(row) for row in self._matrix}) > 1:
            raise ValueError("Every row of the matrix must be the same length.")
        if len(self._offset) != len(self._matrix):
            raise ValueError("The offset must have one entry per row of the matrix.")
        if len(self._matrix) != target.geometric_dimension:
            raise ValueError(
                f"A map into {target!r} lands in {target.geometric_dimension} coordinates, "
                f"so its matrix must have that many rows."
            )
        if len(self._matrix[0]) != source.geometric_dimension:
            raise ValueError(
                f"A map out of {source!r} takes {source.geometric_dimension} coordinates, "
                f"so its matrix must have that many columns."
            )

    @property
    def source(self) -> AbstractCoordinateDomain:
        """The coordinate domain this map starts from."""
        return self._source

    @property
    def target(self) -> AbstractCoordinateDomain:
        """The coordinate domain this map lands in."""
        return self._target

    @property
    def matrix(self) -> tuple[tuple[float, ...], ...]:
        """The rows of the matrix of this map."""
        return self._matrix

    @property
    def offset(self) -> tuple[float, ...]:
        """Where the source's origin lands."""
        return self._offset

    def value(self, point: AbstractVariable) -> AbstractExpression:
        """Multiply the matrix by the point and add the offset."""
        return Vector(
            [
                expression_sum(
                    (RealScalar(a) * point.component(j) for j, a in enumerate(row)),
                    default=RealScalar(0.0),
                )
                + RealScalar(b)
                for row, b in zip(self._matrix, self._offset, strict=True)
            ]
        )

    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """An affine map's derivative is its matrix."""
        return Matrix([[RealScalar(a) for a in row] for row in self._matrix])

    @property
    def is_affine(self) -> bool:
        """An affine map is affine."""
        return True

    @property
    def is_identity(self) -> bool:
        """Whether this map leaves the coordinates it is given untouched."""
        identity = tuple(
            tuple(float(i == j) for j in range(len(self._matrix))) for i in range(len(self._matrix))
        )
        return self._matrix == identity and not any(self._offset)

    def __repr__(self) -> str:
        """Representation."""
        return (
            f"AffineParametrization({self._source!r}, {self._target!r}, "
            f"{self._matrix}, {self._offset})"
        )

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, AffineParametrization)
            and other._source == self._source
            and other._target == self._target
            and other._matrix == self._matrix
            and other._offset == self._offset
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(
            (
                "uflx.AffineParametrization",
                self._source,
                self._target,
                self._matrix,
                self._offset,
            )
        )


class ParametrizedDomain(AbstractCellularDomain):
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

    def on_cell_type(self, cell: AbstractEntity) -> ParametrizedDomain:
        """Keep the element of one cell type, which is a domain of that type alone."""
        if cell not in self._elements:
            raise ValueError(f"{cell} is not a cell type of this domain.")
        return ParametrizedDomain((self._elements[cell],))

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


class ComposedDomain(AbstractCellularDomain):
    """A parametrized domain post-composed with a further map.

    The cells and the topology are the inner domain's; the coordinates
    landed in are the outer map's. The outer map need not be described by
    a finite element, so this is how a mesh is carried onto a surface
    given in closed form.
    """

    def __init__(self, domain: AbstractCellularDomain, outer_map: AbstractParametrization):
        """Initialise.

        Args:
            domain: The domain whose cells are being mapped on
            outer_map: The map applied to the domain's own coordinates
        """
        if outer_map.source.geometric_dimension != domain.geometric_dimension:
            raise ValueError(
                f"Cannot map a domain in {domain.geometric_dimension} coordinates "
                f"with a map starting from {outer_map.source!r}."
            )
        self._domain = domain
        self._outer_map = outer_map

    @property
    def inner_domain(self) -> AbstractCellularDomain:
        """The domain being mapped on.

        What a composite is a composite *of*, which is how a restriction
        asks a facet of a mesh which sub-entity of which cell it is.
        """
        return self._domain

    @property
    def map(self) -> AbstractParametrization:
        """The map applied to that domain's own coordinates."""
        return self._map

    @property
    def geometric_dimension(self) -> int:
        """The composite lands wherever the outer map lands."""
        return self._outer_map.target_dimension

    @property
    def topological_dimension(self) -> int | None:
        """Composing with a map does not change the topology."""
        return self._domain.topological_dimension

    @property
    def cell_types(self) -> tuple[AbstractEntity, ...]:
        """The cells are those of the domain being mapped on."""
        return self._domain.cell_types

    def parametrization(self, cell: AbstractEntity) -> AbstractParametrization:
        """Get the inner map for this cell, followed by the outer map."""
        return composed_parametrization(self._domain.parametrization(cell), self._map)

    def on_cell_type(self, cell: AbstractEntity) -> ComposedDomain:
        """Restrict the domain being mapped on, and keep mapping it the same way."""
        return ComposedDomain(self._domain.on_cell_type(cell), self._outer_map)

    def __repr__(self) -> str:
        """Representation."""
        return f"ComposedDomain({self._domain!r}, {self._outer_map!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, ComposedDomain)
            and self._domain == other._domain
            and self._outer_map == other._outer_map
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.ComposedDomain", self._domain, self._outer_map))


def composed_parametrization(
    inner: AbstractParametrization, outer: AbstractParametrization
) -> AbstractParametrization:
    """Create the map that applies one map and then another.

    Composing with an identity gives the other map back rather than a
    composition of it, so that a parametrization built by composing reads, and
    differentiates, as the map it is.

    Only where the identity is the identity of the same coordinate
    domain. An identity on a cell's coordinates composed with a map out
    of R^d is not a unit: it relabels which domain the composite starts
    from, and dropping it would leave a parametrization that no point of the cell
    is in.

    Args:
        inner: The map applied first
        outer: The map applied to the result of the first

    Returns:
        The composition, or whichever map does the work when the other is
        a unit for it
    """
    if inner.is_identity and inner.source == outer.source:
        return outer
    if outer.is_identity and outer.source == inner.target:
        return inner
    return ComposedParametrization(inner, outer)


def sub_entity_domain(cell: AbstractEntity, dim: int, index: int) -> ComposedDomain:
    """Create the domain one sub-entity of an entity is.

    The sub-entity's own coordinate domain carried into the cell's
    coordinates by the inclusion the cell supplies, which is a
    composition like any other. Facet 0 of a triangle is a domain of
    topological dimension one whose points are named by the triangle's
    two coordinates.

    Args:
        cell: The entity the sub-entity sits in
        dim: The dimension of the sub-entity
        index: Which sub-entity of that dimension

    Returns:
        That sub-entity, as a domain in the cell's coordinates

    Raises:
        ValueError: If the cell has no such sub-entity
    """
    if dim < 0 or dim > cell.topological_dimension:
        raise ValueError(f"{cell} has no sub-entities of dimension {dim}.")
    if index < 0 or index >= cell.sub_entity_count(dim):
        raise ValueError(f"{cell} has no sub-entity {index} of dimension {dim}.")
    sub_entity = cell.sub_entities(dim)[index]
    return ComposedDomain(EntityDomain(sub_entity), cell.sub_entity_inclusion(dim, index))


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
    domain: AbstractCellularDomain, outer_map: AbstractParametrization
) -> ComposedDomain:
    """Create a domain by mapping another domain onward.

    Args:
        domain: The domain whose cells are being mapped on
        outer_map: The map applied to the domain's own coordinates

    Returns:
        A domain whose geometry is the domain's followed by the map
    """
    return ComposedDomain(domain, map)


def affine_parametrization(
    source: AbstractCoordinateDomain,
    target: AbstractCoordinateDomain,
    matrix: Sequence[Sequence[float]],
    offset: Sequence[float],
) -> AffineParametrization:
    """Create a map given by its matrix and its offset.

    Args:
        source: The coordinate domain the map starts from
        target: The coordinate domain it lands in
        matrix: The rows of the matrix of the map
        offset: Where the source's origin lands

    Returns:
        That affine map
    """
    return AffineParametrization(source, target, matrix, offset)


def sub_entities_of(domain: AbstractCellularDomain, codim: int, label: Any = None) -> UnionDomain:
    """Create the domain the sub-entities of a given codimension make up.

    A piece per cell type and sub-entity of that type, each the inclusion
    of the sub-entity into its cell followed by the cell's own map. For
    codimension one and a mesh, that is the domain an exterior facet
    integral is over: one measure over another domain rather than a
    measure of its own.

    The codimension is relative to each cell type, so in a domain of
    mixed topological dimension codimension one is the triangles' edges
    and the tetrahedra's faces at once, and the union reports no
    topological dimension of its own.

    Which sub-entities a mesh actually has, and how many cells each is
    shared by, is external, as its cells are. This says only how such a
    domain is parametrized.

    Args:
        domain: The domain whose cells' sub-entities are taken
        codim: Their codimension in each cell
        label: Anything hashable telling this domain from another of the
            same sub-entities

    Returns:
        The union of those sub-entities, mapped on by the domain

    Raises:
        ValueError: If a cell type has no sub-entities of that codimension
    """
    pieces = []
    for cell in domain.cell_types:
        dim = cell.topological_dimension - codim
        if dim < 0:
            raise ValueError(
                f"{cell} has topological dimension {cell.topological_dimension}, so it "
                f"has no sub-entities of codimension {codim}."
            )
        parametrization = domain.parametrization(cell)
        pieces.extend(
            ComposedDomain(sub_entity_domain(cell, dim, index), parametrization)
            for index in range(cell.sub_entity_count(dim))
        )
    return UnionDomain(tuple(pieces), label)


def facets_of(domain: AbstractCellularDomain, label: Any = None) -> UnionDomain:
    """Create the domain a domain's facets make up.

    The codimension one case of :func:`sub_entities_of`, facet being the
    name of a sub-entity of that codimension.

    Args:
        domain: The domain whose cells' facets are taken
        label: Anything hashable telling this domain from another of the
            same facets

    Returns:
        The union of those facets, mapped on by the domain
    """
    return sub_entities_of(domain, 1, label)
