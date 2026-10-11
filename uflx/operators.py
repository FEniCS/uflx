# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Operators."""

from uflx.complex import conj
from uflx.domains import AbstractCellularDomain
from uflx.expressions import AbstractExpression, BinaryOperator, Sqrt, UnaryOperator
from uflx.functions import AbstractFunction, extract_domain
from uflx.geometry import Jacobian, JacobianInverse, SpatialCoordinate
from uflx.graphs import GraphNode
from uflx.maps import PushedForward
from uflx.tensors import Identity, Tensor, Vector, zero


class Inner(BinaryOperator):
    """Inner product operator.

    NOTE: document what happens here with conjugates.
    """

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class Transpose(UnaryOperator):
    """Transpose of a rank-2 (matrix-shaped) expression."""

    def __init__(self, argument: AbstractExpression):
        """Initialise."""
        if len(argument.value_shape) != 2:
            raise ValueError(
                "transpose() is only defined for rank-2 (matrix-shaped) expressions, got "
                f"shape {argument.value_shape}."
            )
        super().__init__(argument)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.argument.value_shape[::-1]

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        i, j = indices
        return self.argument.component(j, i)


class Tr(UnaryOperator):
    """Trace of a square rank-2 (matrix-shaped) expression."""

    def __init__(self, argument: AbstractExpression):
        """Initialise."""
        shape = argument.value_shape
        if len(shape) != 2 or shape[0] != shape[1]:
            raise ValueError(
                "tr() is only defined for square rank-2 (matrix-shaped) expressions, got "
                f"shape {shape}."
            )
        super().__init__(argument)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")

    def as_complex(self) -> complex:
        """Convert to a complex number."""
        return sum(
            self.argument.component(i, i).as_complex() for i in range(self.argument.value_shape[0])
        )

    def as_float(self) -> float:
        """Convert to a floating point number."""
        return sum(
            self.argument.component(i, i).as_float() for i in range(self.argument.value_shape[0])
        )

    def as_int(self) -> int:
        """Convert to an integer."""
        return sum(
            self.argument.component(i, i).as_int() for i in range(self.argument.value_shape[0])
        )


class Grad(UnaryOperator):
    """Gradient operator.

    Of a field, or of the spatial coordinate. Both are terminals that say
    which domain they are read on, and both have an entity gradient: a
    field's is its basis differentiated, and the spatial coordinate's is
    the Jacobian, the map's derivative being the answer rather than
    something to compute from one.

    Not of an expression. The gradient of a product needs a product rule,
    which is the symbolic differentiation LANGUAGE.md section 4 describes
    and is not this.
    """

    def __init__(self, argument: GraphNode):
        """Initialise.

        Args:
            argument: The field or spatial coordinate to differentiate

        Raises:
            ValueError: If the argument is neither, or is a field already
                pulled back onto a cell
        """
        if isinstance(argument, AbstractFunction):
            if argument.is_on_an_entity:
                raise ValueError(
                    f"{argument!r} is already in a cell's coordinates, so its gradient "
                    f"there is an EntityGrad rather than a Grad."
                )
        elif not isinstance(argument, SpatialCoordinate):
            raise ValueError(
                f"Cannot take the gradient of {argument!r}. A field and the spatial "
                f"coordinate have an entity gradient; an expression of them needs a "
                f"rule for differentiating it, which is not built."
            )
        self._physical_argument = argument
        super().__init__(argument)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        gdim = extract_domain(self._physical_argument).geometric_dimension
        return (*self._physical_argument.value_shape, gdim)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError("Cannot get a 'component' of a Grad")

    def pull_back_to_entity(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the entity's coordinates.

        Raises:
            NotImplementedError: If the field was not pulled back to
                something whose entity gradient is known
        """
        domain = extract_domain(self._physical_argument)
        assert isinstance(domain, AbstractCellularDomain)

        if isinstance(self._physical_argument, SpatialCoordinate):
            # The entity gradient of the coordinates a map lands on is that
            # map's derivative, so this is J K: the identity where the map
            # is square, and the tangential projector on a manifold, which
            # is what simplifying the product gives.
            return Jacobian(domain) @ JacobianInverse(domain)

        if self._physical_argument.is_cellwise_constant:
            # The gradient of a cellwise constant is exactly zero.
            return zero(self.value_shape)

        argument = node_map.get(self.argument, self.argument)
        if isinstance(argument, PushedForward):
            # The last index of the entity gradient is the derivative direction.
            return EntityGrad(argument.function) @ JacobianInverse(domain)
        raise NotImplementedError()


class EntityGrad(UnaryOperator):
    """Gradient operator."""

    def __init__(self, argument: GraphNode):
        """Initialise."""
        assert isinstance(argument, AbstractFunction) and argument.is_on_an_entity
        self._entity_argument = argument
        super().__init__(argument)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (*self._entity_argument.value_shape, self._entity_argument.domain_size)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError(
            "Cannot get a 'component' of a EntityGrad. Try calling expand_geometry first"
        )

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        argument = self._entity_argument
        d = argument.domain_size
        if argument.is_cellwise_constant:
            return zero((*argument.value_shape, d))

        value_shape = argument.value_shape
        diffs = [argument.diff(i) for i in range(d)]
        if value_shape == ():
            # A scalar argument's directional derivatives are themselves
            # scalars -- nothing to index into.
            return Vector(diffs)

        # A non-scalar (eg vector-valued) argument's directional
        # derivatives are each still full value_shape-shaped expressions
        # (EvaluatedBasisFunction.diff keeps the same, un-indexed
        # component), so build the (*value_shape, d) result by indexing
        # into each direction's derivative rather than treating the
        # derivatives themselves as the leaves -- a plain
        # Vector(diffs)/Tensor(diffs) would silently misreport its own
        # shape as (d,), since Tensor's shape inference only looks at the
        # nesting of the Python list it is given, not each leaf's own
        # value_shape.
        def build(indices: tuple[int, ...]):
            if len(indices) == len(value_shape):
                return [diffs[i].component(*indices) for i in range(d)]
            return [build((*indices, k)) for k in range(value_shape[len(indices)])]

        return Tensor(build(()))


def grad(a: AbstractExpression) -> AbstractExpression:
    """The gradient of an expression."""
    if isinstance(a, AbstractFunction) and not a.is_on_an_entity and a.is_cellwise_constant:
        # The Grad of a cellwise constant physical function is zero.
        gdim = a.function_space.domain.geometric_dimension
        return zero((*a.value_shape, gdim))
    return Grad(a)


def inner(a: AbstractExpression, b: AbstractExpression) -> AbstractExpression:
    """Inner product."""
    if a.value_shape != b.value_shape:
        raise ValueError("Incompatible value shapes.")

    if a.value_shape == ():
        return a * conj(b)

    return Inner(a, b)


def sqrt(a: AbstractExpression) -> AbstractExpression:
    """The square root of a scalar expression."""
    return Sqrt(a)


def transpose(a: AbstractExpression) -> AbstractExpression:
    """The transpose of a rank-2 (matrix-shaped) expression."""
    return Transpose(a)


def tr(a: AbstractExpression) -> AbstractExpression:
    """The trace of a square rank-2 (matrix-shaped) expression."""
    return Tr(a)


def sym(a: AbstractExpression) -> AbstractExpression:
    """The symmetric part of a square rank-2 (matrix-shaped) expression."""
    return (a + transpose(a)) / 2


def skew(a: AbstractExpression) -> AbstractExpression:
    """The skew-symmetric (antisymmetric) part of a square rank-2 expression."""
    return (a - transpose(a)) / 2


def dev(a: AbstractExpression) -> AbstractExpression:
    """The deviatoric (trace-free) part of a square rank-2 (matrix-shaped) expression."""
    d = a.value_shape[0]
    return a - (tr(a) / d) * Identity(d)
