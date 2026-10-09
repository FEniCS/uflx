"""Complex number algorithms."""

from itertools import product
from typing import Protocol, runtime_checkable

from uflx.algorithms import replace
from uflx.expressions import AbstractExpression, ComplexScalar, Conj, RealScalar
from uflx.graphs import GraphNode, as_graph
from uflx.tensors import Tensor


@runtime_checkable
class ComplexValued(Protocol):
    """A complex valued node."""

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""


def take_real_part(
    expression: GraphNode,
) -> GraphNode:
    """Take the real part of all complex values."""
    return replace(
        expression,
        {
            node: node.re
            for node in as_graph(expression)
            if isinstance(node, ComplexValued) and isinstance(node, GraphNode)
        },
    )


def take_imaginary_part(
    expression: GraphNode,
) -> GraphNode:
    """Take the imaginary part of all complex values."""
    return replace(
        expression,
        {
            node: node.im
            for node in as_graph(expression)
            if isinstance(node, ComplexValued) and isinstance(node, GraphNode)
        },
    )


def _is_literal_zero(value: AbstractExpression) -> bool:
    """Check if an expression is a literal zero."""
    if isinstance(value, RealScalar):
        return value.as_float() == 0.0
    if isinstance(value, Tensor):
        return all(
            _is_literal_zero(value.component(*i))
            for i in product(*(range(n) for n in value.value_shape))
        )
    return False


def conj(value: AbstractExpression) -> AbstractExpression:
    """Get the complex conjugate."""
    if isinstance(value, ComplexValued):
        if _is_literal_zero(value.im):
            return value.re
        return value.re - ComplexScalar(RealScalar(0.0), RealScalar(1.0)) * value.im
    else:
        return Conj(value)
