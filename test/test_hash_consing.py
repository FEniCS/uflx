"""Test hash-consing of expressions."""

import gc
import time
import weakref
from typing import Any

import numpy as np

from uflx import Coefficient, coordinate_element, function_space
from uflx.algorithms import replace
from uflx.expressions import AbstractExpression, Integer, Rational, RealScalar
from uflx.geometry import SpatialCoordinate
from uflx.graphs import GraphNode, generate_graph


def _build(leaf: AbstractExpression, n: int) -> AbstractExpression:
    """Build e_n = e_{n-1} / e_{n-1} - e_{n-1}, which references e_{n-1} three times."""
    e = leaf
    for _ in range(n):
        e = e / e - e
    return e


def test_structurally_equal_expressions_are_identical():
    """Independently built, structurally equal expressions are the same object."""
    x = SpatialCoordinate(2)
    a = _build(x.component(0), 5)
    b = _build(x.component(0), 5)
    assert a is b
    assert a is not _build(x.component(1), 5)


def test_interning_uses_normalised_init_args():
    """A constructor that normalises its input interns on the normalised value."""
    assert Rational(2, 4) is Rational(1, 2)
    assert RealScalar(1.5) is RealScalar(1.5)
    assert Integer(3) is Integer(3)
    assert Integer(3) == 3


def test_equal_numbers_that_print_differently_are_not_merged():
    """1 and 1.0, and 0.0 and -0.0, compare equal but generate different code."""
    assert RealScalar(-0.0) is not RealScalar(0.0)
    assert str(RealScalar(-0.0).value) == "-0.0"
    assert isinstance(RealScalar(1.0).value, float)
    assert RealScalar(1) is not RealScalar(1.0)


def test_coefficients_with_distinct_labels_stay_distinct(lagrange_element):
    """Coefficients are only merged when their labels match."""
    space = function_space(
        coordinate_element(lagrange_element("triangle", 1, (2,))),
        lagrange_element("triangle", 1),
    )
    assert Coefficient(space) is not Coefficient(space)
    assert Coefficient(space, "f") is Coefficient(space, "f")


def test_replace_returns_interned_expression():
    """replace() rebuilds through the constructors, so its result is interned."""
    x = SpatialCoordinate(2)
    x0, x1 = x.component(0), x.component(1)
    assert replace(_build(x0, 10), {x0: x1}) is _build(x1, 10)


def test_interning_does_not_keep_expressions_alive():
    """The interning table holds expressions weakly."""
    ref = weakref.ref(_build(SpatialCoordinate(3).component(2), 3))
    gc.collect()
    assert ref() is None


class _Unhashable(AbstractExpression):
    """An expression whose init arg cannot be hashed."""

    def __init__(self, data: np.ndarray):
        self.data = data

    @property
    def value_shape(self) -> tuple[int, ...]:
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        return (self.data,)

    def component(self, *indices: int) -> AbstractExpression:
        raise ValueError("Cannot get a component of a scalar expression")


def test_unhashable_init_args_are_not_interned():
    """An expression with an unhashable init arg is constructed but not interned."""
    data = np.array([1.0, 2.0])
    assert _Unhashable(data) is not _Unhashable(data)


def test_shared_subexpressions_scale_linearly():
    """Hash, equality, graph generation and replace are linear in the number of nodes.

    Without hash-consing and a single-expansion graph traversal, every one of these was
    exponential (or quadratic) in the number of levels.
    """
    x = SpatialCoordinate(2)
    x0, x1 = x.component(0), x.component(1)
    n = 1000

    start = time.perf_counter()
    a = _build(x0, n)
    assert a == _build(x0, n)
    assert hash(a) == hash(_build(x0, n))
    graph = generate_graph(a)
    replaced = replace(a, {x0: x1})
    elapsed = time.perf_counter() - start

    # Each level adds a Sum, a Product, a Reciprocal and a Neg; the leaf is a
    # SingleSpatialCoordinate.
    assert len(graph.nodes) == 4 * n + 1
    assert replaced is _build(x1, n)
    assert elapsed < 5.0
