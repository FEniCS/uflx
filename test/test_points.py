"""Test points."""

import pytest

from uflx.expressions import Integer, RealScalar
from uflx.points import RD, Point


@pytest.mark.parametrize("dim", range(5))
def test_rd(dim):
    """Test R^d set of points."""
    points = RD(dim)
    assert points.geometric_dimension == dim


@pytest.mark.parametrize("dim", range(5))
def test_point(dim):
    """Test a point."""
    point = Point([Integer(i) for i in range(dim)])

    assert point.domain.geometric_dimension == dim


def test_points_of_different_dimensions_differ():
    """A point is not equal to a point with more components."""
    a, b = RealScalar(0.5), RealScalar(1.0)
    assert Point([a]) != Point([a, b])
    assert Point([a, b]) != Point([a])
    assert Point([a, b]) == Point([a, b])
    assert hash(Point([a, b])) == hash(Point([a, b]))


def test_reference_and_physical_points_differ():
    """A point on the reference cell is not the physical point with the same coordinates."""
    a, b = RealScalar(0.5), RealScalar(1.0)
    assert Point([a, b], is_reference=True) != Point([a, b])
    assert len({Point([a, b], is_reference=True), Point([a, b])}) == 2
