"""Test points."""

import pytest

from uflx import parametrized_domain
from uflx.algorithms import replace
from uflx.domains import RD, EntityDomain
from uflx.expressions import Integer, RealScalar
from uflx.points import Point


@pytest.mark.parametrize("dim", range(5))
def test_rd(dim):
    """Test R^d set of points."""
    points = RD(dim)
    assert points.geometric_dimension == dim


@pytest.mark.parametrize("dim", range(5))
def test_point(dim):
    """Test a point."""
    p = Point([Integer(i) for i in range(dim)], RD(dim))

    assert p.domain.geometric_dimension == dim
    assert p.domain == RD(dim)
    assert not p.in_entity_coordinates


def test_points_of_different_dimensions_differ():
    """A point is not equal to a point with more components."""
    a, b = RealScalar(0.5), RealScalar(1.0)
    assert Point([a], RD(1)) != Point([a, b], RD(2))
    assert Point([a, b], RD(2)) != Point([a], RD(1))
    assert Point([a, b], RD(2)) == Point([a, b], RD(2))
    assert hash(Point([a, b], RD(2))) == hash(Point([a, b], RD(2)))


def test_entity_and_ambient_points_differ(lagrange_element):
    """The same coordinates in different domains are different points."""
    domain = parametrized_domain(lagrange_element("interval", 1, (1,)))
    a = RealScalar(0.5)

    entity = Point([a], EntityDomain(domain.cell_types[0]))
    ambient = Point([a], RD(1))

    assert entity.in_entity_coordinates
    assert not ambient.in_entity_coordinates
    assert entity != ambient
    assert len({entity, ambient}) == 2


def test_domain_survives_a_rewrite(lagrange_element):
    """A rewrite of a point's components leaves it in the same domain."""
    domain = parametrized_domain(lagrange_element("interval", 1, (1,)))
    a = RealScalar(0.5)
    entity = Point([a], EntityDomain(domain.cell_types[0]))

    rewritten = replace(entity, {a: RealScalar(0.25)})

    assert isinstance(rewritten, Point)
    assert rewritten.domain == entity.domain
    assert rewritten.in_entity_coordinates


def test_point_components_must_match_its_domain(lagrange_element):
    """A point cannot have more components than its domain has dimensions."""
    domain = parametrized_domain(lagrange_element("interval", 1, (1,)))
    a, b = RealScalar(0.5), RealScalar(1.0)

    with pytest.raises(AssertionError):
        Point([a, b], EntityDomain(domain.cell_types[0]))
