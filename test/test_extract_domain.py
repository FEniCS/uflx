"""Test which domain an expression is read on.

A restriction has to know the domain the field under it lives on, and
pulling a gradient back has to know whose inverse Jacobian to use. Both
ask this, and both used to get an answer read off the *functions* alone:
an expression of pure geometry said nothing and asserted.

What says is what has not yet been positioned. A function says until it
is pulled back onto a cell, a geometric quantity until it is told where
it is read, and a restriction says the domain it restricts to and is not
looked inside, that being what a restriction is for.
"""

import pytest
from conftest import LagrangeElement, Triangle
from test_facet_domains import PlacedTriangle

from uflx import (
    Coefficient,
    TestFunction,
    function_space,
    inner,
    parametrized_domain,
    restrict,
    sub_entities_of,
)
from uflx.domains import entity_domain
from uflx.expressions import RealScalar
from uflx.functions import extract_domain
from uflx.geometry import SpatialCoordinate, UnitNormal
from uflx.points import Point
from uflx.tensors import Identity


@pytest.fixture
def mesh():
    """A mesh of placed triangles, so its facets can be taken."""
    return parametrized_domain(LagrangeElement(PlacedTriangle(), 1, (2,)))


@pytest.fixture
def space(mesh):
    """A P1 space on that mesh."""
    return function_space(mesh, LagrangeElement(PlacedTriangle(), 1))


def test_a_function_says_which_domain_it_is_read_on(mesh, space):
    """From its space, which is where a function's domain has always come from."""
    assert extract_domain(Coefficient(space)) == mesh
    assert extract_domain(TestFunction(space)) == mesh


def test_a_geometric_quantity_says_too(mesh):
    """Which is the half that was missing: an integrand of pure geometry.

    A traction g(x) on a boundary has no function in it until it meets a
    test function, and reading it on the facets is a restriction like any
    other.
    """
    assert extract_domain(SpatialCoordinate(mesh)) == mesh
    assert extract_domain(SpatialCoordinate(mesh).component(1)) == mesh
    assert extract_domain(UnitNormal(sub_entities_of(mesh, 1))) == sub_entities_of(mesh, 1)


def test_a_quantity_told_where_it_is_read_says_nothing(mesh):
    """It was put there by whatever knew the point, so it is not a claim.

    This is what keeps a pulled back integrand readable: the quantities
    in one are on the piece and the cell at once, each carrying its own
    point, and none of them is claiming the integrand's domain.
    """
    point = Point([RealScalar(0.2), RealScalar(0.3)], entity_domain(PlacedTriangle()))

    with pytest.raises(ValueError, match="says which domain"):
        extract_domain(SpatialCoordinate(mesh, point))


def test_a_restriction_says_the_domain_it_restricts_to(mesh, space):
    """And is not looked inside, the field under it being on another domain.

    That is the whole of what a restriction is, so the field's domain
    disagreeing with the restriction's is no disagreement at all.
    """
    facets = sub_entities_of(mesh, 1)

    assert extract_domain(restrict(Coefficient(space), facets)) == facets
    assert extract_domain(restrict(SpatialCoordinate(mesh), facets)) == facets


def test_arithmetic_and_tensors_are_asked_of_their_parts(mesh, space):
    """They have no domain of their own, so the question goes through them."""
    u, v = Coefficient(space), TestFunction(space)

    assert extract_domain(u * v) == mesh
    assert extract_domain(inner(SpatialCoordinate(mesh), SpatialCoordinate(mesh))) == mesh
    assert extract_domain(u * SpatialCoordinate(mesh).component(0)) == mesh


def test_an_expression_that_says_nothing_is_refused():
    """Rather than asserting, which is what it used to do."""
    with pytest.raises(ValueError, match="says which domain"):
        extract_domain(RealScalar(2.0))
    with pytest.raises(ValueError, match="says which domain"):
        extract_domain(Identity(2))


def test_an_expression_on_two_domains_is_refused(mesh, space):
    """A second domain is reached through a map, which is what a restriction is."""
    elsewhere = parametrized_domain(LagrangeElement(Triangle(), 2, (2,)))

    spanning = SpatialCoordinate(mesh).component(0) * SpatialCoordinate(elsewhere).component(0)

    with pytest.raises(ValueError, match="is read on 2 domains"):
        extract_domain(spanning)
