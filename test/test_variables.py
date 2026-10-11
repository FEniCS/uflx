"""Test variables."""

import pytest

from uflx import parametrized_domain
from uflx.algorithms import replace, simplify
from uflx.domains import RD, EntityDomain, entity_domain
from uflx.expressions import AbstractExpression, RealScalar
from uflx.points import Point
from uflx.variables import PointComponent, Variable


def test_finite_element_variable_round_trip(lagrange_element):
    """Pulling a variable back and pushing it forward gives the original variable."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    variable = Variable(domain)

    (cell,) = domain.cell_types
    entity = variable.to_entity_coordinates(cell)
    assert entity.is_on_an_entity
    assert not entity.to_domain_coordinates(domain).is_on_an_entity
    assert entity.to_domain_coordinates(domain) == variable


def test_a_variables_coordinates_are_the_coordinates_of_the_point_it_stands_for(
    lagrange_element,
):
    """A coordinate domain's points are tuples of numbers, so a variable has them.

    Nobody has written the numbers down, so each is a terminal rather
    than a value, and it is the one a parametrization given in closed form bottoms
    out in: such a parametrization does arithmetic on the coordinates it is handed.
    """
    cell = parametrized_domain(lagrange_element("triangle", 1, (2,))).cell_types[0]
    variable = Variable(EntityDomain(cell))

    assert variable.value_shape == (2,)
    for index in range(2):
        coordinate = variable.component(index)
        assert isinstance(coordinate, PointComponent)
        assert coordinate.point == variable
        assert coordinate.component_index == index
        assert coordinate.value_shape == ()


def test_a_point_that_is_named_gives_its_own_coordinate():
    """So the terminal stands for a coordinate nobody has written down yet."""
    point = Point([RealScalar(0.3), RealScalar(0.4)], RD(2))

    assert point.component(1) == RealScalar(0.4)
    assert not isinstance(point.component(1), PointComponent)


def test_a_variable_whose_points_are_not_tuples_has_no_coordinates(lagrange_element):
    """Which is the question the domain answers, and not the dimension.

    A point of a mesh is not a tuple of numbers, the mesh having many
    cells, so a pair of them names a point of a cell and not of the mesh.
    The variable is told what to do about it rather than refusing
    blankly. :mod:`test_facet_domains` does the same for a union, whose
    points name no one parametrization either.
    """
    mesh = parametrized_domain(lagrange_element("triangle", 1, (2,)))

    assert mesh.geometric_dimension == 2
    with pytest.raises(ValueError, match="is not a coordinate tuple"):
        Variable(mesh).component(0)


def test_a_variable_has_no_coordinate_its_domain_does_not_have(lagrange_element):
    """The count is the domain's, as it is for a point that is named."""
    cell = parametrized_domain(lagrange_element("interval", 1, (1,))).cell_types[0]
    variable = Variable(EntityDomain(cell))

    assert variable.component(0).value_shape == ()
    with pytest.raises(ValueError, match="no coordinate 1"):
        variable.component(1)
    with pytest.raises(ValueError, match="no coordinate -1"):
        variable.component(-1)


def test_substituting_a_point_for_a_variable_gives_that_points_coordinates(lagrange_element):
    """Which is what turns a parametrization's arithmetic into numbers.

    A coordinate is a node of the graph, so a quadrature point reaches it
    where it cannot reach the variable held by a geometric quantity.
    """
    cell = parametrized_domain(lagrange_element("interval", 1, (1,))).cell_types[0]
    variable = Variable(EntityDomain(cell))
    point = Point([RealScalar(0.25)], entity_domain(cell))

    doubled = variable.component(0) + variable.component(0)
    substituted = simplify(replace(doubled, {variable: point}))

    assert isinstance(substituted, AbstractExpression)
    assert substituted.as_float() == pytest.approx(0.5)
