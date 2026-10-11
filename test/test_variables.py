"""Test variables."""

from uflx import parametrized_domain
from uflx.functions import Variable


def test_finite_element_variable_round_trip(lagrange_element):
    """Pulling a variable back and pushing it forward gives the original variable."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    variable = Variable(domain)

    (cell,) = domain.cell_types
    entity = variable.to_entity_coordinates(cell)
    assert entity.in_entity_coordinates
    assert not entity.to_domain_coordinates(domain).in_entity_coordinates
    assert entity.to_domain_coordinates(domain) == variable
