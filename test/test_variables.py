"""Test variables."""

from uflx import parametrized_domain
from uflx.functions import FiniteElementVariable


def test_finite_element_variable_round_trip(lagrange_element):
    """Pulling a variable back and pushing it forward gives the original variable."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    variable = FiniteElementVariable(domain)

    reference = variable.to_entity_coordinates()
    assert reference.in_entity_coordinates
    assert not reference.to_ambient_coordinates().in_entity_coordinates
    assert reference.to_ambient_coordinates() == variable
