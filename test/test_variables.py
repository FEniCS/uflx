"""Test variables."""

from uflx import coordinate_element
from uflx.functions import FiniteElementVariable


def test_finite_element_variable_round_trip(lagrange_element):
    """Pulling a variable back and pushing it forward gives the original variable."""
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    variable = FiniteElementVariable(domain)

    reference = variable.to_reference()
    assert reference.is_reference
    assert not reference.to_physical().is_reference
    assert reference.to_physical() == variable
