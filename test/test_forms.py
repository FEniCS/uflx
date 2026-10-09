"""Test forms."""

from uflx import TestFunction, TrialFunction, dx, function_space, inner, parametrized_domain
from uflx.integrals import Integral


def test_simple_form(lagrange_element):
    """Test a simple form."""
    element = lagrange_element("triangle", 2)
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx

    print(form)
    assert isinstance(form, Integral)
    form.graph.print()
