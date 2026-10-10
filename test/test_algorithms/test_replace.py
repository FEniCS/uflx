"""Test replace algorithm."""

from uflx import TestFunction, TrialFunction, dx, function_space, parametrized_domain
from uflx.algorithms import replace
from uflx.integrals import Integral


def test_replace(lagrange_element):
    """Test replace algorithm."""
    element = lagrange_element("triangle", 2)
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    form = u * dx(domain)
    assert isinstance(form, Integral)

    replaced_form = replace(form, {u: v})

    assert isinstance(form.integrand, TrialFunction)
    assert isinstance(replaced_form, Integral)
    assert isinstance(replaced_form.integrand, TestFunction)
