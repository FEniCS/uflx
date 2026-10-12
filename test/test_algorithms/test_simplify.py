"""Test forms."""

import os
import subprocess
import sys

import pytest

from uflx import (
    Coefficient,
    TestFunction,
    TrialFunction,
    coordinate_element,
    dx,
    function_space,
    inner,
)
from uflx.algorithms import simplify
from uflx.algorithms.simplify import digest
from uflx.expressions import MatrixProduct, Product
from uflx.geometry import Jacobian, JacobianInverse, JacobianInverseTranspose, JacobianTranspose
from uflx.integrals import Integral
from uflx.operators import Inner


def test_add_and_subtract_integer(lagrange_element):
    """Test that adding 2 and -2 are successfully cancelled."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)

    expression = u + 2 - 2
    simpler_expression = simplify(expression)

    assert not isinstance(expression, TrialFunction)
    assert isinstance(simpler_expression, TrialFunction)


def test_add_and_subtract_more_integers(lagrange_element):
    """Test that adding and subtracting integers are successfully cancelled."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)

    expression = u + 2 + 6 - 1 - 3 - 4
    simpler_expression = simplify(expression)

    assert not isinstance(expression, TrialFunction)
    assert isinstance(simpler_expression, TrialFunction)


def test_multiply_and_divide_integer(lagrange_element):
    """Test that multiplication then division by 2 are successfully cancelled."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)

    expression = u * 2 / 2
    simpler_expression = simplify(expression)

    assert not isinstance(expression, TrialFunction)
    assert isinstance(simpler_expression, TrialFunction)


def test_multiply_and_divide_more_integers(lagrange_element):
    """Test that multiplication then division by integers are successfully cancelled."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)

    expression = u * 4 / 6 * 3 / 2
    simpler_expression = simplify(expression)

    assert not isinstance(expression, TrialFunction)
    assert isinstance(simpler_expression, TrialFunction)


def test_add_and_subtract_function(lagrange_element):
    """Test that Function and -Function are successfully cancelled."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)

    f = Coefficient(space)

    expression = u - f + f
    simpler_expression = simplify(expression)

    assert not isinstance(expression, TrialFunction)
    assert isinstance(simpler_expression, TrialFunction)


def test_add_real_zero(lagrange_element):
    """Test that adding 0.0 is removed."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)

    assert isinstance(simplify(u + 0.0), TrialFunction)


def test_multiply_and_divide_integer_form(lagrange_element):
    """Test that 2 and 1/2 are successfully cancelled."""
    pytest.xfail()

    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    form = inner(2 * u, v / 2) * dx
    simpler_form = simplify(form)

    assert isinstance(form, Integral)
    assert isinstance(simpler_form, Integral)

    assert isinstance(form.integrand, Product)
    assert len(form.integrand._items) > 2

    assert isinstance(simpler_form.integrand, Product)
    assert isinstance(simpler_form.integrand._items[0], TrialFunction)
    assert isinstance(simpler_form.integrand._items[1], TestFunction)


def test_multiply_and_divide_function_form(lagrange_element):
    """Test that Function and 1/Function are successfully cancelled."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    f = Coefficient(space)

    form = (u * f) * (v / f) * dx
    simpler_form = simplify(form)

    assert isinstance(form, Integral)
    assert isinstance(simpler_form, Integral)

    assert isinstance(form.integrand, Product)
    assert len(form.integrand._items) > 2

    assert isinstance(simpler_form.integrand, Product)
    if isinstance(simpler_form.integrand._items[0], TrialFunction):
        assert isinstance(simpler_form.integrand._items[1], TestFunction)
    else:
        assert isinstance(simpler_form.integrand._items[0], TestFunction)
        assert isinstance(simpler_form.integrand._items[1], TrialFunction)


@pytest.mark.parametrize("v_first", [True, False])
@pytest.mark.parametrize("inv_first", [True, False])
@pytest.mark.parametrize("transpose", [True, False])
def test_jacobian_and_inverse_matvec(lagrange_element, v_first, inv_first, transpose):
    """Test that Jacobian and inverse Jacobian are successfully cancelled."""
    element = lagrange_element("triangle", 2, (2,))
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    v = TestFunction(space)

    if transpose:
        first = JacobianTranspose(domain)
        second = JacobianInverseTranspose(domain)
    else:
        first = Jacobian(domain)
        second = JacobianInverse(domain)
    if inv_first:
        first, second = second, first
    if v_first:
        expression = v @ first @ second
    else:
        expression = first @ second @ v

    simpler_expression = simplify(expression)

    assert isinstance(expression, MatrixProduct)
    assert isinstance(simpler_expression, TestFunction)


def test_jacobian_and_inverse_form(lagrange_element):
    """Test that Jacobian and inverse Jacobian are successfully cancelled."""
    pytest.xfail()

    element = lagrange_element("triangle", 2, (2,))
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)

    j = Jacobian(domain)
    j_inv_t = JacobianInverseTranspose(domain)

    form = inner(j @ u, j_inv_t @ v) * dx
    simpler_form = simplify(form)

    assert isinstance(form, Integral)
    assert isinstance(simpler_form, Integral)

    assert isinstance(form.integrand, Inner)
    assert isinstance(form.integrand.first, MatrixProduct)
    assert isinstance(form.integrand.first._items[0], Jacobian)
    assert isinstance(form.integrand.first._items[1], TrialFunction)
    assert isinstance(form.integrand.second, MatrixProduct)
    assert isinstance(form.integrand.second._items[0], JacobianInverseTranspose)
    assert isinstance(form.integrand.second._items[1], TestFunction)

    assert isinstance(simpler_form.integrand, Inner)
    if isinstance(simpler_form.integrand.first, TrialFunction):
        assert isinstance(simpler_form.integrand.second, TestFunction)
    else:
        assert isinstance(simpler_form.integrand.first, TestFunction)
        assert isinstance(simpler_form.integrand.second, TrialFunction)


def test_commutative_operands_are_sorted(lagrange_element):
    """Test that sums and products equal up to the order of operands simplify equally."""
    element = lagrange_element("triangle", 2)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    f = Coefficient(space)
    g = Coefficient(space)

    assert u * f + g != g + f * u
    assert simplify(u * f + g) == simplify(g + f * u)


def test_a_fingerprint_is_the_same_in_every_process():
    """What the digest is for: ordering that does not depend on the hash seed.

    Python salts the hashes of strings per process, so two processes
    with different seeds disagree on ``hash``. The fingerprint of an
    expression must not, or two MPI ranks would sort the operands of a
    sum differently.
    """
    script = (
        "from uflx.algorithms.simplify import Fingerprints\n"
        "from uflx.expressions import RealScalar\n"
        "print(Fingerprints()(RealScalar(1.0) * RealScalar(2.0) + RealScalar(3.0)).hex())"
    )
    printed = {
        seed: subprocess.run(
            [sys.executable, "-c", script],
            env={**os.environ, "PYTHONHASHSEED": seed},
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        for seed in ("1", "2")
    }
    assert printed["1"] == printed["2"]
    assert len(printed["1"]) == 16


def test_the_digest_is_eight_bytes_and_unambiguous():
    """Sized for ordering and not for an adversary, with parts that cannot run together."""
    assert len(digest("a")) == 8
    assert digest("ab", "c") != digest("a", "bc")
    assert digest("a", "") != digest("a")
