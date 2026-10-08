"""Test complex values."""

import numpy as np

from uflx import Coefficient, coordinate_element, function_space, inner
from uflx.complex import conj, take_imaginary_part, take_real_part
from uflx.expressions import ComplexScalar, Integer, RealScalar, to_scalar


def test_real_part():
    """Test taking real part."""
    assert take_real_part(ComplexScalar(Integer(4), Integer(6))) == 4


def test_imaginary_part():
    """Test taking imaginary part."""
    assert take_imaginary_part(ComplexScalar(Integer(4), Integer(6))) == 6


def test_complex_scalar_add():
    """Test addition with complex scalars."""
    five = to_scalar(5)

    z = five + 2j
    assert np.isclose(z.re.as_float(), 5)
    assert np.isclose(z.im.as_float(), 2)

    z = 2j + five
    assert np.isclose(z.re.as_float(), 5)
    assert np.isclose(z.im.as_float(), 2)


def test_complex_scalar_sub():
    """Test subtraction with complex scalars."""
    five = to_scalar(5)

    z = five - 2j
    assert np.isclose(z.re.as_float(), 5)
    assert np.isclose(z.im.as_float(), -2)

    z = 2j - five
    assert np.isclose(z.re.as_float(), -5)
    assert np.isclose(z.im.as_float(), 2)


def test_complex_scalar_mult():
    """Test multiplication with complex scalars."""
    five = to_scalar(5)

    z = five * 2j
    assert np.isclose(z.re.as_float(), 0)
    assert np.isclose(z.im.as_float(), 10)

    z = 2j * five
    assert np.isclose(z.re.as_float(), 0)
    assert np.isclose(z.im.as_float(), 10)


def test_complex_scalar_div():
    """Test division with complex scalars."""
    five = to_scalar(5)

    z = five / 2j
    assert np.isclose(z.re.as_float(), 0)
    assert np.isclose(z.im.as_float(), -2.5)

    z = 2j / five
    assert np.isclose(z.re.as_float(), 0)
    assert np.isclose(z.im.as_float(), 0.4)


def test_complex_scalar_neg():
    """Test negation with complex scalars."""
    z = to_scalar(3 - 2j)

    assert np.isclose(z.re.as_float(), 3)
    assert np.isclose(z.im.as_float(), -2)

    assert np.isclose((-z).re.as_float(), -3)
    assert np.isclose((-z).im.as_float(), 2)


def test_conj_of_complex_scalar():
    """conj(1 + 2j) = 1 - 2j."""
    z = ComplexScalar(RealScalar(1.0), RealScalar(2.0))
    assert conj(z).as_complex() == 1 - 2j


def test_inner_of_complex_scalars():
    """inner(z, z) = |z|^2."""
    z = ComplexScalar(RealScalar(1.0), RealScalar(2.0))
    assert inner(z, z).as_complex() == 5
    i = ComplexScalar(RealScalar(0.0), RealScalar(1.0))
    assert inner(i, i).as_complex() == 1


def test_conj_of_real_function(lagrange_element):
    """A real-valued function is its own conjugate."""
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    u = Coefficient(function_space(domain, lagrange_element("triangle", 1)))
    assert conj(u) == u
