"""Test the square root operator."""

import pytest

from uflx import sqrt
from uflx.expressions import ComplexScalar, RealScalar, Sqrt
from uflx.tensors import Vector


def test_sqrt_real():
    """sqrt() of a real scalar."""
    s = sqrt(RealScalar(2.25))
    assert isinstance(s, Sqrt)
    assert s.value_shape == ()
    assert s.as_float() == pytest.approx(1.5, rel=1e-15)


def test_sqrt_complex():
    """sqrt() of a complex scalar."""
    s = sqrt(ComplexScalar(RealScalar(-4.0), RealScalar(0.0)))
    assert s.as_complex() == pytest.approx(2j, rel=1e-15)


def test_sqrt_rejects_non_scalar():
    """sqrt() is only defined for scalar expressions."""
    with pytest.raises(ValueError):
        sqrt(Vector([RealScalar(1.0), RealScalar(4.0)]))


def test_sqrt_equality_and_hash():
    """Two square roots of the same expression are equal and hash equally."""
    assert sqrt(RealScalar(2.0)) == sqrt(RealScalar(2.0))
    assert hash(sqrt(RealScalar(2.0))) == hash(sqrt(RealScalar(2.0)))
    assert sqrt(RealScalar(2.0)) != sqrt(RealScalar(3.0))
