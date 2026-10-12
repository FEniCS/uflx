"""Test expressions."""

import pytest

from uflx.expressions import Div, Integer, RealScalar


@pytest.mark.parametrize("zero", [Integer(0), RealScalar(0.0)])
def test_dividing_by_a_zero_scalar_is_refused(zero):
    """Division by a zero scalar raises, whether the zero is an integer or real."""
    with pytest.raises(ZeroDivisionError):
        Div(Integer(1), zero)
