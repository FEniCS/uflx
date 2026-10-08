"""Test equality of measures."""

from uflx import dx
from uflx.integrals import Measure


def test_equal_measures():
    """Measures with the same attributes are equal and hash equally."""
    assert Measure(codim=0) == dx
    assert hash(Measure(codim=0)) == hash(dx)
    assert Measure(codim=1) != dx
    assert Measure(codim=1) != Measure(codim=1, boundary_only=True)
