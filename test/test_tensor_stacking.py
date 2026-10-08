"""Test tensors with non-scalar entries."""

import pytest

from uflx import Coefficient, coordinate_element, function_space, grad
from uflx.expressions import RealScalar
from uflx.tensors import Matrix, Tensor, Vector


def test_stacked_shape(lagrange_element):
    """A tensor of entries of shape (2,) has shape (n, 2)."""
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    p = Coefficient(function_space(domain, lagrange_element("triangle", 1)))

    assert Tensor([grad(p), grad(p), grad(p)]).value_shape == (3, 2)


def test_stacked_components():
    """Component (i, j) of a stack of vectors is component j of entry i."""
    a, b, c, d = (RealScalar(float(i)) for i in range(4))
    t = Tensor([Vector([a, b]), [c, d]])
    assert t.value_shape == (2, 2)
    assert t.component(0, 1) == b
    assert t.component(1, 0) == c


def test_vector_and_matrix_entries_are_scalars(lagrange_element):
    """Vector and Matrix keep scalar entries, and raise for stacked ones."""
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    p = Coefficient(function_space(domain, lagrange_element("triangle", 1)))

    with pytest.raises(ValueError, match="scalars"):
        Vector([grad(p), grad(p)])
    with pytest.raises(ValueError, match="scalars"):
        Matrix([[grad(p)], [grad(p)]])


def test_entries_of_different_shapes(lagrange_element):
    """Entries of different shapes cannot be stacked."""
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    p = Coefficient(function_space(domain, lagrange_element("triangle", 1)))

    with pytest.raises(ValueError, match="cannot be stacked"):
        Tensor([p, grad(p)])
