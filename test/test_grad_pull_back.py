"""Test pulling back the gradient of a function to the reference cell."""

import pytest

from uflx import Coefficient, coordinate_element, function_space, grad
from uflx.algorithms import pull_back_to_reference
from uflx.expressions import MatrixProduct
from uflx.geometry import JacobianInverse
from uflx.operators import ReferenceGrad

cells_and_gdims = [
    ("interval", 1),
    ("interval", 2),
    ("triangle", 2),
    ("triangle", 3),
    ("quadrilateral", 3),
    ("tetrahedron", 3),
]


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
@pytest.mark.parametrize("vector", [False, True])
def test_grad_pull_back(cell, gdim, vector, lagrange_element):
    """grad(w) = ReferenceGrad(w_ref) K, for scalar and vector functions."""
    domain = coordinate_element(lagrange_element(cell, 1, (gdim,)))
    shape = (gdim,) if vector else None
    w = Coefficient(function_space(domain, lagrange_element(cell, 1, shape)))

    pulled_back = pull_back_to_reference(grad(w))
    assert isinstance(pulled_back, MatrixProduct)
    assert pulled_back.value_shape == grad(w).value_shape
    first, second = pulled_back.init_args[0]
    assert isinstance(first, ReferenceGrad)
    assert isinstance(second, JacobianInverse)
