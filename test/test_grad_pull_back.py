"""Test pulling back the gradient of a function to the reference cell."""

import pytest

from uflx import Coefficient, function_space, grad, parametrized_domain
from uflx.algorithms import pull_back_to_entity
from uflx.expressions import MatrixProduct
from uflx.geometry import JacobianInverse
from uflx.operators import EntityGrad

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
    """grad(w) = EntityGrad(w_ref) K, for scalar and vector functions."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    shape = (gdim,) if vector else None
    w = Coefficient(function_space(domain, lagrange_element(cell, 1, shape)))

    pulled_back = pull_back_to_entity(grad(w))
    assert isinstance(pulled_back, MatrixProduct)
    assert pulled_back.value_shape == grad(w).value_shape
    first, second = pulled_back.init_args[0]
    assert isinstance(first, EntityGrad)
    assert isinstance(second, JacobianInverse)
