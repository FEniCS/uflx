"""Test pulling back the gradient of a function to the reference cell."""

import pytest
from test_parametrizations import Paraboloidal

from uflx import Coefficient, composed_domain, function_space, grad, parametrized_domain
from uflx.algorithms import pull_back_to_entity, simplify
from uflx.expressions import MatrixProduct
from uflx.geometry import JacobianInverse, SpatialCoordinate, TangentialProjector
from uflx.operators import EntityGrad
from uflx.tensors import Identity

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


def test_the_gradient_of_the_spatial_coordinate_is_the_identity(lagrange_element):
    """``grad(x) = J K``, and on a cell that is the identity.

    The entity gradient of the coordinates a map lands on is that map's
    derivative, so nothing is differentiated here: the answer is the
    Jacobian, and the pull back multiplies its inverse in as it does for
    any gradient. That the product collapses is the simplifier's, and it
    was already there.
    """
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    x = SpatialCoordinate(domain)

    assert grad(x).value_shape == (2, 2)
    assert simplify(pull_back_to_entity(grad(x))) == Identity(2)


def test_the_gradient_of_the_spatial_coordinate_on_a_manifold_is_the_projector(
    lagrange_element,
):
    """Where the map is not square, ``J K`` is the tangential projector instead.

    This is the one closed form check that the pseudo-inverse ``grad``
    uses on a manifold is the right object. The gradient of the
    coordinates cannot be the identity on a surface -- a surface has two
    directions to move in and three coordinates -- and what it is instead
    is the projection onto the tangent plane, which is what ``J K``
    comes to.
    """
    flat = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    surface = composed_domain(flat, Paraboloidal())
    x = SpatialCoordinate(surface)

    assert grad(x).value_shape == (3, 3)

    pulled_back = simplify(pull_back_to_entity(grad(x)))

    assert pulled_back == TangentialProjector(surface)
    assert pulled_back != Identity(3)


def test_the_gradient_of_an_expression_is_refused(lagrange_element):
    """A field and the spatial coordinate have an entity gradient; a product does not.

    Differentiating one needs a product rule, which is the symbolic
    differentiation LANGUAGE.md section 4 describes and is a different
    piece of work. Refusing says so rather than asserting.
    """
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    u = Coefficient(function_space(domain, lagrange_element("triangle", 1)))

    with pytest.raises(ValueError, match="needs a rule for differentiating it"):
        grad(u * u)
    with pytest.raises(ValueError, match="needs a rule for differentiating it"):
        grad(SpatialCoordinate(domain).component(0))
