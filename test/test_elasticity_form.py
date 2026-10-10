"""Test that the full isotropic linear-elasticity weak form composes end to end.

This module instead builds the real bilinear form

    a(u, v) = inner(sigma(u), sym(grad(v))) * dx

with real vector-valued TrialFunction/TestFunction arguments and applies the
pull_back_to_entity, apply_push_forwards, and simplify algorithms.
"""

import pytest

from uflx import TestFunction, TrialFunction, dx, function_space, grad, inner, parametrized_domain
from uflx.algorithms import pull_back_to_entity, simplify
from uflx.functions import Argument
from uflx.graphs import as_graph
from uflx.integrals import Integral
from uflx.maps import apply_push_forwards
from uflx.operators import Grad, Tr, Transpose, sym, tr
from uflx.tensors import Identity


def _sigma(displacement, lambda_, mu):
    """Isotropic Hooke's law, with hard-wired (literal float) Lame parameters."""
    strain = sym(grad(displacement))
    d = displacement.value_shape[0]
    return lambda_ * tr(strain) * Identity(d) + 2 * mu * strain


@pytest.mark.parametrize(("cell", "dim"), [("triangle", 2), ("tetrahedron", 3)])
def test_elasticity_bilinear_form_composes(lagrange_element, cell, dim):
    """inner(sigma(u), sym(grad(v))) * dx should build without error and stay scalar."""
    domain = parametrized_domain(lagrange_element(cell, 1, (dim,)))
    space = function_space(domain, lagrange_element(cell, 1, (dim,)))

    u = TrialFunction(space)
    v = TestFunction(space)
    assert u.value_shape == (dim,)

    lambda_, mu = 1.7, 0.8
    sigma_u = _sigma(u, lambda_, mu)
    assert sigma_u.value_shape == (dim, dim)

    form = inner(sigma_u, sym(grad(v))) * dx(domain)
    assert isinstance(form, Integral)
    assert form.integrand.value_shape == ()


@pytest.mark.parametrize(("cell", "dim"), [("triangle", 2), ("tetrahedron", 3)])
def test_elasticity_bilinear_form_pulls_back_to_reference(lagrange_element, cell, dim):
    """The whole form must pull back to the entity's coordinates and stay well-shaped."""
    domain = parametrized_domain(lagrange_element(cell, 1, (dim,)))
    space = function_space(domain, lagrange_element(cell, 1, (dim,)))

    u = TrialFunction(space)
    v = TestFunction(space)
    lambda_, mu = 1.7, 0.8

    form = inner(_sigma(u, lambda_, mu), sym(grad(v))) * dx(domain)

    pulled_back = pull_back_to_entity(form)
    assert isinstance(pulled_back, Integral)
    assert pulled_back.integrand.value_shape == ()

    pushed_forward = apply_push_forwards(pulled_back)
    simplified = simplify(pushed_forward)
    assert isinstance(simplified, Integral)
    assert simplified.integrand.value_shape == ()

    nodes = list(as_graph(simplified))

    # Grad works only on physical (non-reference) arguments; by this point in
    # the real pipeline every Grad must have been replaced by a reference-cell
    # equivalent, and every Argument must be reference-valued.
    assert not any(isinstance(n, Grad) for n in nodes)
    assert not any(isinstance(n, Argument) and not n.in_entity_coordinates for n in nodes)

    # The tensor-op nodes themselves must have survived the round trip -- not
    # been silently dropped or left wrapping a stale (physical) Grad.
    assert any(isinstance(n, Tr) for n in nodes)
    assert any(isinstance(n, Transpose) for n in nodes)
