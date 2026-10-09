"""Test forms."""

from uflx import function_space, parametrized_domain
from uflx.basis_functions import (
    AbstractEvaluatedBasisFunction,
    EvaluatedBasisFunction,
)
from uflx.expressions import RealScalar
from uflx.points import Point


def test_physical_basis_function(lagrange_element):
    """Test physical basis function."""
    element = lagrange_element("triangle", 1)
    domain = parametrized_domain(lagrange_element("triangle", 1, (3,)))
    space = function_space(domain, element)

    point = Point([RealScalar(1.0)] * 3, in_entity_coordinates=False)

    phys_f = EvaluatedBasisFunction(space, 0, point)

    assert phys_f.derivative == (0, 0)
    assert phys_f.domain_size == 2
    assert phys_f.value_shape == ()

    d1 = phys_f.diff(1)
    d11 = phys_f.diff(1).diff(1)
    d101 = phys_f.diff(1).diff(0).diff(1)
    assert isinstance(d1, AbstractEvaluatedBasisFunction) and not d1.in_entity_coordinates
    assert isinstance(d11, AbstractEvaluatedBasisFunction) and not d11.in_entity_coordinates
    assert isinstance(d101, AbstractEvaluatedBasisFunction) and not d101.in_entity_coordinates
    assert d1.derivative == (0, 1)
    assert d11.derivative == (0, 2)
    assert d101.derivative == (1, 2)


def test_reference_basis_function(lagrange_element):
    """Test reference basis function."""
    element = lagrange_element("triangle", 1)
    domain = parametrized_domain(lagrange_element("triangle", 1, (3,)))
    space = function_space(domain, element)

    point = Point([RealScalar(1.0)] * 3, in_entity_coordinates=True)

    ref_f = EvaluatedBasisFunction(space, 0, point)

    assert ref_f.derivative == (0, 0)
    assert ref_f.domain_size == 2
    assert ref_f.value_shape == ()

    d1 = ref_f.diff(1)
    d11 = ref_f.diff(1).diff(1)
    d101 = ref_f.diff(1).diff(0).diff(1)
    assert isinstance(d1, AbstractEvaluatedBasisFunction) and d1.in_entity_coordinates
    assert isinstance(d11, AbstractEvaluatedBasisFunction) and d11.in_entity_coordinates
    assert isinstance(d101, AbstractEvaluatedBasisFunction) and d101.in_entity_coordinates
    assert d1.derivative == (0, 1)
    assert d11.derivative == (0, 2)
    assert d101.derivative == (1, 2)
