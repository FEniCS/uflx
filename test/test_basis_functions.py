"""Test forms."""

from uflx import function_space, parametrized_domain
from uflx.basis_functions import (
    AbstractEvaluatedBasisFunction,
    EvaluatedBasisFunction,
)
from uflx.domains import RD, EntityDomain
from uflx.expressions import RealScalar
from uflx.points import Point


def test_ambient_basis_function(lagrange_element):
    """Test a basis function evaluated at a point in ambient coordinates."""
    element = lagrange_element("triangle", 1)
    domain = parametrized_domain(lagrange_element("triangle", 1, (3,)))
    space = function_space(domain, element)

    ambient_point = Point([RealScalar(1.0)] * 3, RD(3))

    phys_f = EvaluatedBasisFunction(space, 0, ambient_point)

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


def test_entity_basis_function(lagrange_element):
    """Test a basis function evaluated at a point in the entity's coordinates."""
    element = lagrange_element("triangle", 1)
    domain = parametrized_domain(lagrange_element("triangle", 1, (3,)))
    space = function_space(domain, element)

    entity_point = Point([RealScalar(0.25)] * 2, EntityDomain(domain.cell_types[0]))

    ref_f = EvaluatedBasisFunction(space, 0, entity_point)

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
