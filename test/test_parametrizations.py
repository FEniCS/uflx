"""Test parametrizations of a domain."""

import pytest
from conftest import Interval

from uflx import parametrized_domain
from uflx.basis_functions import EvaluatedBasisFunction
from uflx.domains import RD, AbstractParametrizedDomain, entity_domain
from uflx.entities import AbstractEntity
from uflx.expressions import AbstractExpression, Integer, RealScalar
from uflx.geometry import Jacobian, PushedForwardPoint
from uflx.graphs import as_graph
from uflx.points import AbstractPoint, Point
from uflx.tensors import zero


class Parabola(AbstractParametrizedDomain):
    """The unit interval mapped to the parabola (X, X^2), with no finite element.

    Nothing here describes a basis or a degree of freedom. It exists to
    check that a domain can give its parametrization in closed form.
    """

    def __init__(self, cell: AbstractEntity):
        """Initialise."""
        self._cell = cell

    @property
    def geometric_dimension(self) -> int:
        """The parabola sits in the plane."""
        return 2

    @property
    def topological_dimension(self) -> int:
        """The parabola is a curve."""
        return 1

    @property
    def cell_types(self) -> tuple[AbstractEntity, ...]:
        """The one cell type."""
        return (self._cell,)

    def parametrization_component(
        self,
        cell: AbstractEntity,
        point: AbstractPoint,
        component: int,
        derivative: tuple[int, ...],
    ) -> AbstractExpression:
        """Differentiate (X, X^2) by hand."""
        (order,) = derivative
        coordinate = point.component(0)
        if component == 0:
            derivatives = [coordinate, Integer(1)]
        else:
            derivatives = [coordinate * coordinate, Integer(2) * coordinate, Integer(2)]
        if order >= len(derivatives):
            return zero(())
        return derivatives[order]

    def __repr__(self) -> str:
        """Representation."""
        return f"Parabola({self._cell!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, Parabola) and self._cell == other._cell

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.test.Parabola", self._cell))


@pytest.fixture
def parabola():
    """A parametrized domain with no element behind it."""
    return Parabola(Interval())


@pytest.fixture
def entity_point(parabola):
    """A point of the parabola's cell's coordinate domain."""
    (cell,) = parabola.cell_types
    return Point([RealScalar(0.25)], entity_domain(cell))


def test_push_forward_needs_no_element(parabola, entity_point):
    """Pushing a point forward reads the parametrization, not an element."""
    x = PushedForwardPoint(entity_point, parabola).expand_geometry()

    coordinate = entity_point.component(0)
    assert x == Point([coordinate, coordinate * coordinate], RD(2))


def test_jacobian_needs_no_element(parabola, entity_point):
    """The Jacobian differentiates the parametrization, whatever describes it."""
    j = Jacobian(parabola, entity_point).expand_geometry()

    coordinate = entity_point.component(0)
    assert j.value_shape == (2, 1)
    assert j.component(0, 0) == Integer(1)
    assert j.component(1, 0) == Integer(2) * coordinate


def test_an_analytic_parametrization_tabulates_nothing(parabola, entity_point):
    """No basis function evaluation appears, because there is no basis."""
    for expression in [
        PushedForwardPoint(entity_point, parabola).expand_geometry(),
        Jacobian(parabola, entity_point).expand_geometry(),
    ]:
        assert not any(isinstance(n, EvaluatedBasisFunction) for n in as_graph(expression))


def test_a_parametrization_is_not_affine_unless_it_says_so(parabola, lagrange_element):
    """The default is conservative, and the element flavour answers from its element."""
    assert not parabola.has_affine_parametrization
    assert parametrized_domain(lagrange_element("triangle", 1, (2,))).has_affine_parametrization


def test_an_element_parametrization_still_exposes_its_element(lagrange_element):
    """Describing the map with an element is one option, and it stays reachable."""
    element = lagrange_element("triangle", 1, (2,))
    domain = parametrized_domain(element)
    (cell,) = domain.cell_types

    assert domain.parametrization_element(cell) == element
