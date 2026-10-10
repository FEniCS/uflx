"""Test parametrizations of a domain."""

import pytest
from conftest import Interval

from uflx import parametrized_domain
from uflx.algorithms import simplify
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


class Parabolic:
    """The analytic map y -> (y, y^2) from the line into the plane."""

    @property
    def target_dimension(self) -> int:
        """The map lands in the plane."""
        return 2

    def value(self, component: int, y: AbstractExpression) -> AbstractExpression:
        """One component of the map at y."""
        return y if component == 0 else y * y

    def derivative(self, component: int, y: AbstractExpression) -> AbstractExpression:
        """The derivative of one component of the map at y."""
        return Integer(1) if component == 0 else Integer(2) * y


class ComposedDomain(AbstractParametrizedDomain):
    """A parametrized domain post-composed with an analytic map.

    The inner domain is a mesh: its parametrization is a finite element
    basis summed against coordinate dofs. The outer map is a closed form
    chart out of the ambient coordinates of that mesh. Both appear in the
    composite's parametrization, combined by the chain rule.
    """

    def __init__(self, inner: AbstractParametrizedDomain, outer: Parabolic):
        """Initialise."""
        if inner.geometric_dimension != 1:
            raise ValueError("This map composes with a mesh of the line.")
        self._inner = inner
        self._outer = outer

    @property
    def geometric_dimension(self) -> int:
        """The composite lands wherever the outer map lands."""
        return self._outer.target_dimension

    @property
    def topological_dimension(self) -> int | None:
        """Composing with a chart does not change the topology."""
        return self._inner.topological_dimension

    @property
    def cell_types(self) -> tuple[AbstractEntity, ...]:
        """The cells are the inner mesh's cells."""
        return self._inner.cell_types

    def parametrization_component(
        self,
        cell: AbstractEntity,
        point: AbstractPoint,
        component: int,
        derivative: tuple[int, ...],
    ) -> AbstractExpression:
        """Apply the chain rule to the inner parametrization and the outer map."""
        (order,) = derivative
        inner_value = self._inner.parametrization_component(cell, point, 0, (0,))
        if order == 0:
            return self._outer.value(component, inner_value)
        if order == 1:
            inner_derivative = self._inner.parametrization_component(cell, point, 0, (1,))
            return self._outer.derivative(component, inner_value) * inner_derivative
        raise NotImplementedError("The chain rule here stops at first derivatives.")

    def __repr__(self) -> str:
        """Representation."""
        return f"ComposedDomain({self._inner!r}, {self._outer!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, ComposedDomain)
            and self._inner == other._inner
            and type(self._outer) is type(other._outer)
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.test.ComposedDomain", self._inner, type(self._outer)))


@pytest.fixture
def mesh_on_a_parabola(lagrange_element):
    """A finite element mesh of the line, mapped onto a parabola in the plane."""
    return ComposedDomain(parametrized_domain(lagrange_element("interval", 1, (1,))), Parabolic())


def test_a_composed_domain_has_the_dimensions_of_both_halves(mesh_on_a_parabola):
    """The topology comes from the mesh and the ambient dimension from the chart."""
    assert mesh_on_a_parabola.topological_dimension == 1
    assert mesh_on_a_parabola.geometric_dimension == 2


def test_a_composed_push_forward_squares_the_inner_map(mesh_on_a_parabola, lagrange_element):
    """The second ambient coordinate is the first one squared."""
    (cell,) = mesh_on_a_parabola.cell_types
    X = Point([RealScalar(0.25)], entity_domain(cell))

    x = PushedForwardPoint(X, mesh_on_a_parabola).expand_geometry()

    assert isinstance(x, Point)
    assert x.domain == RD(2)
    s = x.component(0)
    assert x.component(1) == s * s


def test_a_composed_jacobian_obeys_the_chain_rule(mesh_on_a_parabola):
    """The Jacobian is (s', 2 s s') for the inner map s."""
    (cell,) = mesh_on_a_parabola.cell_types
    X = Point([RealScalar(0.25)], entity_domain(cell))
    inner = mesh_on_a_parabola._inner

    j = Jacobian(mesh_on_a_parabola, X).expand_geometry()

    s = inner.parametrization_component(cell, X, 0, (0,))
    ds = inner.parametrization_component(cell, X, 0, (1,))
    assert j.value_shape == (2, 1)
    # The chain rule multiplies the first component by one, which simplify removes.
    assert simplify(j.component(0, 0)) == simplify(ds)
    assert simplify(j.component(1, 0)) == simplify(Integer(2) * s * ds)


def test_a_composed_domain_still_tabulates_its_inner_basis(mesh_on_a_parabola):
    """Basis functions and a closed form map appear in the same expression."""
    (cell,) = mesh_on_a_parabola.cell_types
    X = Point([RealScalar(0.25)], entity_domain(cell))

    j = Jacobian(mesh_on_a_parabola, X).expand_geometry()
    nodes = list(as_graph(j))

    assert any(isinstance(n, EvaluatedBasisFunction) for n in nodes)


def test_a_composed_domain_refuses_derivatives_it_cannot_take(mesh_on_a_parabola):
    """The chain rule is only implemented as far as the Jacobian needs."""
    (cell,) = mesh_on_a_parabola.cell_types
    X = Point([RealScalar(0.25)], entity_domain(cell))

    with pytest.raises(NotImplementedError):
        mesh_on_a_parabola.parametrization_component(cell, X, 1, (2,))
