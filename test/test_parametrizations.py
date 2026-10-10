"""Test parametrizations of a domain, and compositions of them."""

import pytest
from conftest import Interval

from uflx import composed_domain, parametrized_domain
from uflx.algorithms import simplify
from uflx.basis_functions import EvaluatedBasisFunction
from uflx.domains import (
    RD,
    AbstractCoordinateDomain,
    AbstractParametrization,
    IdentityParametrization,
    entity_domain,
)
from uflx.expressions import AbstractExpression, Integer, MatrixProduct, RealScalar
from uflx.functions import AbstractVariable
from uflx.geometry import Jacobian, PushedForwardPoint
from uflx.graphs import as_graph
from uflx.parametrizations import ComposedParametrization, FiniteElementParametrization
from uflx.points import Point
from uflx.tensors import Matrix, Vector


class Parabolic(AbstractParametrization):
    """The analytic map y -> (y, y^2) from the line into the plane.

    Nothing here describes a basis or a degree of freedom. It exists to
    check that a map can be given in closed form.
    """

    @property
    def source(self) -> AbstractCoordinateDomain:
        """This map starts on the line."""
        return RD(1)

    @property
    def target_dimension(self) -> int:
        """This map lands in the plane."""
        return 2

    def value(self, point: AbstractVariable) -> AbstractExpression:
        """Square the coordinate to get the second component."""
        y = point.component(0)
        return Vector([y, y * y])

    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """Differentiate (y, y^2) by hand."""
        y = point.component(0)
        return Matrix([[Integer(1)], [Integer(2) * y]])

    def __repr__(self) -> str:
        """Representation."""
        return "Parabolic()"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, Parabolic)

    def __hash__(self) -> int:
        """Hash."""
        return hash("uflx.test.Parabolic")


@pytest.fixture
def parabola():
    """The analytic map, on its own."""
    return Parabolic()


@pytest.fixture
def line(lagrange_element):
    """A finite element mesh of the line."""
    return parametrized_domain(lagrange_element("interval", 1, (1,)))


@pytest.fixture
def line_map(line):
    """The finite element map of the mesh of the line."""
    return line.parametrization(Interval())


@pytest.fixture
def mesh_on_a_parabola(line, parabola):
    """A finite element mesh of the line, mapped onto a parabola in the plane."""
    return composed_domain(line, parabola)


@pytest.fixture
def entity_point():
    """A point of an interval's coordinate domain."""
    return Point([RealScalar(0.25)], entity_domain(Interval()))


def test_a_closed_form_map_needs_no_element(parabola, entity_point):
    """A map can be evaluated and differentiated without a basis."""
    y = entity_point.component(0)
    point = Point([y], RD(1))

    assert parabola.value(point) == Vector([y, y * y])
    assert parabola.jacobian(point).value_shape == (2, 1)
    assert not parabola.is_affine
    assert not parabola.is_identity


def test_an_element_map_knows_its_source_and_target(line_map):
    """A finite element map starts on its cell and lands in the ambient coordinates."""
    parametrization = line_map

    assert isinstance(parametrization, FiniteElementParametrization)
    assert parametrization.source == entity_domain(Interval())
    assert parametrization.target_dimension == 1
    assert parametrization.is_affine


def test_composing_dimensions_takes_one_from_each_half(mesh_on_a_parabola):
    """The topology comes from the mesh and the ambient dimension from the map."""
    assert mesh_on_a_parabola.topological_dimension == 1
    assert mesh_on_a_parabola.geometric_dimension == 2
    assert mesh_on_a_parabola.cell_types == (Interval(),)


def test_composing_incompatible_maps_is_rejected(line, line_map):
    """A map must start where the one before it lands."""

    class FromThePlane(Parabolic):
        @property
        def source(self):
            """Start somewhere the line does not land."""
            return RD(2)

    with pytest.raises(ValueError, match="Cannot compose"):
        ComposedParametrization(line_map, FromThePlane())

    with pytest.raises(ValueError, match="Cannot map a domain"):
        composed_domain(line, FromThePlane())


def test_a_composed_push_forward_squares_the_inner_map(mesh_on_a_parabola, entity_point):
    """The second ambient coordinate is the first one squared."""
    x = PushedForwardPoint(entity_point, mesh_on_a_parabola)
    expanded = x.expand_geometry()

    assert isinstance(expanded, Point)
    assert expanded.domain == RD(2)
    s = expanded.component(0)
    assert expanded.component(1) == s * s


def test_a_composed_jacobian_is_the_chain_rule(mesh_on_a_parabola, line_map, entity_point):
    """The composite's Jacobian is the outer Jacobian times the inner one."""
    parametrization = mesh_on_a_parabola

    j = Jacobian(parametrization, entity_point).expand_geometry()

    assert isinstance(j, MatrixProduct)
    assert j.value_shape == (2, 1)

    inner = line_map
    s = inner.value(entity_point).component(0)
    ds = inner.jacobian(entity_point).component(0, 0)
    assert simplify(j.component(0, 0)) == simplify(ds)
    assert simplify(j.component(1, 0)) == simplify(Integer(2) * s * ds)


def test_a_composed_domain_still_tabulates_its_inner_basis(mesh_on_a_parabola, entity_point):
    """Basis functions and a closed form map appear in the same expression."""
    j = Jacobian(mesh_on_a_parabola, entity_point).expand_geometry()

    assert any(isinstance(n, EvaluatedBasisFunction) for n in as_graph(j))


def test_composing_with_the_identity_changes_nothing(parabola, entity_point):
    """The identity map is a unit for composition, at least up to simplification."""
    composed = ComposedParametrization(IdentityParametrization(RD(1)), parabola)
    point = Point([entity_point.component(0)], RD(1))

    assert composed.target_dimension == parabola.target_dimension
    assert composed.source == RD(1)
    assert composed.value(point) == parabola.value(point)
    assert simplify(composed.jacobian(point)) == simplify(parabola.jacobian(point))


def test_affineness_and_identity_conjoin(line, line_map, parabola):
    """A composite is affine only if both halves are, and likewise the identity."""
    identity = IdentityParametrization(RD(1))

    assert ComposedParametrization(identity, identity).is_identity
    assert ComposedParametrization(identity, identity).is_affine
    assert not ComposedParametrization(identity, parabola).is_affine
    assert not ComposedParametrization(line_map, parabola).is_affine
    assert not composed_domain(line, parabola).has_affine_parametrization


def test_parametrizations_compare_by_value(line, parabola):
    """Geometric quantities hold a map, so rewriting and simplifying compare them."""
    cell = Interval()
    assert line.parametrization(cell) == line.parametrization(cell)
    assert hash(line.parametrization(cell)) == hash(line.parametrization(cell))

    composed = composed_domain(line, parabola)
    assert composed.parametrization(cell) == composed.parametrization(cell)
    assert hash(composed.parametrization(cell)) == hash(composed.parametrization(cell))
    assert composed.parametrization(cell) != line.parametrization(cell)


def test_an_element_map_still_exposes_its_element(lagrange_element):
    """Describing a map with an element is one option, and it stays reachable."""
    element = lagrange_element("triangle", 1, (2,))
    domain = parametrized_domain(element)
    (cell,) = domain.cell_types

    assert domain.parametrization_element(cell) == element
    assert domain.parametrization(cell).element == element
