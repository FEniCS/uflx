"""Test parametrizations of a domain, and compositions of them."""

import math

import pytest
from conftest import Interval, Triangle

from uflx import (
    Coefficient,
    TestFunction,
    composed_domain,
    dx,
    function_space,
    inner,
    parametrized_domain,
)
from uflx.algorithms import pull_back_to_entity, simplify
from uflx.basis_functions import EvaluatedBasisFunction
from uflx.domains import (
    RD,
    AbstractCoordinateDomain,
    AbstractParametrization,
    EntityDomain,
    IdentityParametrization,
    entity_domain,
)
from uflx.expressions import (
    AbstractExpression,
    Integer,
    MatrixProduct,
    RealScalar,
    expression_sum,
)
from uflx.functions import AbstractVariable
from uflx.geometry import (
    AbstractGeometricQuantity,
    ImagePoint,
    Jacobian,
    JacobianDeterminant,
    JacobianInverse,
    MetricTensor,
    SpatialCoordinate,
    TangentialProjector,
    UnitNormal,
    VolumeElement,
    _as_dense_matrix,
    expand_geometry,
)
from uflx.graphs import as_graph
from uflx.integrals import Integral
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


class Paraboloidal(AbstractParametrization):
    """The analytic map (u, v) -> (u, v, u^2 + v^2) from the plane into space.

    The surface analogue of :class:`Parabolic`, and the graph of a
    function as that one is, so its area element is
    ``sqrt(1 + |grad f|^2)`` and its normal is ``(-grad f, 1)`` divided by
    the same thing. Unlike a curve it has a metric with off-diagonal
    entries, a tangent plane rather than a tangent line, and three minors
    in its normal rather than two.
    """

    @property
    def source(self) -> AbstractCoordinateDomain:
        """This map starts in the plane."""
        return RD(2)

    @property
    def target_dimension(self) -> int:
        """This map lands in space."""
        return 3

    def value(self, point: AbstractVariable) -> AbstractExpression:
        """Sum the squares of the coordinates to get the third component."""
        u, v = point.component(0), point.component(1)
        return Vector([u, v, u * u + v * v])

    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """Differentiate (u, v, u^2 + v^2) by hand."""
        u, v = point.component(0), point.component(1)
        return Matrix(
            [
                [Integer(1), Integer(0)],
                [Integer(0), Integer(1)],
                [Integer(2) * u, Integer(2) * v],
            ]
        )

    def __repr__(self) -> str:
        """Representation."""
        return "Paraboloidal()"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, Paraboloidal)

    def __hash__(self) -> int:
        """Hash."""
        return hash("uflx.test.Paraboloidal")


class LinearMap(AbstractParametrization):
    """A linear map of the coordinates it starts in, given by its matrix.

    A closed form map to apply after another one, so that a curve or a
    surface can be composed with something downstream of it. A rotation
    and a reflection are isometries and a uniform scaling is not, which
    is what makes them worth composing with: the metric a map induces has
    to notice the difference, and the measure that follows from it has to
    notice the same difference.
    """

    def __init__(self, entries: tuple[tuple[float, ...], ...]):
        """Initialise.

        Args:
            entries: The rows of the square matrix of the map
        """
        self._entries = tuple(tuple(row) for row in entries)
        if any(len(row) != len(self._entries) for row in self._entries):
            raise ValueError("A linear map of its own coordinates has a square matrix.")

    @property
    def source(self) -> AbstractCoordinateDomain:
        """This map starts where it lands."""
        return RD(len(self._entries))

    @property
    def target_dimension(self) -> int:
        """This map lands where it starts."""
        return len(self._entries)

    @property
    def determinant(self) -> float:
        """The determinant of the matrix, which says if it reverses orientation.

        Expanded along the first row in ordinary arithmetic, this being
        what the test expects of UFLx rather than something it asks UFLx
        for.
        """

        def expand(rows: tuple[tuple[float, ...], ...]) -> float:
            if len(rows) == 1:
                return rows[0][0]
            total = 0.0
            for j, entry in enumerate(rows[0]):
                minor = tuple(row[:j] + row[j + 1 :] for row in rows[1:])
                total += (-1.0) ** j * entry * expand(minor)
            return total

        return expand(self._entries)

    def apply(self, coordinates: list[float]) -> list[float]:
        """Apply the map to ordinary numbers, to say what is expected of it."""
        return [sum(a * x for a, x in zip(row, coordinates, strict=True)) for row in self._entries]

    def value(self, point: AbstractVariable) -> AbstractExpression:
        """Multiply the matrix by the point."""
        return Vector(
            [
                expression_sum(RealScalar(a) * point.component(j) for j, a in enumerate(row))
                for row in self._entries
            ]
        )

    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """A linear map is its own derivative."""
        return Matrix([[RealScalar(a) for a in row] for row in self._entries])

    @property
    def is_affine(self) -> bool:
        """A linear map is affine."""
        return True

    def __repr__(self) -> str:
        """Representation."""
        return f"LinearMap({self._entries})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, LinearMap) and other._entries == self._entries

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.test.LinearMap", self._entries))


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


@pytest.fixture
def curve(parabola):
    """The parabola alone, as a domain: one cell, mapped in closed form.

    A mesh of the line carries coordinate dofs into every Jacobian, which
    no simplification turns into a number. Starting from the reference
    interval instead keeps the geometry concrete, so it can be evaluated
    and checked rather than matched against the expression it should be.
    """
    return composed_domain(entity_domain(Interval()), parabola)


@pytest.fixture
def paraboloid():
    """The analytic surface map, on its own."""
    return Paraboloidal()


@pytest.fixture
def surface(paraboloid):
    """The paraboloid alone, as a domain: one cell, mapped in closed form.

    The surface counterpart of ``curve``, and concrete for the same
    reason.
    """
    return composed_domain(entity_domain(Triangle()), paraboloid)


@pytest.fixture
def rotation():
    """A quarter turn of the plane, which is an isometry."""
    return LinearMap(((0.0, -1.0), (1.0, 0.0)))


@pytest.fixture
def reflection():
    """A reflection of the plane, an isometry that reverses orientation."""
    return LinearMap(((1.0, 0.0), (0.0, -1.0)))


@pytest.fixture
def scaling():
    """A uniform scaling of the plane by three, which is not an isometry."""
    return LinearMap(((3.0, 0.0), (0.0, 3.0)))


@pytest.fixture
def space_rotation():
    """A quarter turn of space about the third axis, which is an isometry."""
    return LinearMap(((0.0, -1.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)))


@pytest.fixture
def space_reflection():
    """A reflection of space, an isometry that reverses orientation."""
    return LinearMap(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)))


@pytest.fixture
def space_scaling():
    """A uniform scaling of space by three, which is not an isometry."""
    return LinearMap(((3.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0)))


def arc_length_element(y: float) -> float:
    """The arc length element of y -> (y, y^2), whose tangent is (1, 2y)."""
    return math.sqrt(1.0 + 4.0 * y * y)


def area_element(u: float, v: float) -> float:
    """The area element of (u, v) -> (u, v, u^2 + v^2).

    The graph of a function has ``sqrt(1 + |grad f|^2)``, and here
    ``grad f`` is ``(2u, 2v)``.
    """
    return math.sqrt(1.0 + 4.0 * u * u + 4.0 * v * v)


def on_the_interval(y: float) -> Point:
    """A point of an interval's coordinate domain."""
    return Point([RealScalar(y)], entity_domain(Interval()))


def on_the_triangle(u: float, v: float) -> Point:
    """A point of a triangle's coordinate domain."""
    return Point([RealScalar(u), RealScalar(v)], entity_domain(Triangle()))


def columns_of(quantity, rows: int, cols: int) -> list[list[float]]:
    """Read a matrix valued quantity column by column, as numbers."""
    matrix = _as_dense_matrix(quantity.expand_geometry())
    return [[matrix.component(i, j).as_float() for i in range(rows)] for j in range(cols)]


def entries_of(quantity, rows: int, cols: int) -> list[list[float]]:
    """Read a matrix valued quantity row by row, as numbers."""
    matrix = _as_dense_matrix(quantity.expand_geometry())
    return [[matrix.component(i, j).as_float() for j in range(cols)] for i in range(rows)]


def flat(matrix: list[list[float]]) -> list[float]:
    """Flatten a matrix, pytest.approx not comparing nested sequences."""
    return [entry for row in matrix for entry in row]


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
    x = ImagePoint(entity_point, mesh_on_a_parabola)
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


def test_restricting_a_composed_domain_keeps_its_map(mesh_on_a_parabola, parabola):
    """Restriction reaches the mesh underneath; the chart is unchanged by it."""
    (cell,) = mesh_on_a_parabola.cell_types

    restricted = mesh_on_a_parabola.restricted_to(cell)

    assert restricted.cell_types == (cell,)
    assert restricted.geometric_dimension == 2
    assert restricted.topological_dimension == 1
    assert isinstance(restricted.parametrization(cell), ComposedParametrization)


def test_a_composed_domain_has_a_measure(mesh_on_a_parabola, entity_point):
    """A composed map's Jacobian is a matrix product, which still has to reduce.

    The volume element, the metric and the projector all have to write the
    Jacobian out before working on it, and a chain rule does not hand
    them a Matrix.
    """
    for quantity in [VolumeElement, MetricTensor, TangentialProjector]:
        expanded = quantity(mesh_on_a_parabola, entity_point).expand_geometry()
        assert expanded.value_shape == quantity(mesh_on_a_parabola, entity_point).value_shape

    assert VolumeElement(mesh_on_a_parabola, entity_point).value_shape == ()
    assert MetricTensor(mesh_on_a_parabola, entity_point).value_shape == (1, 1)
    assert TangentialProjector(mesh_on_a_parabola, entity_point).value_shape == (2, 2)


def test_a_form_over_a_composed_domain_expands(mesh_on_a_parabola, lagrange_element):
    """Pulling a form back onto a composed domain leaves no geometry behind."""
    space = function_space(mesh_on_a_parabola, lagrange_element("interval", 1))
    form = inner(Coefficient(space), TestFunction(space)) * dx(mesh_on_a_parabola)

    expanded = expand_geometry(pull_back_to_entity(form))

    assert not any(isinstance(n, AbstractGeometricQuantity) for n in as_graph(expanded))


def test_the_normal_of_a_curve_is_a_unit_vector_across_the_tangent(parabola, entity_point):
    """Check the normal numerically, on a map whose Jacobian is concrete.

    An element map's Jacobian sums over coordinate dofs, which no amount
    of simplification turns into a number. A closed form map's does not,
    so here the normal can be evaluated and checked rather than matched
    against the expression it is expected to be.
    """
    curve = composed_domain(entity_domain(Interval()), parabola)

    jacobian = _as_dense_matrix(Jacobian(curve, entity_point).expand_geometry())
    tangent = [jacobian.component(i, 0).as_float() for i in range(2)]
    normal = [UnitNormal(curve, entity_point).component(i).as_float() for i in range(2)]

    # y = x^2 at x = 0.25, so the tangent is (1, 1/2).
    assert tangent == pytest.approx([1.0, 0.5])
    assert sum(c * c for c in normal) == pytest.approx(1.0)
    assert sum(a * b for a, b in zip(tangent, normal, strict=True)) == pytest.approx(0.0)


def test_the_spatial_coordinate_of_a_curve_is_the_point_on_it(parabola, entity_point):
    """On y = x^2 at x = 1/4, x is (1/4, 1/16)."""
    curve = composed_domain(entity_domain(Interval()), parabola)
    coordinates = SpatialCoordinate(curve, entity_point)

    assert coordinates.value_shape == (2,)
    assert coordinates[0].expand_geometry().as_float() == pytest.approx(0.25)
    assert coordinates[1].expand_geometry().as_float() == pytest.approx(0.0625)


@pytest.mark.parametrize("y", [0.0, 0.25, 1.0])
def test_the_volume_element_of_a_curve_is_its_arc_length_element(curve, y):
    """On y -> (y, y^2) the measure is sqrt(1 + 4y^2), and nothing is a determinant.

    The map is not square, so it has no determinant to take. What an
    integral over the curve picks up is the length the map gives a unit
    length of the interval, which is the norm of its one tangent.
    """
    point = on_the_interval(y)

    jacobian = _as_dense_matrix(Jacobian(curve, point).expand_geometry())
    tangent = [jacobian.component(i, 0).as_float() for i in range(2)]
    metric = MetricTensor(curve, point).expand_geometry()

    assert VolumeElement(curve, point).expand_geometry().as_float() == pytest.approx(
        arc_length_element(y)
    )
    assert math.sqrt(sum(c * c for c in tangent)) == pytest.approx(arc_length_element(y))
    assert math.sqrt(metric.component(0, 0).as_float()) == pytest.approx(arc_length_element(y))

    with pytest.raises(ValueError, match="has no determinant"):
        JacobianDeterminant(curve, point).expand_geometry()


@pytest.mark.parametrize("isometry", ["rotation", "reflection"])
def test_an_ambient_isometry_leaves_a_curves_measure_alone(curve, isometry, request):
    """Composing with an isometry of the plane cannot change a length on the curve.

    The metric is the ambient one pulled back, so a map that preserves
    the ambient one preserves it, and the measure that follows from it.
    The normal is not so indifferent: the generalised cross product
    tracks orientation, so a map that reverses orientation flips the
    normal it gives while leaving the measure untouched. That difference
    is the whole of what a density discards.
    """
    isometry = request.getfixturevalue(isometry)
    point = on_the_interval(0.25)
    moved = composed_domain(curve, isometry)

    assert VolumeElement(moved, point).expand_geometry().as_float() == pytest.approx(
        VolumeElement(curve, point).expand_geometry().as_float()
    )
    assert MetricTensor(moved, point).expand_geometry().component(0, 0).as_float() == pytest.approx(
        MetricTensor(curve, point).expand_geometry().component(0, 0).as_float()
    )

    normal = [UnitNormal(curve, point).component(i).as_float() for i in range(2)]
    moved_normal = [UnitNormal(moved, point).component(i).as_float() for i in range(2)]
    carried = isometry.apply(normal)
    sign = 1.0 if isometry.determinant > 0 else -1.0

    assert moved_normal == pytest.approx([sign * c for c in carried])


def test_an_ambient_scaling_scales_a_curves_measure(curve, scaling):
    """A uniform scaling by c stretches a length by c and leaves a direction alone.

    The metric picks up c^2, being quadratic in the map, and the measure
    picks up c^tdim, which on a curve is c itself.
    """
    point = on_the_interval(0.25)
    scaled = composed_domain(curve, scaling)

    assert VolumeElement(scaled, point).expand_geometry().as_float() == pytest.approx(
        3.0 * arc_length_element(0.25)
    )
    assert MetricTensor(scaled, point).expand_geometry().component(
        0, 0
    ).as_float() == pytest.approx(9.0 * (1.0 + 4.0 * 0.25 * 0.25))

    normal = [UnitNormal(curve, point).component(i).as_float() for i in range(2)]
    scaled_normal = [UnitNormal(scaled, point).component(i).as_float() for i in range(2)]

    assert scaled_normal == pytest.approx(normal)


def test_composing_twice_applies_the_maps_in_turn(curve, rotation):
    """A composite composes again, so a curve can be carried onward.

    The reference interval, then the parabola, then a quarter turn: three
    maps, and the point of the plane they land on is the parabola's point
    rotated.
    """
    point = on_the_interval(0.25)
    rotated = composed_domain(curve, rotation)
    (cell,) = rotated.cell_types

    assert isinstance(rotated.parametrization(cell), ComposedParametrization)
    assert rotated.topological_dimension == 1
    assert rotated.geometric_dimension == 2
    assert Jacobian(rotated, point).value_shape == (2, 1)

    on_the_curve = SpatialCoordinate(curve, point).expand_geometry()
    on_the_curve = [on_the_curve.component(i).as_float() for i in range(2)]
    rotated_coordinates = SpatialCoordinate(rotated, point).expand_geometry()

    assert on_the_curve == pytest.approx([0.25, 0.0625])
    assert [rotated_coordinates.component(i).as_float() for i in range(2)] == pytest.approx(
        rotation.apply(on_the_curve)
    )


def test_pulling_a_form_over_a_parabola_back_uses_the_volume_element(
    mesh_on_a_parabola, lagrange_element
):
    """A manifold's measure is its volume element, which is not a determinant.

    A form over a mesh carried onto the parabola has a non-square
    Jacobian, so pulling it back cannot reach for a determinant. The
    measure follows the integral onto the cell it is pulled back to.
    """
    space = function_space(mesh_on_a_parabola, lagrange_element("interval", 1))
    form = inner(Coefficient(space), TestFunction(space)) * dx(mesh_on_a_parabola)

    pulled = pull_back_to_entity(form)

    assert isinstance(pulled, Integral)
    assert pulled.measure.domain == EntityDomain(Interval())
    nodes = list(as_graph(pulled))
    assert any(isinstance(n, VolumeElement) for n in nodes)
    assert not any(isinstance(n, JacobianDeterminant) for n in nodes)


surface_points = [(0.0, 0.0), (0.25, 0.5), (0.5, 0.5)]


@pytest.mark.parametrize(("u", "v"), surface_points)
def test_the_volume_element_of_a_surface_is_its_area_element(surface, u, v):
    """On (u, v) -> (u, v, u^2 + v^2) the measure is sqrt(1 + 4u^2 + 4v^2).

    A surface has a metric with off-diagonal entries, where a curve's is
    a single number, so this is the first check that the Gram
    determinant of a two by two metric is taken as a number and not only
    as the right shape.
    """
    point = on_the_triangle(u, v)
    metric = entries_of(MetricTensor(surface, point), 2, 2)

    assert VolumeElement(surface, point).expand_geometry().as_float() == pytest.approx(
        area_element(u, v)
    )
    assert flat(metric) == pytest.approx(
        [1.0 + 4.0 * u * u, 4.0 * u * v, 4.0 * u * v, 1.0 + 4.0 * v * v]
    )
    determinant = metric[0][0] * metric[1][1] - metric[0][1] * metric[1][0]
    assert math.sqrt(determinant) == pytest.approx(area_element(u, v))

    with pytest.raises(ValueError, match="has no determinant"):
        JacobianDeterminant(surface, point).expand_geometry()


@pytest.mark.parametrize(("u", "v"), surface_points)
def test_the_normal_of_a_surface_is_the_graphs_gradient_formula(surface, u, v):
    """The graph of a function has normal (-grad f, 1), over the area element.

    Three minors rather than a curve's two, and the middle one is the
    negated term, which is where an alternating sign is easiest to lose.
    """
    point = on_the_triangle(u, v)
    normal = [UnitNormal(surface, point).component(i).as_float() for i in range(3)]
    tangents = columns_of(Jacobian(surface, point), 3, 2)

    assert normal == pytest.approx(
        [-2.0 * u / area_element(u, v), -2.0 * v / area_element(u, v), 1.0 / area_element(u, v)]
    )
    assert sum(c * c for c in normal) == pytest.approx(1.0)
    for tangent in tangents:
        assert sum(a * b for a, b in zip(normal, tangent, strict=True)) == pytest.approx(0.0)


def test_the_tangential_projector_keeps_a_surfaces_tangent_plane(surface):
    """It keeps the two directions the surface can move in and discards the third.

    Its rank is the topological dimension, which for a projection is its
    trace, so on a surface in space the trace is two and not three.
    """
    point = on_the_triangle(0.25, 0.5)
    projector = entries_of(TangentialProjector(surface, point), 3, 3)
    normal = [UnitNormal(surface, point).component(i).as_float() for i in range(3)]
    tangents = columns_of(Jacobian(surface, point), 3, 2)

    def apply(matrix, vector):
        return [sum(row[k] * vector[k] for k in range(3)) for row in matrix]

    assert sum(projector[i][i] for i in range(3)) == pytest.approx(2.0)
    assert apply(projector, normal) == pytest.approx([0.0, 0.0, 0.0])
    for tangent in tangents:
        assert apply(projector, tangent) == pytest.approx(tangent)
    for i in range(3):
        assert apply(projector, projector[i]) == pytest.approx(projector[i])
        for j in range(3):
            assert projector[i][j] == pytest.approx(projector[j][i])


def test_the_pseudo_inverse_inverts_a_surface_map_on_its_tangent_plane(surface):
    """J+ J is the identity on the cell, and J J+ is the projector, not the identity.

    The two products differ in which space they act on, and only the
    first is an identity. Taking the second for one is the manifold
    mistake worth a number rather than a shape.
    """
    point = on_the_triangle(0.25, 0.5)
    jacobian = entries_of(Jacobian(surface, point), 3, 2)
    pseudo_inverse = entries_of(JacobianInverse(surface, point), 2, 3)
    projector = entries_of(TangentialProjector(surface, point), 3, 3)

    on_the_cell = [
        [sum(pseudo_inverse[i][k] * jacobian[k][j] for k in range(3)) for j in range(2)]
        for i in range(2)
    ]
    in_space = [
        [sum(jacobian[i][k] * pseudo_inverse[k][j] for k in range(2)) for j in range(3)]
        for i in range(3)
    ]

    assert flat(on_the_cell) == pytest.approx([1.0, 0.0, 0.0, 1.0])
    assert flat(in_space) == pytest.approx(flat(projector))


@pytest.mark.parametrize("isometry", ["space_rotation", "space_reflection"])
def test_an_ambient_isometry_leaves_a_surfaces_measure_alone(surface, isometry, request):
    """The same law as for a curve, with three components in the normal.

    An isometry of space preserves the metric the surface induces and so
    its area element, while the normal picks up the sign of the map's
    determinant.
    """
    isometry = request.getfixturevalue(isometry)
    point = on_the_triangle(0.25, 0.5)
    moved = composed_domain(surface, isometry)

    assert VolumeElement(moved, point).expand_geometry().as_float() == pytest.approx(
        area_element(0.25, 0.5)
    )
    assert flat(entries_of(MetricTensor(moved, point), 2, 2)) == pytest.approx(
        flat(entries_of(MetricTensor(surface, point), 2, 2))
    )

    normal = [UnitNormal(surface, point).component(i).as_float() for i in range(3)]
    moved_normal = [UnitNormal(moved, point).component(i).as_float() for i in range(3)]
    carried = isometry.apply(normal)
    sign = 1.0 if isometry.determinant > 0 else -1.0

    assert moved_normal == pytest.approx([sign * c for c in carried])


def test_an_ambient_scaling_scales_a_surfaces_measure(surface, space_scaling):
    """A scaling by c multiplies the measure by c to the topological dimension.

    On a curve that was c itself. On a surface it is c squared, the
    metric picking up c^2 in each of its two directions and the area
    element being the square root of their determinant.
    """
    point = on_the_triangle(0.25, 0.5)
    scaled = composed_domain(surface, space_scaling)

    assert scaled.topological_dimension == 2
    assert VolumeElement(scaled, point).expand_geometry().as_float() == pytest.approx(
        3.0**2 * area_element(0.25, 0.5)
    )
    assert flat(entries_of(MetricTensor(scaled, point), 2, 2)) == pytest.approx(
        [9.0 * e for e in flat(entries_of(MetricTensor(surface, point), 2, 2))]
    )

    normal = [UnitNormal(surface, point).component(i).as_float() for i in range(3)]
    scaled_normal = [UnitNormal(scaled, point).component(i).as_float() for i in range(3)]

    assert scaled_normal == pytest.approx(normal)
