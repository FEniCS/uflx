"""Test the exterior boundary of a domain, built as a domain in its own right.

A boundary is not a new kind of object. The exterior boundary of a cell is
a set of domains of codimension one, each the image of a reference facet
under its inclusion into the cell followed by the cell's own map. That is
a composition, which UFLx already has, so every quantity a domain carries
-- its metric, its volume element, its normal, its tangential projector --
is available on a facet with no new machinery.

Nothing here involves a finite element. The inclusions are closed form
affine maps written out by hand, which is what makes the geometry
concrete: a facet's length comes out as a number that can be checked
against the number it should be. A mesh would supply these inclusions from
its reference cell, and the point of doing it by hand first is to see
exactly what it has to supply.

What is missing is in :func:`test_a_facets_normal_is_normal_but_not_outward`.
"""

import math
from collections.abc import Sequence

import pytest
from conftest import Interval, Triangle
from conftest import Point as PointEntity

from uflx import composed_domain
from uflx.domains import (
    RD,
    AbstractCoordinateDomain,
    AbstractParametrization,
    entity_domain,
)
from uflx.expressions import AbstractExpression, RealScalar, expression_sum
from uflx.functions import AbstractVariable
from uflx.geometry import (
    Jacobian,
    MetricTensor,
    SpatialCoordinate,
    TangentialProjector,
    UnitNormal,
    VolumeElement,
    _as_dense_matrix,
)
from uflx.points import Point
from uflx.tensors import Matrix, Vector


class AffineMap(AbstractParametrization):
    """A closed form affine map, given by its matrix and its offset.

    Enough to describe the inclusion of a reference facet into a cell, or
    a linear map of the coordinates a domain lives in, without a finite
    element anywhere near it.
    """

    def __init__(
        self,
        source: AbstractCoordinateDomain,
        matrix: Sequence[Sequence[float]],
        offset: Sequence[float],
    ):
        """Initialise.

        Args:
            source: The region this map starts from
            matrix: The rows of the matrix of the map
            offset: Where the source's origin lands
        """
        self._source = source
        self._matrix = tuple(tuple(row) for row in matrix)
        self._offset = tuple(offset)
        if len({len(row) for row in self._matrix}) != 1:
            raise ValueError("Every row of the matrix must be the same length.")
        if len(self._offset) != len(self._matrix):
            raise ValueError("The offset must have one entry per row of the matrix.")
        if len(self._matrix[0]) != source.geometric_dimension:
            raise ValueError("The matrix must have one column per coordinate of the source.")

    @property
    def source(self) -> AbstractCoordinateDomain:
        """The region this map starts from."""
        return self._source

    @property
    def target_dimension(self) -> int:
        """This map lands in as many coordinates as the matrix has rows."""
        return len(self._matrix)

    def apply(self, coordinates: Sequence[float]) -> list[float]:
        """Apply the map to ordinary numbers, to say what is expected of it."""
        return [
            sum(a * x for a, x in zip(row, coordinates, strict=True)) + b
            for row, b in zip(self._matrix, self._offset, strict=True)
        ]

    def value(self, point: AbstractVariable) -> AbstractExpression:
        """Multiply the matrix by the point and add the offset."""
        return Vector(
            [
                expression_sum(
                    (RealScalar(a) * point.component(j) for j, a in enumerate(row)),
                    default=RealScalar(0.0),
                )
                + RealScalar(b)
                for row, b in zip(self._matrix, self._offset, strict=True)
            ]
        )

    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """An affine map's derivative is its matrix."""
        return Matrix([[RealScalar(a) for a in row] for row in self._matrix])

    @property
    def is_affine(self) -> bool:
        """An affine map is affine."""
        return True

    def __repr__(self) -> str:
        """Representation."""
        return f"AffineMap({self._source!r}, {self._matrix}, {self._offset})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, AffineMap)
            and other._source == self._source
            and other._matrix == self._matrix
            and other._offset == self._offset
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.test.AffineMap", self._source, self._matrix, self._offset))


# The reference triangle's vertices. AbstractEntity says which vertices a
# facet has and not where they are, so these coordinates come from outside,
# as they would from a mesh's reference cell.
reference_vertices = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))

# Facet k of the reference triangle, its length, and the outward normal of
# the cell along it. Facet k's vertices are Triangle().sub_entity_vertices(1)[k].
facet_lengths = (1.0, 1.0, math.sqrt(2.0))
outward_normals = ((0.0, -1.0), (-1.0, 0.0), (1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0)))


def facet_inclusion(facet: int) -> AffineMap:
    """The inclusion of the reference interval into facet ``facet`` of a triangle.

    Affine, taking 0 to the facet's first vertex and 1 to its second, so
    the reference interval covers the facet once.
    """
    first, second = Triangle().sub_entity_vertices(1)[facet]
    start, end = reference_vertices[first], reference_vertices[second]
    return AffineMap(
        entity_domain(Interval()), tuple((end[i] - start[i],) for i in range(2)), start
    )


@pytest.fixture
def cell():
    """The reference triangle, as a domain."""
    return entity_domain(Triangle())


@pytest.fixture
def facets():
    """The three facets of the reference triangle, each a domain of its own."""
    return [composed_domain(entity_domain(Interval()), facet_inclusion(k)) for k in range(3)]


@pytest.fixture
def on_the_facet():
    """A point of the reference interval, away from either end."""
    return Point([RealScalar(0.3)], entity_domain(Interval()))


@pytest.fixture
def in_the_cell():
    """A point of the reference triangle's coordinate domain."""
    return Point([RealScalar(0.2), RealScalar(0.3)], entity_domain(Triangle()))


@pytest.fixture
def stretch():
    """A map of the plane that scales the two axes differently."""
    return AffineMap(RD(2), ((2.0, 0.0), (0.0, 1.0)), (0.0, 0.0))


@pytest.fixture
def scaling():
    """A uniform scaling of the plane by three."""
    return AffineMap(RD(2), ((3.0, 0.0), (0.0, 3.0)), (0.0, 0.0))


def tangent_of(domain, point) -> list[float]:
    """The one tangent of a domain of topological dimension one, as numbers."""
    jacobian = _as_dense_matrix(Jacobian(domain, point).expand_geometry())
    return [jacobian.component(i, 0).as_float() for i in range(2)]


def normal_of(domain, point) -> list[float]:
    """A domain's unit normal, as numbers."""
    return [UnitNormal(domain, point).component(i).as_float() for i in range(2)]


def measure_of(domain, point) -> float:
    """A domain's volume element, as a number."""
    return VolumeElement(domain, point).expand_geometry().as_float()


def test_a_facet_is_a_domain_of_codimension_one(cell, facets, on_the_facet):
    """Each facet lives in the cell's coordinates while being one dimensional."""
    assert cell.topological_dimension == 2
    assert cell.geometric_dimension == 2

    for facet in facets:
        assert facet.topological_dimension == 1
        assert facet.geometric_dimension == 2
        assert Jacobian(facet, on_the_facet).value_shape == (2, 1)


def test_a_facet_covers_the_edge_between_its_two_vertices(facets, on_the_facet):
    """At 0.3 along, a facet is three tenths of the way from one vertex to the other."""
    for k, facet in enumerate(facets):
        first, second = Triangle().sub_entity_vertices(1)[k]
        start, end = reference_vertices[first], reference_vertices[second]
        expected = [start[i] + 0.3 * (end[i] - start[i]) for i in range(2)]

        coordinates = SpatialCoordinate(facet, on_the_facet).expand_geometry()

        assert [coordinates.component(i).as_float() for i in range(2)] == pytest.approx(expected)


@pytest.mark.parametrize("k", range(3))
def test_the_measure_of_a_facet_is_its_length(facets, on_the_facet, k):
    """Two edges of the reference triangle have length one and the hypotenuse root two.

    The volume element is a density rather than an integral, so this is
    the length the inclusion gives a unit length of the reference
    interval, which here is constant because the inclusion is affine.
    """
    facet = facets[k]

    assert measure_of(facet, on_the_facet) == pytest.approx(facet_lengths[k])
    assert MetricTensor(facet, on_the_facet).expand_geometry().component(0, 0).as_float() == (
        pytest.approx(facet_lengths[k] ** 2)
    )
    assert math.sqrt(sum(c * c for c in tangent_of(facet, on_the_facet))) == pytest.approx(
        facet_lengths[k]
    )


def test_the_measure_of_the_cell_itself_is_one(cell, in_the_cell):
    """The reference cell's own map is the identity, so it rescales nothing.

    Its area is a half, which is not this: a volume element is the factor
    an integral picks up, and integrating one over the reference triangle
    is what gives the half.
    """
    assert measure_of(cell, in_the_cell) == pytest.approx(1.0)


@pytest.mark.parametrize("k", range(3))
def test_a_facets_normal_is_normal_but_not_outward(facets, on_the_facet, k):
    """Every facet gets a unit normal across it, with a sign nothing can fix yet.

    The tangent space of a facet of codimension one fixes the normal up
    to sign with no further information, and that much works. Which of
    the two directions points out of the cell does not follow from the
    facet alone: it needs to know which side the cell is on, and a facet
    built this way is not told.

    What is actually returned is the generalised cross product of the
    inclusion's columns, whose sign follows the facet's vertex ordering.
    Facets 0 and 2 of the reference triangle happen to come out outward
    and facet 1 comes out inward, because the vertex pairs
    AbstractEntity gives are not consistently oriented around the cell.
    Fixing this is what an orientation would be for.
    """
    facet = facets[k]
    normal = normal_of(facet, on_the_facet)

    assert sum(c * c for c in normal) == pytest.approx(1.0)
    assert sum(a * b for a, b in zip(normal, tangent_of(facet, on_the_facet))) == pytest.approx(0.0)

    agreement = sum(a * b for a, b in zip(normal, outward_normals[k]))
    assert abs(agreement) == pytest.approx(1.0)
    assert agreement == pytest.approx(1.0 if k in (0, 2) else -1.0)


def test_the_normals_of_a_cells_facets_do_not_all_point_the_same_way_out(facets, on_the_facet):
    """Stated once over the cell, rather than per facet: the signs disagree.

    This is the one thing an exterior facet integral needs that is not
    here. A form that wants the outward normal cannot be written until
    something says which side of each facet its cell is on.
    """
    agreements = [
        sum(a * b for a, b in zip(normal_of(facet, on_the_facet), outward_normals[k]))
        for k, facet in enumerate(facets)
    ]

    assert all(abs(a) == pytest.approx(1.0) for a in agreements)
    assert not all(a > 0 for a in agreements)


@pytest.mark.parametrize("k", range(3))
def test_a_facets_projector_is_the_identity_less_its_normal(facets, on_the_facet, k):
    """On a domain of codimension one the tangential projector is I - n (x) n.

    The quantity's docstring says so for a hypersurface, and a facet of a
    triangle is the smallest hypersurface there is.
    """
    facet = facets[k]
    projector = _as_dense_matrix(TangentialProjector(facet, on_the_facet).expand_geometry())
    normal = normal_of(facet, on_the_facet)

    expected = [
        [(1.0 if i == j else 0.0) - normal[i] * normal[j] for j in range(2)] for i in range(2)
    ]
    entries = [[projector.component(i, j).as_float() for j in range(2)] for i in range(2)]

    assert [e for row in entries for e in row] == pytest.approx(
        [e for row in expected for e in row]
    )


def test_scaling_a_cell_scales_its_area_by_c_squared_and_its_facets_by_c(
    cell, facets, in_the_cell, on_the_facet, scaling
):
    """The measure of a domain scales as c to its own topological dimension.

    Which is the whole difference between a cell integral and an
    exterior facet integral over the same geometry, and it falls out of
    the one measure rather than being built into two.
    """
    scaled_cell = composed_domain(cell, scaling)

    assert measure_of(scaled_cell, in_the_cell) == pytest.approx(3.0**2 * 1.0)

    for k, facet in enumerate(facets):
        scaled_facet = composed_domain(facet, scaling)

        assert measure_of(scaled_facet, on_the_facet) == pytest.approx(3.0**1 * facet_lengths[k])


@pytest.mark.parametrize("k", range(3))
def test_a_facets_normal_follows_the_cofactor_rule(facets, on_the_facet, stretch, k):
    """Carrying a facet through a map of the plane carries its normal by J^-T.

    A normal is not a tangent, so it does not move with the map. Under a
    scaling that treats the axes differently the two behave visibly
    differently, which a uniform scaling would hide.
    """
    facet = facets[k]
    moved = composed_domain(facet, stretch)
    normal = normal_of(facet, on_the_facet)

    carried = stretch.apply(normal)
    inverse_transpose = [normal[0] / 2.0, normal[1] / 1.0]
    length = math.sqrt(sum(c * c for c in inverse_transpose))
    expected = [c / length for c in inverse_transpose]

    assert normal_of(moved, on_the_facet) == pytest.approx(expected)
    if k == 2:
        # The hypotenuse's normal is not along an axis, so the two rules differ.
        assert normal_of(moved, on_the_facet) != pytest.approx(carried)


def test_reparametrizing_a_facet_keeps_its_measure_and_flips_its_normal(facets, on_the_facet):
    """The hypotenuse covered backwards is the same facet, measured the same way.

    Two charts of one facet differ by a reparametrization of the
    reference interval, and that reparametrization is an isometry of it,
    so the density they induce agrees. The normal does not: the
    generalised cross product follows the chart. This is what an
    interior facet integral will have to reckon with, where the two
    charts are the two cells' own.
    """
    backwards = composed_domain(
        entity_domain(Interval()),
        AffineMap(entity_domain(Interval()), ((1.0,), (-1.0,)), (0.0, 1.0)),
    )

    assert measure_of(backwards, on_the_facet) == pytest.approx(facet_lengths[2])
    assert normal_of(backwards, on_the_facet) == pytest.approx(
        [-c for c in normal_of(facets[2], on_the_facet)]
    )


def test_the_boundary_of_an_interval_is_two_points_with_no_normal():
    """A facet of an interval has topological dimension zero, and no normal.

    The normal to a point is a sign, which is a convention rather than
    something to compute, so the quantity says so rather than guessing.
    """
    end = composed_domain(
        entity_domain(PointEntity()), AffineMap(entity_domain(PointEntity()), ((),), (1.0,))
    )
    origin = Point([], entity_domain(PointEntity()))

    assert end.topological_dimension == 0
    assert end.geometric_dimension == 1

    with pytest.raises(NotImplementedError, match="topological dimension zero is a sign"):
        UnitNormal(end, origin).value_shape


@pytest.mark.xfail(raises=AssertionError, strict=True, reason="Matrix([[]]) cannot be built")
def test_the_measure_of_a_point_is_the_counting_measure():
    """A zero dimensional domain has nothing to measure, so its density is one.

    The volume element of a point should be 1, the empty Gram
    determinant being 1. Instead the empty matrix its Jacobian needs
    cannot be built, and tensors.py's shape inference asserts. Harmless
    until a measure over a point domain is wanted, which is what a point
    evaluation is.
    """
    end = composed_domain(
        entity_domain(PointEntity()), AffineMap(entity_domain(PointEntity()), ((),), (1.0,))
    )
    origin = Point([], entity_domain(PointEntity()))

    assert measure_of(end, origin) == pytest.approx(1.0)
