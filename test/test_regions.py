"""Test a domain that is a region of R^d, with no cells in it at all.

Everything else in the suite starts from a reference cell. A cell is a
finite element idea -- it is what an element is attached to -- and nothing
about integrating over a set requires one. A chart is a map out of a
parameter region, and a reference cell is one kind of parameter region.

So here the domain is a box, charted by itself, and its boundary is its
faces, each charted by the box with one coordinate dropped. No
AbstractEntity, no sub_entity_vertices, no element, no reference cell.
Box and BoxFace are in conftest, consumers being where concrete domains
live.

The interesting result is in
:func:`test_opposite_faces_are_indistinguishable_to_the_normal`.
"""

import pytest
from conftest import Box

from uflx.domains import RD, AbstractCellularDomain, AbstractCoordinateDomain
from uflx.expressions import RealScalar
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


@pytest.fixture
def rectangle():
    """The region [0, 2] x [0, 3] of the plane."""
    return Box(((0.0, 2.0), (0.0, 3.0)))


@pytest.fixture
def brick():
    """The region [0, 2] x [0, 3] x [0, 5] of space."""
    return Box(((0.0, 2.0), (0.0, 3.0), (0.0, 5.0)))


def point_in(region: AbstractCoordinateDomain, *coordinates: float) -> Point:
    """A point of a region, named by its coordinates."""
    return Point([RealScalar(c) for c in coordinates], region)


def entries_of(quantity, rows: int, cols: int) -> list[float]:
    """Read a matrix valued quantity as a flat list of numbers."""
    matrix = _as_dense_matrix(quantity.expand_geometry())
    return [matrix.component(i, j).as_float() for i in range(rows) for j in range(cols)]


def normal_of(domain, point) -> list[float]:
    """A domain's unit normal, as numbers."""
    gdim = domain.geometric_dimension
    return [UnitNormal(domain, point).component(i).as_float() for i in range(gdim)]


def test_a_region_is_a_domain_with_no_cells_in_it(rectangle):
    """A box is charted, not cellular. Nothing asks it for a cell type."""
    assert not isinstance(rectangle, AbstractCellularDomain)
    assert not hasattr(rectangle, "cell_types")

    assert rectangle.geometric_dimension == 2
    assert rectangle.topological_dimension == 2
    assert rectangle.chart_sources == (rectangle,)


def test_a_regions_own_chart_is_the_identity(rectangle):
    """Its points are coordinate tuples already, so nothing maps them anywhere.

    The density is one and not the area: a volume element is the factor
    an integral picks up, and the area of six comes from integrating that
    one over the region, which is a quadrature's business and not a
    domain's.
    """
    point = point_in(rectangle, 0.5, 1.5)

    assert rectangle.chart(rectangle).is_identity
    assert VolumeElement(rectangle, point).expand_geometry().as_float() == pytest.approx(1.0)
    assert entries_of(MetricTensor(rectangle, point), 2, 2) == pytest.approx([1.0, 0.0, 0.0, 1.0])


def test_a_region_refuses_a_chart_it_does_not_have(rectangle):
    """Geometry at a point of the wrong region is not this domain's to give."""
    with pytest.raises(ValueError, match="is charted by itself"):
        VolumeElement(rectangle, point_in(RD(2), 0.5, 1.5)).expand_geometry()


def test_a_face_is_a_domain_of_codimension_one(rectangle):
    """Each face lives in the region's coordinates while being one lower."""
    for face in rectangle.faces:
        assert face.geometric_dimension == 2
        assert face.topological_dimension == 1
        assert face.chart_sources == (face.parameter_region,)
        assert face.parameter_region.geometric_dimension == 1


def test_a_faces_chart_puts_the_held_coordinate_back(rectangle):
    """The face where x = 2 is the points (2, y), and its chart says so."""
    face = rectangle.face(axis=0, upper=True)
    region = face.parameter_region

    assert region.bounds == ((0.0, 3.0),)
    assert face.held_at == 2.0

    coordinates = SpatialCoordinate(face, point_in(region, 1.5)).expand_geometry()

    assert [coordinates.component(i).as_float() for i in range(2)] == pytest.approx([2.0, 1.5])


def test_the_measure_of_a_face_is_one_its_region_carrying_the_length(rectangle):
    """A face's chart is an inclusion, which stretches nothing.

    With a reference cell the parameter region is always the unit one, so
    the inclusion's Jacobian carries the facet's length. With a region
    the extent is in the region, so the chart is unit and the length
    comes from integrating over it. Both are right, and it is worth
    knowing which convention is in play.
    """
    for face in rectangle.faces:
        point = point_in(face.parameter_region, 1.0)

        assert VolumeElement(face, point).expand_geometry().as_float() == pytest.approx(1.0)
        assert Jacobian(face, point).value_shape == (2, 1)


def test_a_region_states_its_outward_normals_with_no_convention(rectangle, brick):
    """A face knows which way is out of the box from its own description.

    No vertex ordering, no cross product, no orientation left open: the
    held coordinate increases out of the box at its upper end and into
    it at its lower.
    """
    assert [face.outward_normal for face in rectangle.faces] == [
        (-1.0, 0.0),
        (1.0, 0.0),
        (0.0, -1.0),
        (0.0, 1.0),
    ]
    assert len(brick.faces) == 6
    for face in brick.faces:
        assert sum(abs(c) for c in face.outward_normal) == pytest.approx(1.0)


@pytest.mark.parametrize("gdim", [2, 3])
def test_opposite_faces_are_indistinguishable_to_the_normal(gdim):
    """The outward normal is not a function of what UnitNormal is given.

    A face's chart is an inclusion whose offset holds one coordinate
    fixed. A Jacobian does not see an offset, so the lower and upper
    faces of one axis have the same Jacobian at every point -- while
    their outward normals are opposite. No function of the Jacobian can
    tell them apart, so no amount of work inside UnitNormal can produce
    an outward normal.

    This is not about cells. The reference triangle's facets made the
    same point less sharply, where one facet of three came out inward
    and it looked like a question of ordering its vertices. Here two
    faces are given identical data and must give opposite answers.

    What it needs is for the domain to say, as BoxFace.outward_normal
    does, or for the chart to be oriented so that the cross product
    comes out the wanted way. Either way the information belongs to the
    domain and not to the Jacobian.
    """
    box = Box(tuple((0.0, 1.0 + axis) for axis in range(gdim)))

    for axis in range(gdim):
        lower, upper = box.face(axis, False), box.face(axis, True)
        region = lower.parameter_region
        point = point_in(region, *[0.5] * region.geometric_dimension)

        assert lower.parameter_region == upper.parameter_region
        assert entries_of(Jacobian(lower, point), gdim, gdim - 1) == pytest.approx(
            entries_of(Jacobian(upper, point), gdim, gdim - 1)
        )
        assert normal_of(lower, point) == pytest.approx(normal_of(upper, point))

        assert lower.outward_normal == pytest.approx([-c for c in upper.outward_normal])
        agreements = [
            sum(a * b for a, b in zip(normal_of(face, point), face.outward_normal))
            for face in (lower, upper)
        ]
        assert sorted(agreements) == pytest.approx([-1.0, 1.0])


@pytest.mark.parametrize("axis", [0, 1])
def test_a_faces_normal_is_still_normal_to_it(rectangle, axis):
    """What the geometry does give is the normal direction, which is correct.

    The sign is the only thing missing, so a face's normal is plus or
    minus the outward one and never anything else.
    """
    face = rectangle.face(axis, upper=False)
    point = point_in(face.parameter_region, 1.0)

    normal = normal_of(face, point)
    jacobian = _as_dense_matrix(Jacobian(face, point).expand_geometry())
    tangent = [jacobian.component(i, 0).as_float() for i in range(2)]

    assert sum(c * c for c in normal) == pytest.approx(1.0)
    assert sum(a * b for a, b in zip(normal, tangent)) == pytest.approx(0.0)
    assert abs(sum(a * b for a, b in zip(normal, face.outward_normal))) == pytest.approx(1.0)


def test_a_faces_projector_is_the_identity_less_its_normal(brick):
    """A face of a brick is a plane, and the projector onto it is I - n (x) n.

    Which the stated outward normal and the computed one agree on, the
    projector being quadratic in the normal and so blind to its sign.
    """
    face = brick.face(axis=2, upper=True)
    point = point_in(face.parameter_region, 0.5, 1.0)

    projector = entries_of(TangentialProjector(face, point), 3, 3)
    stated = face.outward_normal
    expected = [
        (1.0 if i == j else 0.0) - stated[i] * stated[j] for i in range(3) for j in range(3)
    ]

    assert projector == pytest.approx(expected)
