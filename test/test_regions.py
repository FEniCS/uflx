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

from __future__ import annotations

from collections.abc import Sequence

import pytest

from uflx import dx
from uflx.algorithms import pull_back_to_entity
from uflx.domains import (
    RD,
    AbstractCellularDomain,
    AbstractChartedDomain,
    AbstractCoordinateDomain,
    AbstractParametrization,
    IdentityParametrization,
)
from uflx.expressions import AbstractExpression, Integer, RealScalar
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
from uflx.integrals import Integral
from uflx.points import Point
from uflx.tensors import Matrix, Vector


class FaceInclusion(AbstractParametrization):
    """The map that puts a face's held coordinate back.

    A point of a face's parameter region has one coordinate fewer than
    the box it bounds. This writes the held value in at its own axis and
    the region's coordinates around it, which is a face's chart.
    """

    def __init__(self, region: AbstractCoordinateDomain, axis: int, held_at: float, dim: int):
        """Initialise.

        Args:
            region: The box with the held coordinate dropped
            axis: Which coordinate is held
            held_at: The value the held coordinate takes
            dim: How many coordinates the box has
        """
        self._region = region
        self._axis = axis
        self._held_at = held_at
        self._dim = dim

    @property
    def source(self) -> AbstractCoordinateDomain:
        """A face's chart starts from the box with the held coordinate dropped."""
        return self._region

    @property
    def target_dimension(self) -> int:
        """It lands in the box's coordinates."""
        return self._dim

    @property
    def _kept(self) -> list[int]:
        """Which coordinate of the box each of the region's coordinates becomes."""
        return [i for i in range(self._dim) if i != self._axis]

    def value(self, point: AbstractVariable) -> AbstractExpression:
        """Write the held value at its axis and the region's coordinates around it."""
        kept = self._kept
        return Vector(
            [
                RealScalar(self._held_at) if i == self._axis else point.component(kept.index(i))
                for i in range(self._dim)
            ]
        )

    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """The held coordinate does not move, and the rest pass straight through."""
        kept = self._kept
        return Matrix(
            [
                [
                    Integer(0) if i == self._axis else Integer(int(kept.index(i) == j))
                    for j in range(len(kept))
                ]
                for i in range(self._dim)
            ]
        )

    @property
    def is_affine(self) -> bool:
        """An inclusion is affine."""
        return True

    def __repr__(self) -> str:
        """Representation."""
        return f"FaceInclusion({self._region!r}, {self._axis}, {self._held_at}, {self._dim})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, FaceInclusion)
            and other._region == self._region
            and other._axis == self._axis
            and other._held_at == self._held_at
            and other._dim == self._dim
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("FaceInclusion", self._region, self._axis, self._held_at, self._dim))


class Box(AbstractChartedDomain, AbstractCoordinateDomain):
    """An axis aligned box of R^d, as a domain in its own right.

    A region whose points are coordinate tuples, charted by itself
    through the identity. There are no cells in it, no entities and no
    elements: a box is named by its extent and that is the whole of its
    description.
    """

    def __init__(self, bounds: Sequence[tuple[float, float]]):
        """Initialise.

        Args:
            bounds: The lower and upper end of each coordinate's range
        """
        self._bounds = tuple((float(lower), float(upper)) for lower, upper in bounds)
        if any(upper <= lower for lower, upper in self._bounds):
            raise ValueError("Each coordinate of a box must run from a lower end to a higher one.")

    @property
    def bounds(self) -> tuple[tuple[float, float], ...]:
        """The lower and upper end of each coordinate's range."""
        return self._bounds

    @property
    def geometric_dimension(self) -> int:
        """One coordinate per pair of bounds."""
        return len(self._bounds)

    @property
    def chart_sources(self) -> tuple[AbstractCoordinateDomain, ...]:
        """A box is charted by itself."""
        return (self,)

    def chart(self, source: AbstractCoordinateDomain) -> AbstractParametrization:
        """A box's own chart is the identity.

        Raises:
            ValueError: If the region is not this box
        """
        if source != self:
            raise ValueError(f"{self!r} is charted by itself, not by {source!r}.")
        return IdentityParametrization(self)

    def face(self, axis: int, upper: bool) -> BoxFace:
        """Get the face where one coordinate is held at an end of its range.

        Args:
            axis: Which coordinate is held
            upper: Whether it is held at the upper end rather than the lower

        Returns:
            That face, as a domain
        """
        return BoxFace(self, axis, upper)

    @property
    def faces(self) -> tuple[BoxFace, ...]:
        """Every face of this box, two per coordinate."""
        return tuple(
            self.face(axis, upper)
            for axis in range(self.geometric_dimension)
            for upper in (False, True)
        )

    def __repr__(self) -> str:
        """Representation."""
        return f"Box({self._bounds})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, Box) and other._bounds == self._bounds

    def __hash__(self) -> int:
        """Hash."""
        return hash(("Box", self._bounds))


class BoxFace(AbstractChartedDomain):
    """One face of a box, a domain of codimension one in the box's coordinates.

    Its parameter region is the box with the held coordinate dropped, and
    its chart is the inclusion that puts that coordinate back.

    A face knows which way is out of the box from its own description:
    the held coordinate decreases out of the box at the lower end and
    increases out of it at the upper, so the outward normal is minus or
    plus that axis. Nothing is inferred from a vertex ordering, and
    nothing needs a convention. What the geometry cannot do with it is
    the subject of the tests.
    """

    def __init__(self, box: Box, axis: int, upper: bool):
        """Initialise.

        Args:
            box: The box this is a face of
            axis: Which coordinate is held at an end of its range
            upper: Whether it is held at the upper end rather than the lower
        """
        if axis < 0 or axis >= box.geometric_dimension:
            raise ValueError(f"{box!r} has no coordinate {axis}.")
        self._box = box
        self._axis = axis
        self._upper = upper

    @property
    def box(self) -> Box:
        """The box this is a face of."""
        return self._box

    @property
    def axis(self) -> int:
        """Which coordinate is held."""
        return self._axis

    @property
    def held_at(self) -> float:
        """The value the held coordinate takes on this face."""
        return self._box.bounds[self._axis][1 if self._upper else 0]

    @property
    def outward_normal(self) -> tuple[float, ...]:
        """The direction out of the box, as ordinary numbers.

        Minus or plus the held axis, which the face's own description
        gives with nothing to compute and no sign left open.
        """
        sign = 1.0 if self._upper else -1.0
        return tuple(sign if i == self._axis else 0.0 for i in range(self._box.geometric_dimension))

    @property
    def geometric_dimension(self) -> int:
        """A face lives in the coordinates of the box it bounds."""
        return self._box.geometric_dimension

    @property
    def topological_dimension(self) -> int:
        """One less than the box, a coordinate being held."""
        return self._box.geometric_dimension - 1

    @property
    def parameter_region(self) -> Box:
        """The box with the held coordinate dropped."""
        bounds = self._box.bounds
        return Box(tuple(b for i, b in enumerate(bounds) if i != self._axis))

    @property
    def chart_sources(self) -> tuple[AbstractCoordinateDomain, ...]:
        """A face is charted by the box with the held coordinate dropped."""
        return (self.parameter_region,)

    def chart(self, source: AbstractCoordinateDomain) -> AbstractParametrization:
        """The inclusion that puts the held coordinate back.

        Raises:
            ValueError: If the region is not this face's parameter region
        """
        region = self.parameter_region
        if source != region:
            raise ValueError(f"{self!r} is charted by {region!r}, not by {source!r}.")
        return FaceInclusion(region, self._axis, self.held_at, self.geometric_dimension)

    def __repr__(self) -> str:
        """Representation."""
        end = "upper" if self._upper else "lower"
        return f"BoxFace({self._box!r}, {self._axis}, {end})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, BoxFace)
            and other._box == self._box
            and other._axis == self._axis
            and other._upper == self._upper
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("BoxFace", self._box, self._axis, self._upper))


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

    This is not about cells. Two faces are given identical data here and
    must give opposite answers, which is sharper than a facet of a
    reference triangle whose normal comes out inward, where the vertex
    ordering is also in play.

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


def test_r_d_is_charted_by_itself(rectangle):
    """R^d is a domain a measure can be put on, and the simplest one.

    Its points are coordinate tuples already, so its chart is the
    identity and its density is one: the measure of R^d is the Lebesgue
    measure. It is the only coordinate domain with no boundary, which is
    the difference between it and a box.
    """
    space = RD(2)
    point = point_in(space, 0.5, 1.5)

    assert space.chart_sources == (space,)
    assert space.chart(space).is_identity
    assert VolumeElement(space, point).expand_geometry().as_float() == pytest.approx(1.0)
    assert entries_of(MetricTensor(space, point), 2, 2) == pytest.approx([1.0, 0.0, 0.0, 1.0])

    with pytest.raises(ValueError, match="is charted by itself"):
        space.chart(rectangle)


@pytest.mark.parametrize("dim", [1, 2, 3])
def test_a_measure_on_r_d_is_the_lebesgue_one(dim):
    """A measure needs a chart and not a cell, and R^d has one."""
    point = point_in(RD(dim), *[0.5] * dim)
    measure = dx(RD(dim))

    assert measure.domain == RD(dim)
    assert measure.density == VolumeElement(RD(dim))
    assert VolumeElement(RD(dim), point).expand_geometry().as_float() == pytest.approx(1.0)


def test_an_integral_over_r_d_builds(rectangle):
    """With the domain stated there is nothing to infer, so an integrand of none works.

    No finite element function can live on R^d -- an element is attached
    to a cell -- so the integrand here has no functions in it at all,
    which is what an integral whose domain is named allows.
    """
    for domain in (RD(2), rectangle, rectangle.face(axis=0, upper=True)):
        integral = RealScalar(1.0) * dx(domain)
        assert isinstance(integral, Integral)

        assert integral.domain == domain
        assert integral.measure.density == VolumeElement(domain)


def test_a_measure_on_a_regions_boundary_is_an_exterior_facet_measure(rectangle):
    """ds, with no cells and no new kind of measure.

    An exterior boundary integral is the one measure over a domain of
    codimension one, which a box's face is. What is still missing is the
    outward normal's sign, not the measure.
    """
    face = rectangle.face(axis=1, upper=True)

    integral = RealScalar(1.0) * dx(face)
    assert isinstance(integral, Integral)

    assert integral.domain == face
    assert integral.domain.topological_dimension == 1
    assert integral.measure != dx(rectangle)


@pytest.mark.parametrize("domain_of", ["plane", "region", "face"])
def test_a_charted_domain_that_is_not_cellular_cannot_be_pulled_back(rectangle, domain_of):
    """Pulling back is still entity shaped, and says so rather than failing oddly.

    Change of variables onto a chart's parameter region is meaningful for
    any charted domain -- for R^d and for a region it is the identity, and
    for a face it is the inclusion -- but the variable it would hand the
    integrand is a cell's, so this is guarded rather than generalised.
    """
    domain = {
        "plane": RD(2),
        "region": rectangle,
        "face": rectangle.face(axis=0, upper=False),
    }[domain_of]
    integral = RealScalar(1.0) * dx(domain)
    assert isinstance(integral, Integral)

    assert integral.split_by_cell_type() is None

    with pytest.raises(ValueError, match="is not made of cells"):
        integral.cellular_domain
    with pytest.raises(ValueError, match="is not made of cells"):
        pull_back_to_entity(integral)
