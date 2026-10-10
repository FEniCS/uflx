"""Test a sub-entity of a cell as a domain, from the inclusion its entity supplies.

:mod:`test_boundaries` writes the inclusions of a reference triangle's
facets out by hand to see what a mesh has to supply. Here the entity
supplies them instead, through ``sub_entity_inclusion``, and the geometry
that follows is the same: the facet lengths come out as ``1, 1, sqrt(2)``
again, from an object built by asking the cell rather than by writing the
map down.

Where the vertices of a reference cell are is not something UFLx holds --
``sub_entity_vertices`` says which vertices a facet has and not where they
are -- so the entities that place their vertices are defined here, in the
module that uses them, as a mesh's element library would define them.

Nothing in this is specific to a facet. The last tests take the
codimension two sub-entities of a tetrahedron, whose measure is an edge
length, and check that the one thing which does not generalise says so:
an edge in space has no single normal direction.
"""

import math

import pytest
from conftest import Entity, Interval, LagrangeElement, Tetrahedron, Triangle

from uflx import composed_domain, parametrized_domain
from uflx.domains import (
    RD,
    EntityDomain,
    SubEntityDomain,
    entity_domain,
    sub_entity_domain,
)
from uflx.expressions import RealScalar
from uflx.geometry import (
    ImagePoint,
    MetricTensor,
    TangentialProjector,
    UnitNormal,
    VolumeElement,
    _as_dense_matrix,
)
from uflx.parametrizations import (
    AffineParametrization,
    ComposedParametrization,
    facets_of,
    sub_entities_of,
)
from uflx.points import Point


class Placed(Entity):
    """An entity that says where its vertices are.

    The reference geometry an element library owns. With it the inclusion
    of a sub-entity is the affine map taking the sub-entity's own
    coordinates to the vertices it has in this cell, which is what
    ``sub_entity_inclusion`` is asked for.
    """

    vertices: tuple[tuple[float, ...], ...]

    def sub_entity_inclusion(self, dim: int, index: int) -> AffineParametrization:
        """Get the affine map onto the vertices this cell gives the sub-entity.

        Raises:
            ValueError: If the sub-entity is not a simplex, which needs its
                own reference coordinates to place rather than only its
                vertices here.
        """
        sub_entity = self.sub_entities(dim)[index]
        points = [self.vertices[v] for v in self.sub_entity_vertices(dim)[index]]
        if len(points) != dim + 1:
            raise ValueError(f"{sub_entity} is not a simplex, so its vertices do not place it.")
        origin = points[0]
        matrix = [[points[j + 1][i] - origin[i] for j in range(dim)] for i in range(len(origin))]
        return AffineParametrization(entity_domain(sub_entity), matrix, origin)


class PlacedTriangle(Placed, Triangle):
    """The reference triangle, with its vertices where Basix puts them."""

    vertices = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))


class PlacedTetrahedron(Placed, Tetrahedron):
    """The reference tetrahedron, with its vertices where Basix puts them."""

    vertices = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))


# Facet k of the reference triangle and its length. Facet k's vertices are
# Triangle().sub_entity_vertices(1)[k], so facet 2 is the hypotenuse.
facet_lengths = (1.0, 1.0, math.sqrt(2.0))

# Edge k of the reference tetrahedron and its length, likewise: the three
# edges from the origin have length one and the three between the unit
# vectors have length root two.
edge_lengths = (1.0, 1.0, 1.0, math.sqrt(2.0), math.sqrt(2.0), math.sqrt(2.0))


@pytest.fixture
def facets():
    """The three facets of the reference triangle, each a domain of its own."""
    return [sub_entity_domain(PlacedTriangle(), 1, k) for k in range(3)]


@pytest.fixture
def on_the_facet():
    """A point of the reference interval, away from either end."""
    return Point([RealScalar(0.3)], entity_domain(Interval()))


def measure_of(domain, point) -> float:
    """A domain's volume element, as a number."""
    return VolumeElement(domain, point).expand_geometry().as_float()


def normal_of(domain, point) -> list[float]:
    """A domain's unit normal, as numbers."""
    (gdim,) = UnitNormal(domain, point).value_shape
    return [UnitNormal(domain, point).component(i).as_float() for i in range(gdim)]


def test_a_sub_entity_of_a_cell_is_a_domain_in_the_cells_coordinates(facets):
    """One dimensional, named by the two coordinates of the triangle it sits in."""
    for facet in facets:
        assert facet.topological_dimension == 1
        assert facet.geometric_dimension == 2
        assert facet.cell_types == (Interval(),)
        assert facet.coordinate_domain == EntityDomain(PlacedTriangle())


def test_a_sub_entity_domain_compares_and_hashes_by_value():
    """It is a description of which sub-entity of which cell, and nothing more."""
    cell = PlacedTriangle()

    assert sub_entity_domain(cell, 1, 0) == sub_entity_domain(cell, 1, 0)
    assert hash(sub_entity_domain(cell, 1, 0)) == hash(sub_entity_domain(cell, 1, 0))
    assert sub_entity_domain(cell, 1, 0) != sub_entity_domain(cell, 1, 1)
    assert sub_entity_domain(cell, 1, 0) != sub_entity_domain(cell, 0, 0)


def test_a_cell_is_asked_for_a_sub_entity_it_has():
    """The dimension and the index are checked against the cell's own count."""
    cell = PlacedTriangle()

    with pytest.raises(ValueError, match="no sub-entities of dimension 3"):
        sub_entity_domain(cell, 3, 0)
    with pytest.raises(ValueError, match="no sub-entity 3 of dimension 1"):
        sub_entity_domain(cell, 1, 3)


def test_an_entity_that_does_not_place_its_vertices_says_so():
    """Reference geometry is the consumer's, so an entity without it refuses.

    The combinatorics are still there: an unplaced triangle knows facet 0
    has vertices 1 and 2, and only where those are is missing.
    """
    unplaced = sub_entity_domain(Triangle(), 1, 0)

    assert unplaced.sub_entity == Interval()
    with pytest.raises(NotImplementedError, match="does not say where its sub-entities are"):
        unplaced.parametrization(Interval())


@pytest.mark.parametrize("k", range(3))
def test_the_measure_of_a_facet_is_its_length(facets, on_the_facet, k):
    """The same numbers test_boundaries gets, from the inclusion the entity gave."""
    facet = facets[k]

    assert measure_of(facet, on_the_facet) == pytest.approx(facet_lengths[k])
    assert MetricTensor(facet, on_the_facet).expand_geometry().component(0, 0).as_float() == (
        pytest.approx(facet_lengths[k] ** 2)
    )


@pytest.mark.parametrize("k", range(3))
def test_a_facet_of_a_triangle_has_a_unit_normal_across_it(facets, on_the_facet, k):
    """Codimension one fixes the normal up to a sign, which nothing here fixes."""
    facet = facets[k]
    normal = normal_of(facet, on_the_facet)
    jacobian = _as_dense_matrix(facet.parametrization(Interval()).jacobian(on_the_facet))
    tangent = [jacobian.component(i, 0).as_float() for i in range(2)]

    assert sum(c * c for c in normal) == pytest.approx(1.0)
    assert sum(a * b for a, b in zip(normal, tangent)) == pytest.approx(0.0)

    projector = _as_dense_matrix(TangentialProjector(facet, on_the_facet).expand_geometry())
    for i in range(2):
        for j in range(2):
            expected = float(i == j) - normal[i] * normal[j]
            assert projector.component(i, j).as_float() == pytest.approx(expected)


@pytest.mark.parametrize("k", range(3))
def test_a_point_of_a_facet_is_included_in_the_cells_coordinates(facets, on_the_facet, k):
    """Three tenths along the facet, in the coordinates of the triangle.

    This is the point a cell's basis is read at on one of its facets, so
    it has to arrive in the cell's coordinates and say that it has.
    """
    first, second = Triangle().sub_entity_vertices(1)[k]
    start, end = PlacedTriangle.vertices[first], PlacedTriangle.vertices[second]
    expected = [start[i] + 0.3 * (end[i] - start[i]) for i in range(2)]

    image = ImagePoint(on_the_facet, facets[k])

    assert image.domain == EntityDomain(PlacedTriangle())
    assert image.in_entity_coordinates
    expanded = image.expand_geometry()
    assert [expanded.component(i).as_float() for i in range(2)] == pytest.approx(expected)


@pytest.mark.parametrize("k", range(6))
def test_the_measure_of_an_edge_of_a_tetrahedron_is_its_length(k):
    """Codimension two, and the same construction: the volume element is sqrt(det g)."""
    edge = sub_entity_domain(PlacedTetrahedron(), 1, k)
    point = Point([RealScalar(0.3)], entity_domain(Interval()))

    assert edge.topological_dimension == 1
    assert edge.geometric_dimension == 3
    assert measure_of(edge, point) == pytest.approx(edge_lengths[k])


def test_an_edge_in_space_has_no_single_normal():
    """What does not generalise above codimension one is a claim, not a construction.

    The directions orthogonal to an edge of a tetrahedron are a plane, so
    there is no normal to return and the quantity says so rather than
    picking one.
    """
    edge = sub_entity_domain(PlacedTetrahedron(), 1, 0)
    point = Point([RealScalar(0.3)], entity_domain(Interval()))

    with pytest.raises(ValueError, match="single direction only where"):
        UnitNormal(edge, point).value_shape


def test_a_facet_of_a_tetrahedron_is_a_triangle_with_a_normal():
    """Codimension one in three dimensions, to check the geometry is not planar only."""
    facet = SubEntityDomain(PlacedTetrahedron(), 2, 3)
    point = Point([RealScalar(0.2), RealScalar(0.3)], entity_domain(Triangle()))

    # Facet 3 is the one between the three unit vectors, of area sqrt(3)/2
    # against the reference triangle's half, so its volume element is sqrt(3).
    assert measure_of(facet, point) == pytest.approx(math.sqrt(3.0))
    assert normal_of(facet, point) == pytest.approx([1.0 / math.sqrt(3.0)] * 3)


@pytest.fixture
def mesh():
    """A mesh of placed triangles in the plane, geometry described by an element."""
    return parametrized_domain(LagrangeElement(PlacedTriangle(), 1, (2,)))


@pytest.fixture
def stretched_cell():
    """The reference triangle with the first axis doubled, in closed form.

    A domain of one cell whose map expands to numbers, which a mesh's
    cannot: a mesh's coordinate dofs are external. So this is what lets a
    facet's measure be checked against the length it should be.
    """
    return composed_domain(
        entity_domain(PlacedTriangle()),
        AffineParametrization(RD(2), ((2.0, 0.0), (0.0, 1.0)), (0.0, 0.0)),
    )


def test_the_facets_of_a_mesh_are_its_cells_facets_mapped_on(mesh):
    """One piece per cell type and local facet, and each is a composition.

    The inclusion of the reference facet into the cell followed by the
    cell's own map, which is the whole of what an exterior facet integral
    is over.
    """
    cell = PlacedTriangle()

    facets = facets_of(mesh)

    assert len(facets.pieces) == 3
    assert facets.geometric_dimension == 2
    assert facets.topological_dimension == 1
    for k, piece in enumerate(facets.pieces):
        assert piece.inner_domain == sub_entity_domain(cell, 1, k)
        assert piece.map == mesh.parametrization(cell)
        assert piece.cell_types == (Interval(),)
        assert piece.parametrization_out_of(entity_domain(Interval())) == ComposedParametrization(
            cell.sub_entity_inclusion(1, k), mesh.parametrization(cell)
        )


def test_the_facets_of_a_mesh_of_two_cell_types_are_both_cells_facets(mesh):
    """The fan out is over cell types and local facets together."""
    elements = (
        LagrangeElement(PlacedTriangle(), 1, (3,)),
        LagrangeElement(PlacedTetrahedron(), 1, (3,)),
    )
    mixed = parametrized_domain(elements)

    facets = facets_of(mixed)

    assert len(facets.pieces) == 3 + 4
    # A triangle's facets are one dimensional and a tetrahedron's are two.
    assert facets.topological_dimension is None
    assert facets.geometric_dimension == 3


def test_a_label_tells_one_part_of_a_boundary_from_another(mesh):
    """The consumer says which facets it means; the measure never carries an id."""
    assert facets_of(mesh) == facets_of(mesh)
    assert facets_of(mesh, 1) != facets_of(mesh, 2)
    assert facets_of(mesh, 1).label == 1
    assert facets_of(mesh).label is None


def test_sub_entities_of_a_codimension_a_cell_does_not_have_are_refused(mesh):
    """A triangle has no sub-entities of codimension three."""
    with pytest.raises(ValueError, match="no sub-entities of codimension 3"):
        sub_entities_of(mesh, 3)


def test_the_codimension_two_sub_entities_of_a_tetrahedral_mesh_are_its_edges():
    """Nothing in the construction is specific to a facet."""
    mesh = parametrized_domain(LagrangeElement(PlacedTetrahedron(), 1, (3,)))

    edges = sub_entities_of(mesh, 2)

    assert len(edges.pieces) == 6
    assert edges.topological_dimension == 1
    assert edges.geometric_dimension == 3
    assert facets_of(mesh).topological_dimension == 2


@pytest.mark.parametrize(("k", "length"), [(0, 2.0), (1, 1.0), (2, math.sqrt(5.0))])
def test_a_facet_of_a_mapped_cell_measures_the_mapped_length(
    stretched_cell, on_the_facet, k, length
):
    """Doubling the first axis of the reference triangle, facet by facet.

    Facet 0 runs along that axis and doubles, facet 1 is across it and is
    unchanged, and the hypotenuse becomes the diagonal of a 2 by 1 box.
    The measure is the one measure: nothing about it is a facet measure.
    """
    facets = facets_of(stretched_cell)

    assert measure_of(facets.pieces[k], on_the_facet) == pytest.approx(length)


def test_the_normal_of_a_facet_of_a_mapped_cell_is_not_the_mapped_normal(
    stretched_cell, on_the_facet
):
    """A normal is carried by the cofactor rule, not by the map, as test_boundaries says."""
    hypotenuse = facets_of(stretched_cell).pieces[2]

    normal = normal_of(hypotenuse, on_the_facet)

    assert sum(c * c for c in normal) == pytest.approx(1.0)
    assert abs(normal[0]) == pytest.approx(1.0 / math.sqrt(5.0))
    assert abs(normal[1]) == pytest.approx(2.0 / math.sqrt(5.0))
