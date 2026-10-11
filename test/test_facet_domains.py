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
from conftest import Point as PointEntity

# The same analytic surface the surface measure tests use, rather than a
# second copy of a map differentiated by hand, so the two stay in step.
from test_parametrizations import Paraboloidal

from uflx import (
    Coefficient,
    TestFunction,
    composed_domain,
    dx,
    function_space,
    grad,
    inner,
    parametrized_domain,
)
from uflx.algorithms import pull_back_to_entity, replace, simplify
from uflx.domains import (
    RD,
    EntityDomain,
    IdentityParametrization,
    LabelledDomain,
    entity_domain,
    union_domain,
)
from uflx.entities import AbstractEntity
from uflx.expressions import AbstractExpression, RealScalar
from uflx.geometry import (
    AbstractGeometricQuantity,
    ImagePoint,
    Jacobian,
    JacobianInverse,
    JacobianInverseTranspose,
    MetricTensor,
    SingleSpatialCoordinate,
    SpatialCoordinate,
    TangentialProjector,
    UnitNormal,
    VolumeElement,
    _as_dense_matrix,
    expand_geometry,
)
from uflx.graphs import as_graph
from uflx.integrals import Integral, IntegralSum
from uflx.parametrizations import (
    AffineParametrization,
    ComposedDomain,
    ComposedParametrization,
    SubEntityInclusion,
    included_sub_entity,
    sub_entities_of,
    sub_entity_domain,
)
from uflx.points import Point
from uflx.restrictions import Restriction, restrict


class Placed(Entity):
    """An entity that says where its vertices are.

    The reference geometry an element library owns. With it the inclusion
    of a sub-entity is the affine map taking the sub-entity's own
    coordinates to the vertices it has in this cell, which is what
    ``sub_entity_inclusion`` is asked for.
    """

    vertices: tuple[tuple[float, ...], ...]

    def sub_entity_inclusion(self, dim: int, index: int) -> SubEntityInclusion:
        """Get the affine map onto the vertices this cell gives the sub-entity.

        It says which sub-entity it includes, which the cell knows because
        it was asked for that one, and which is what identifies a kernel
        along with the cell type.

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
        return SubEntityInclusion(
            entity_domain(sub_entity), entity_domain(self), matrix, origin, dim, index
        )


class PlacedInterval(Placed, Interval):
    """The reference interval, whose facets are its two end points."""

    vertices = ((0.0,), (1.0,))


class PlacedTriangle(Placed, Triangle):
    """The reference triangle, with its vertices where Basix puts them."""

    vertices = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))


class PlacedTetrahedron(Placed, Tetrahedron):
    """The reference tetrahedron, with its vertices where Basix puts them."""

    vertices = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))


class NestedTetrahedron(PlacedTetrahedron):
    """A tetrahedron whose faces place their own vertices as well.

    An element library does this at every level -- Basix gives reference
    coordinates for a cell, its faces and its edges alike -- so a face of
    a cell is a cell with facets of its own, and the boundary of a
    boundary is a domain. The entities above stop at the cell, which is
    why this is here.
    """

    def sub_entities(self, dim: int) -> list[AbstractEntity]:
        """Get the sub-entities, with the faces placed."""
        if dim == 2:
            return [PlacedTriangle() for _ in super().sub_entities(dim)]
        return super().sub_entities(dim)


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


def entries_of(quantity, rows: int, cols: int) -> list[list[float]]:
    """A matrix valued quantity, as numbers."""
    expanded = expand_geometry(quantity)
    assert isinstance(expanded, AbstractExpression)
    matrix = _as_dense_matrix(expanded)
    return [[matrix.component(i, j).as_float() for j in range(cols)] for i in range(rows)]


def normalised(vector: list[float]) -> list[float]:
    """A vector of numbers, scaled to unit length."""
    length = math.sqrt(sum(c * c for c in vector))
    return [c / length for c in vector]


def transported_normal(cell, domain, index, point) -> list[float]:
    """The reference facet's normal carried onto a domain by ``J^-T``.

    INTEGRATION.md section 4's transport rule, written in the language:
    the reference facet is a domain of codimension one in the cell's
    coordinates, so it has a normal with no reference geometry needed
    beyond the inclusion, and the cell's inverse Jacobian at the included
    point carries it.
    """
    reference_facet = sub_entity_domain(cell, cell.topological_dimension - 1, index)
    nhat = normal_of(reference_facet, point)
    # Written out, since a parametrization whose Jacobian depends on its point reads
    # the coordinates and an image point only has them once expanded.
    included = ImagePoint(point, reference_facet).expand_geometry()
    assert isinstance(included, Point)
    jacobian_inverse_transpose = entries_of(
        JacobianInverseTranspose(domain, included), domain.geometric_dimension, len(nhat)
    )
    return normalised(
        [sum(row[j] * nhat[j] for j in range(len(nhat))) for row in jacobian_inverse_transpose]
    )


def test_a_sub_entity_of_a_cell_is_a_domain_in_the_cells_coordinates(facets):
    """One dimensional, named by the two coordinates of the triangle it sits in."""
    for facet in facets:
        assert facet.topological_dimension == 1
        assert facet.geometric_dimension == 2
        assert facet.cell_types == (Interval(),)
        assert facet.coordinate_domain == EntityDomain(PlacedTriangle())
        assert facet.inner_domain == EntityDomain(Interval())


def test_a_facet_of_a_reference_cell_is_one_domain_however_it_is_reached():
    """Mapping a domain on by its own identity gives that domain back.

    An entity domain's map is the identity, so the piece an enumeration
    builds from it is the sub-entity domain itself and not a composite of
    it with nothing. Domain identity is structural, and it is what a
    restriction and the split onto a piece key on, so two spellings of
    one set had better be one object.

    Only where the identity is of the coordinates the domain's points are
    named in: an identity of R^2 relabels which coordinates a facet of a
    cell is named in, which is a change and stays a composite.
    """
    cell = PlacedTriangle()

    facets = sub_entities_of(entity_domain(cell), 1)

    for k, piece in enumerate(facets.pieces):
        assert piece == sub_entity_domain(cell, 1, k)
        assert hash(piece) == hash(sub_entity_domain(cell, 1, k))

    relabelled = composed_domain(sub_entity_domain(cell, 1, 0), IdentityParametrization(RD(2)))
    assert isinstance(relabelled, ComposedDomain)
    assert relabelled.coordinate_domain == RD(2)


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
    has vertices 1 and 2, and only where those are is missing. The domain
    is the inclusion carrying a sub-entity into the cell, so there is no
    such domain at all until the inclusion is supplied.
    """
    assert Triangle().sub_entities(1)[0] == Interval()
    assert Triangle().sub_entity_vertices(1)[0] == [0, 1]

    with pytest.raises(NotImplementedError, match="does not say where its sub-entities are"):
        sub_entity_domain(Triangle(), 1, 0)


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
    (parametrization,) = facet.parametrizations
    jacobian = _as_dense_matrix(parametrization.jacobian(on_the_facet))
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
    assert image.is_on_an_entity
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
    facet = sub_entity_domain(PlacedTetrahedron(), 2, 3)
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
        AffineParametrization(RD(2), RD(2), ((2.0, 0.0), (0.0, 1.0)), (0.0, 0.0)),
    )


def test_the_facets_of_a_mesh_are_its_cells_facets_mapped_on(mesh):
    """One piece per cell type and local facet, and each is a composition.

    The inclusion of the reference facet into the cell followed by the
    cell's own map, which is the whole of what an exterior facet integral
    is over.
    """
    cell = PlacedTriangle()

    facets = sub_entities_of(mesh, 1)

    assert len(facets.pieces) == 3
    assert facets.geometric_dimension == 2
    assert facets.topological_dimension == 1
    for k, piece in enumerate(facets.pieces):
        assert isinstance(piece, ComposedDomain)
        assert piece.inner_domain == sub_entity_domain(cell, 1, k)
        assert piece.inner_domain.coordinate_domain == EntityDomain(cell)
        assert piece.map == mesh.parametrization(cell)
        assert piece.cell_types == (Interval(),)
        assert piece.parametrizations == (
            ComposedParametrization(cell.sub_entity_inclusion(1, k), mesh.parametrization(cell)),
        )


def test_the_facets_of_a_mesh_of_two_cell_types_are_both_cells_facets(mesh):
    """The fan out is over cell types and local facets together."""
    elements = (
        LagrangeElement(PlacedTriangle(), 1, (3,)),
        LagrangeElement(PlacedTetrahedron(), 1, (3,)),
    )
    mixed = parametrized_domain(elements)

    facets = sub_entities_of(mixed, 1)

    assert len(facets.pieces) == 3 + 4
    # A triangle's facets are one dimensional and a tetrahedron's are two.
    assert facets.topological_dimension is None
    assert facets.geometric_dimension == 3


def test_the_meshs_own_coordinate_in_a_facet_integral_is_rejected(mesh, space):
    """Which is what a boundary traction looks like written without a restriction.

    `x` here is the mesh's, and a point of the facets is not a point of
    the mesh, so reading it there needs the map that carries one to the
    other. It used to build, split and pull back, and fail only when
    something tried to expand it. The facets' own coordinate is a
    different expression, and is read on them.
    """
    facets = sub_entities_of(mesh, 1)
    v = TestFunction(space)
    x = SpatialCoordinate(mesh)

    with pytest.raises(ValueError, match="Cannot read a quantity of"):
        inner(x, x) * dx(facets)

    on_the_facets = SpatialCoordinate(facets)
    assert isinstance(inner(on_the_facets, on_the_facets) * dx(facets), Integral)
    assert isinstance(restrict(x.component(0) * v, facets) * dx(facets), Integral)


def test_a_label_marks_each_piece_so_a_kernel_keeps_it(mesh, space):
    """The consumer says which facets it means, and a kernel still says it.

    The label is on the pieces and not on the union, so splitting the
    integral keeps it: the pull back leaves the piece named in the
    quantities it reads, and the label is there to be read off.
    """
    v = TestFunction(space)
    marked = sub_entities_of(mesh, 1, "inflow")

    assert marked != sub_entities_of(mesh, 1)
    assert marked != sub_entities_of(mesh, 1, "outflow")
    assert all(
        isinstance(piece, LabelledDomain) and piece.label == "inflow" for piece in marked.pieces
    )
    assert all(not isinstance(piece, LabelledDomain) for piece in sub_entities_of(mesh, 1).pieces)

    pulled = pull_back_to_entity(restrict(v, marked) * dx(marked))
    assert isinstance(pulled, IntegralSum)
    for term in pulled.terms:
        read_on = quantities_in(term.integrand)["VolumeElement"].domain
        assert isinstance(read_on, LabelledDomain)
        assert read_on.label == "inflow"


def test_the_sub_entities_of_a_marked_domain_are_refused(mesh):
    """A label on the cells is not a label on their facets, and does not become one.

    Everything geometric about a marked domain is the domain's, the
    parametrization included, so a piece built from it is the unlabelled
    domain's piece: equal to it, hash for hash. The selection has nowhere
    to go, and the boundary of a marked region is a different question
    from a marked selection of facets, so this refuses rather than
    quietly giving the facets of the whole mesh.
    """
    marked = LabelledDomain(mesh, "steel")

    assert marked.parametrization(PlacedTriangle()) == mesh.parametrization(PlacedTriangle())

    with pytest.raises(ValueError, match="selection of cells"):
        sub_entities_of(marked, 1)
    with pytest.raises(ValueError, match="selection of cells"):
        sub_entities_of(marked, 1, "boundary")


def test_the_edges_of_a_cells_faces_are_the_cells_edges():
    """The boundary of a boundary, measured, agrees with taking the codimension.

    Reaching an edge through a face composes two inclusions where taking
    codimension two composes one, so the parametrizations differ and the geometry
    had better not: these are the same six edges of the reference
    tetrahedron, each face contributing three and each edge shared by two
    faces. Only the identity differs, which the test below is about.
    """
    cell = NestedTetrahedron()
    at_the_end = Point([RealScalar(0.3)], entity_domain(Interval()))

    through_the_faces = []
    for face in sub_entities_of(entity_domain(cell), 1).pieces:
        assert isinstance(face, ComposedDomain)
        for edge in sub_entities_of(face, 1).pieces:
            through_the_faces.append(measure_of(edge, at_the_end))

    directly = [
        measure_of(edge, at_the_end) for edge in sub_entities_of(entity_domain(cell), 2).pieces
    ]

    assert sorted(directly) == pytest.approx(sorted(edge_lengths))
    assert len(through_the_faces) == 2 * len(directly)
    assert sorted(set(round(length, 12) for length in through_the_faces)) == pytest.approx(
        sorted(set(round(length, 12) for length in directly))
    )


def test_a_sub_entity_of_a_sub_entity_names_no_index_of_the_cell():
    """An edge of a face is included twice, and neither index is the cell's.

    The geometry is right either way: the edges of a tetrahedron's faces
    are its edges, measured the same. What a parametrization built this way cannot
    say is which edge of the cell it is, since the inner index counts the
    face's edges and nothing here turns the pair into the cell's
    numbering. Two different edges would answer the same, so this refuses
    and says to take the codimension of the cell instead.
    """
    cell = NestedTetrahedron()
    mesh = parametrized_domain(LagrangeElement(cell, 1, (3,)))
    face = sub_entities_of(mesh, 1).pieces[1]
    assert isinstance(face, ComposedDomain)

    assert included_sub_entity(face) == (2, 1)
    assert [included_sub_entity(piece) for piece in sub_entities_of(mesh, 2).pieces] == [
        (1, k) for k in range(6)
    ]

    for edge in sub_entities_of(face, 1).pieces:
        with pytest.raises(ValueError, match="sub-entity of a sub-entity"):
            included_sub_entity(edge)


def test_a_kernel_says_which_sub_entity_it_is_for(mesh, space):
    """The local index, which identifies a kernel along with the cell type.

    The inclusion an entity supplies says which sub-entity it is, so the
    index can be read off a pulled back integral rather than recovered by
    matching the inclusion's matrix against each of the cell's.
    """
    v = TestFunction(space)
    facets = sub_entities_of(mesh, 1)

    assert [included_sub_entity(piece) for piece in facets.pieces] == [(1, 0), (1, 1), (1, 2)]

    pulled = pull_back_to_entity(restrict(v, facets) * dx(facets))
    assert isinstance(pulled, IntegralSum)
    indices = [
        included_sub_entity(quantities_in(term.integrand)["VolumeElement"].domain)
        for term in pulled.terms
    ]
    assert indices == [(1, 0), (1, 1), (1, 2)]

    assert included_sub_entity(mesh) is None


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
    assert sub_entities_of(mesh, 1).topological_dimension == 2


@pytest.mark.parametrize(("k", "length"), [(0, 2.0), (1, 1.0), (2, math.sqrt(5.0))])
def test_a_facet_of_a_mapped_cell_measures_the_mapped_length(
    stretched_cell, on_the_facet, k, length
):
    """Doubling the first axis of the reference triangle, facet by facet.

    Facet 0 runs along that axis and doubles, facet 1 is across it and is
    unchanged, and the hypotenuse becomes the diagonal of a 2 by 1 box.
    The measure is the one measure: nothing about it is a facet measure.
    """
    facets = sub_entities_of(stretched_cell, 1)

    assert measure_of(facets.pieces[k], on_the_facet) == pytest.approx(length)


@pytest.mark.parametrize(("k", "length"), [(0, 2.0), (1, 1.0), (2, math.sqrt(5.0))])
def test_a_pulled_back_facet_integral_expands_to_that_facets_density(
    stretched_cell, on_the_facet, k, length
):
    """The measure an exterior facet integral picks up, as an expression and then a number.

    Pulling back multiplies the piece's volume element into the
    integrand, and expanding writes it out in the coordinates of the
    variable the integral hands over: a facet's parametrization is the inclusion
    followed by the cell's map, so the inclusion is asked for the point
    it passes on. Substituting a quadrature point then gives the length
    the facet is, which is what `measure_of` reads off the piece
    directly.
    """
    facets = sub_entities_of(stretched_cell, 1)

    pulled = pull_back_to_entity(RealScalar(1.0) * dx(facets))

    assert isinstance(pulled, IntegralSum)
    term = pulled.terms[k]
    expanded = expand_geometry(term.integrand)
    substituted = simplify(replace(expanded, {term.variable: on_the_facet}))

    assert isinstance(substituted, AbstractExpression)
    assert substituted.as_float() == pytest.approx(length)
    assert measure_of(facets.pieces[k], on_the_facet) == pytest.approx(length)


def test_an_integral_over_the_facets_measures_the_boundary(stretched_cell, on_the_facet):
    """So the perimeter is the integral of one over the facets, end to end.

    The density is constant on each facet here, the cell's map being
    affine, so the sum of the terms is the boundary's measure with no
    quadrature rule needed.
    """
    pulled = pull_back_to_entity(RealScalar(1.0) * dx(sub_entities_of(stretched_cell, 1)))

    assert isinstance(pulled, IntegralSum)
    densities = []
    for term in pulled.terms:
        substituted = simplify(
            replace(expand_geometry(term.integrand), {term.variable: on_the_facet})
        )
        assert isinstance(substituted, AbstractExpression)
        densities.append(substituted.as_float())

    assert sum(densities) == pytest.approx(3.0 + math.sqrt(5.0))


def test_the_variable_of_an_unsplit_facet_integral_names_no_coordinates(space, mesh):
    """A point of a union is a point of one of its pieces, and which is not said.

    So the coordinates a parametrization would be evaluated at are not there to be
    had until the integral has been split, and the variable says so.
    :mod:`test_variables` does the same for a point of a mesh.
    """
    facets = sub_entities_of(mesh, 1)
    integral = restrict(TestFunction(space), facets) * dx(facets)

    assert isinstance(integral, Integral)
    with pytest.raises(ValueError, match="is not a coordinate tuple"):
        integral.variable.component(0)


@pytest.mark.parametrize("k", range(3))
def test_a_facets_normal_is_its_reference_normal_carried_by_the_cofactor_rule(
    stretched_cell, on_the_facet, k
):
    """INTEGRATION.md section 4's transport rule, which costs no new information.

    A reference facet is a domain of codimension one in the cell's
    coordinates, so it has a normal already, and ``J^-T`` at the included
    point carries it to the mapped facet's. Up to sign, which is the
    thing still missing: outwardness is about the interior and the sign
    here is the parametrization's.

    What it buys is an outward normal for free once a reference normal is
    oriented, since the transport is exact rather than approximate.
    """
    facet = sub_entities_of(stretched_cell, 1).pieces[k]

    transported = transported_normal(PlacedTriangle(), stretched_cell, k, on_the_facet)
    own = normal_of(facet, on_the_facet)

    assert transported == pytest.approx(own) or transported == pytest.approx([-c for c in own])


def test_the_conormal_of_a_facet_of_a_surface_is_the_same_transport(patch, on_the_facet):
    """Where a normal does not exist at all, the transport still gives the right vector.

    A facet of a curved patch is a curve in space, codimension two, so it
    has a normal plane and ``UnitNormal`` refuses. What a surface form
    wants there is the conormal: the direction in the surface, across the
    facet. ``J^-T nhat`` is it, with no cross product and no new
    quantity -- it lies in the tangent plane and is orthogonal to the
    facet, which is the whole of what a conormal is.
    """
    cell = PlacedTriangle()

    for k in range(3):
        facet = sub_entities_of(patch, 1).pieces[k]
        with pytest.raises(ValueError, match="codimension one"):
            UnitNormal(facet, on_the_facet).value_shape

        conormal = transported_normal(cell, patch, k, on_the_facet)
        tangent = [row[0] for row in entries_of(Jacobian(facet, on_the_facet), 3, 1)]
        included = ImagePoint(on_the_facet, sub_entity_domain(cell, 1, k)).expand_geometry()
        assert isinstance(included, Point)
        projector = entries_of(TangentialProjector(patch, included), 3, 3)
        projected = [sum(projector[i][j] * conormal[j] for j in range(3)) for i in range(3)]

        assert sum(a * b for a, b in zip(conormal, tangent, strict=True)) == pytest.approx(
            0.0, abs=1e-12
        )
        assert projected == pytest.approx(conormal)


@pytest.mark.parametrize("k", range(3))
def test_a_facets_density_is_nansons_formula(stretched_cell, on_the_facet, k):
    """``ds = abs(det J) abs(J^-T nhat) dshat``, exactly and by construction.

    The facet's volume element is never computed as a facet measure. It
    is ``sqrt(det g)`` of the facet's own parametrization, which is the inclusion
    followed by the cell's map -- and that comes out as the cell's volume
    scaling times how much the cofactor rule stretches the reference
    normal, which is what a continuum mechanics text calls Nanson's
    formula.
    """
    cell = PlacedTriangle()
    facet = sub_entities_of(stretched_cell, 1).pieces[k]
    reference_facet = sub_entity_domain(cell, 1, k)
    included = ImagePoint(on_the_facet, reference_facet).expand_geometry()
    assert isinstance(included, Point)

    nhat = normal_of(reference_facet, on_the_facet)
    cofactor = entries_of(JacobianInverseTranspose(stretched_cell, included), 2, 2)
    stretch = math.sqrt(sum(sum(row[j] * nhat[j] for j in range(2)) ** 2 for row in cofactor))

    nanson = (
        measure_of(stretched_cell, included) * stretch * measure_of(reference_facet, on_the_facet)
    )

    assert measure_of(facet, on_the_facet) == pytest.approx(nanson)


@pytest.mark.parametrize("k", range(3))
def test_a_facets_metric_is_the_cells_metric_pulled_back(stretched_cell, on_the_facet, k):
    """``g_facet = Diota^T g_cell Diota``, which is what composing the parametrizations means.

    A facet's parametrization is the cell's map after the inclusion, so its first
    fundamental form is the cell's pulled back along the inclusion. Not a
    second computation of a metric: the chain rule, read off the one the
    composed parametrization already induces.
    """
    cell = PlacedTriangle()
    facet = sub_entities_of(stretched_cell, 1).pieces[k]
    reference_facet = sub_entity_domain(cell, 1, k)
    included = ImagePoint(on_the_facet, reference_facet).expand_geometry()
    assert isinstance(included, Point)

    inclusion = entries_of(Jacobian(reference_facet, on_the_facet), 2, 1)
    cell_metric = entries_of(MetricTensor(stretched_cell, included), 2, 2)
    pulled_back = sum(
        inclusion[i][0] * cell_metric[i][j] * inclusion[j][0] for i in range(2) for j in range(2)
    )

    assert entries_of(MetricTensor(facet, on_the_facet), 1, 1)[0][0] == pytest.approx(pulled_back)


def test_the_normal_of_a_facet_of_a_mapped_cell_is_not_the_mapped_normal(
    stretched_cell, on_the_facet
):
    """A normal is carried by the cofactor rule, not by the map, as test_boundaries says."""
    hypotenuse = sub_entities_of(stretched_cell, 1).pieces[2]

    normal = normal_of(hypotenuse, on_the_facet)

    assert sum(c * c for c in normal) == pytest.approx(1.0)
    assert abs(normal[0]) == pytest.approx(1.0 / math.sqrt(5.0))
    assert abs(normal[1]) == pytest.approx(2.0 / math.sqrt(5.0))


@pytest.fixture
def space(mesh):
    """A P1 space on the mesh of placed triangles."""
    return function_space(mesh, LagrangeElement(PlacedTriangle(), 1))


@pytest.fixture
def facet(mesh):
    """One facet of that mesh, a piece of the domain its facets make up."""
    piece = sub_entities_of(mesh, 1).pieces[0]
    assert isinstance(piece, ComposedDomain)
    return piece


def quantities_in(expression):
    """The geometric quantities an expression reads, by class name."""
    return {
        type(node).__name__: node
        for node in as_graph(expression).ordered_nodes()
        if isinstance(node, AbstractGeometricQuantity)
    }


def test_a_restriction_is_the_field_it_restricts(space, facet):
    """A pullback of a function leaves the values alone, so only where changes."""
    u = Coefficient(space)

    restricted = restrict(u, facet)

    assert isinstance(restricted, Restriction)
    assert restricted.field == u
    assert restricted.domain == facet
    assert restricted.value_shape == u.value_shape
    assert restricted.successors == {u}
    assert restricted == restrict(u, facet)
    assert hash(restricted) == hash(restrict(u, facet))


def test_a_field_is_restricted_only_to_pieces_of_its_own_domain(space, facet):
    """The parametrization has to be the one the field's own cells are mapped by."""
    elsewhere = parametrized_domain(LagrangeElement(PlacedTriangle(), 2, (2,)))
    u = Coefficient(space)

    with pytest.raises(ValueError, match="is not a sub-entity of a cell of"):
        restrict(u, sub_entities_of(elsewhere, 1).pieces[0])


def test_a_traction_in_the_spatial_coordinate_is_restricted_like_any_field(mesh, space):
    """The ordinary Neumann term, which had no native spelling.

    `x` here is the mesh's, so reading it on the facets is a map like any
    other -- and the restriction asks which domain its field is on, which
    used to be read off the functions alone. An integrand of pure
    geometry has none, so this asserted.
    """
    facets = sub_entities_of(mesh, 1)
    v = TestFunction(space)
    x = SpatialCoordinate(mesh)

    traction = inner(restrict(x, facets), restrict(x, facets)) * restrict(v, facets) * dx(facets)

    assert isinstance(traction, Integral)
    pulled = pull_back_to_entity(traction)
    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == 3
    for term in pulled.terms:
        assert isinstance(expand_geometry(term.integrand), AbstractExpression)


def test_a_field_is_read_through_one_map_and_so_is_never_restricted_twice():
    """Restricting along a composition is what a second restriction would be.

    ``iota_2^*(iota_1^* u) = (iota_1 . iota_2)^* u``, so iterating adds no
    expressiveness, and it costs identity: these twelve parametrizations are the six
    edges of the tetrahedron, and none of them says which edge it is,
    where taking the codimension of the cell gives each once in the cell's
    own numbering. So the mesh's field is refused here and the message
    says what to write.
    """
    cell = NestedTetrahedron()
    mesh = parametrized_domain(LagrangeElement(cell, 1, (3,)))
    u = Coefficient(function_space(mesh, LagrangeElement(cell, 1)))
    face = sub_entities_of(mesh, 1).pieces[1]
    assert isinstance(face, ComposedDomain)
    edge = sub_entities_of(face, 1).pieces[0]

    with pytest.raises(ValueError, match="sub_entities_of"):
        restrict(u, edge)
    with pytest.raises(ValueError, match="sub_entities_of"):
        restrict(restrict(u, face), edge)

    ridges = sub_entities_of(mesh, 2)
    assert len(ridges.pieces) == 6
    assert isinstance(restrict(u, ridges), Restriction)


def test_a_field_living_on_a_facet_is_restricted_to_that_facets_boundary():
    """Which is one map out of the domain the field is on, and not an iteration.

    A multiplier on a facet is a field whose domain is that facet, so its
    own boundary is a sub-entity of its own cells like any other. This is
    what the rule above leaves available, and it is the case that wanted
    it.
    """
    cell = NestedTetrahedron()
    mesh = parametrized_domain(LagrangeElement(cell, 1, (3,)))
    face = sub_entities_of(mesh, 1).pieces[1]
    assert isinstance(face, ComposedDomain)
    (face_cell,) = face.cell_types
    lam = Coefficient(function_space(face, LagrangeElement(face_cell, 1)))
    edge = sub_entities_of(face, 1).pieces[0]

    integral = restrict(lam, edge) * dx(edge)

    assert isinstance(integral, Integral)
    pulled = pull_back_to_entity(integral)
    assert isinstance(pulled, Integral)
    assert isinstance(expand_geometry(pulled.integrand), AbstractExpression)


def test_a_field_is_restricted_along_an_inclusion_and_not_along_any_map(space, mesh):
    """Restriction is reading a field on a sub-entity of its own cells."""
    u = Coefficient(space)

    with pytest.raises(ValueError, match="is not a sub-entity of a cell"):
        restrict(u, mesh)


def test_an_integral_over_a_facet_takes_a_restricted_field(space, facet):
    """A field on the mesh is integrated over a facet by being read there.

    Without the restriction the integral refuses, as it should: a
    function reaches another domain through a map and not by being
    integrated there.
    """
    u = Coefficient(space)
    v = TestFunction(space)

    integral = inner(restrict(u, facet), restrict(v, facet)) * dx(facet)

    assert isinstance(integral, Integral)
    assert integral.domain == facet
    assert integral.variable.domain == facet

    with pytest.raises(ValueError, match="Cannot integrate a function on"):
        inner(u, v) * dx(facet)


def test_pulling_a_facet_integral_back_lands_on_the_reference_facet(space, facet):
    """The measure moves onto the facet's own coordinate domain and leaves its density.

    That density is the facet's volume element, read at the facet's own
    point: the factor an exterior facet integral picks up, from the one
    measure rather than from a measure of its own.
    """
    v = TestFunction(space)

    pulled = pull_back_to_entity(restrict(v, facet) * dx(facet))

    assert isinstance(pulled, Integral)
    assert pulled.domain == EntityDomain(Interval())
    assert pulled.variable.domain == EntityDomain(Interval())
    volume_element = quantities_in(pulled.integrand)["VolumeElement"]
    assert volume_element == VolumeElement(facet, pulled.variable)


def test_a_restricted_fields_basis_is_read_at_the_included_point(space, facet):
    """The cell's basis, at the point the facet's inclusion gives in the cell."""
    v = TestFunction(space)

    pulled = pull_back_to_entity(restrict(v, facet) * dx(facet))
    assert isinstance(pulled, Integral)

    (function,) = [
        node
        for node in as_graph(pulled.integrand).ordered_nodes()
        if isinstance(node, TestFunction)
    ]
    variable = function.variable
    assert variable == ImagePoint(
        pulled.variable.to_entity_coordinates(Interval()), facet.inner_domain
    )
    assert variable is not None
    assert variable.domain == EntityDomain(PlacedTriangle())
    assert function.is_on_an_entity


def test_a_gradient_inside_a_restriction_is_the_cells_gradient(space, facet, mesh):
    """Not the facet's: a restriction is in the integrand precisely so this holds.

    The inverse Jacobian carrying the entity gradient is the mesh's, read
    at the included point, where the facet's own would be the tangential
    one and would lose the derivative across the facet.
    """
    v = TestFunction(space)
    normal = UnitNormal(facet)

    pulled = pull_back_to_entity(inner(restrict(grad(v), facet), normal) * dx(facet))
    assert isinstance(pulled, Integral)

    quantities = quantities_in(pulled.integrand)
    assert quantities["JacobianInverse"] == JacobianInverse(
        mesh, ImagePoint(pulled.variable.to_entity_coordinates(Interval()), facet.inner_domain)
    )
    assert quantities["UnitNormal"] == UnitNormal(facet, pulled.variable)
    assert quantities["VolumeElement"] == VolumeElement(facet, pulled.variable)


def test_a_restriction_to_a_union_has_no_one_inclusion_to_read_through(space, mesh):
    """Several pieces, several inclusions, so the point is not known until the split."""
    v = TestFunction(space)
    facets = sub_entities_of(mesh, 1)

    integral = restrict(v, facets) * dx(facets)
    assert isinstance(integral, Integral)
    (restriction,) = [
        node
        for node in as_graph(integral.integrand).ordered_nodes()
        if isinstance(node, Restriction)
    ]

    assert restriction.variable == integral.variable
    with pytest.raises(ValueError, match="parametrization per piece"):
        restriction.image_point


def test_an_integral_over_the_facets_splits_into_one_per_facet(space, mesh):
    """The second fan out: a piece names a parametrization, as a cell type does for a mesh."""
    v = TestFunction(space)
    facets = sub_entities_of(mesh, 1)

    integral = restrict(v, facets) * dx(facets)
    assert isinstance(integral, Integral)
    split = integral.split_by_piece()

    assert isinstance(split, IntegralSum)
    assert len(split.terms) == 3
    for piece, term in zip(facets.pieces, split.terms, strict=True):
        assert term.domain == piece
        assert term.measure == dx(piece)
        assert term.variable.domain == piece


def test_a_union_of_one_piece_is_still_split(space, mesh):
    """The union is not the piece, and only a piece names the parametrization to read through."""
    v = TestFunction(space)
    one = union_domain(sub_entities_of(mesh, 1).pieces[:1])

    integral = restrict(v, one) * dx(one)
    assert isinstance(integral, Integral)
    split = integral.split_by_piece()

    assert isinstance(split, IntegralSum)
    assert len(split.terms) == 1
    assert split.terms[0].domain == one.pieces[0]


def test_an_integral_over_a_mesh_has_no_pieces_to_split_by(space, mesh):
    """A cell integral is over one domain, not over a union of them."""
    v = TestFunction(space)

    integral = v * dx(mesh)
    assert isinstance(integral, Integral)

    assert integral.split_by_piece() is None


def test_splitting_reads_a_quantity_off_the_piece_it_lands_on(space, mesh):
    """A normal written against the facets is each facet's own normal.

    The union has a parametrization per piece, so a quantity on it says nothing
    until the split says which piece. Where on the piece it is read is a
    later question, answered when the integral is pulled back onto
    coordinates. The coordinate goes the same way, which also says the
    retargeting carries a quantity's other arguments.
    """
    v = TestFunction(space)
    facets = sub_entities_of(mesh, 1)

    integral = (
        inner(SpatialCoordinate(facets), UnitNormal(facets))
        * SpatialCoordinate(facets).component(1)
        * restrict(v, facets)
        * dx(facets)
    )
    assert isinstance(integral, Integral)
    split = integral.split_by_piece()
    assert isinstance(split, IntegralSum)

    for piece, term in zip(facets.pieces, split.terms, strict=True):
        quantities = quantities_in(term.integrand)
        assert quantities["UnitNormal"] == UnitNormal(piece)
        assert quantities["SpatialCoordinate"] == SpatialCoordinate(piece)
        coordinate = quantities["SingleSpatialCoordinate"]
        assert isinstance(coordinate, SingleSpatialCoordinate)
        assert coordinate == SingleSpatialCoordinate(piece, 1)


def test_pulling_back_an_integral_over_the_facets_gives_one_kernel_per_facet(space, mesh):
    """The whole of an exterior facet integral, from the one measure.

    Each term is over a reference facet, carries that facet's volume
    element as its density, reads the normal off that facet, and reads
    the cell's inverse Jacobian at the point the facet's inclusion gives.
    That is what UFL's ds means, with no measure of its own anywhere in
    it.
    """
    v = TestFunction(space)
    u = Coefficient(space)
    facets = sub_entities_of(mesh, 1)

    pulled = pull_back_to_entity(
        inner(restrict(grad(u), facets), UnitNormal(facets)) * restrict(v, facets) * dx(facets)
    )

    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == 3
    for piece, term in zip(facets.pieces, pulled.terms, strict=True):
        assert isinstance(term, Integral)
        assert isinstance(piece, ComposedDomain)
        assert term.domain == EntityDomain(Interval())
        included = ImagePoint(term.variable.to_entity_coordinates(Interval()), piece.inner_domain)
        quantities = quantities_in(term.integrand)
        assert quantities["VolumeElement"] == VolumeElement(piece, term.variable)
        assert quantities["UnitNormal"] == UnitNormal(piece, term.variable)
        assert quantities["JacobianInverse"] == JacobianInverse(mesh, included)


def test_the_facets_of_a_one_dimensional_mesh_are_its_end_points():
    """Ds in one dimension, which is the counting measure over two points.

    A point has nothing to measure, so its Jacobian is a map out of no
    coordinates and the empty Gram determinant is one. Stretching the
    cell does not change that, where it would change a length: the
    density of a zero dimensional domain is one however it is mapped,
    which is what makes the measure a counting one.
    """
    stretched = composed_domain(
        entity_domain(PlacedInterval()),
        AffineParametrization(RD(1), RD(1), ((3.0,),), (0.0,)),
    )
    at_the_end = Point([], entity_domain(PointEntity()))

    ends = sub_entities_of(stretched, 1)

    assert len(ends.pieces) == 2
    assert ends.topological_dimension == 0
    assert all(measure_of(end, at_the_end) == pytest.approx(1.0) for end in ends.pieces)

    pulled = pull_back_to_entity(RealScalar(1.0) * dx(ends))

    assert isinstance(pulled, IntegralSum)
    densities = []
    for term in pulled.terms:
        substituted = simplify(
            replace(expand_geometry(term.integrand), {term.variable: at_the_end})
        )
        assert isinstance(substituted, AbstractExpression)
        densities.append(substituted.as_float())

    # Two end points, each counted once.
    assert sum(densities) == pytest.approx(2.0)


def test_the_facets_of_a_mesh_of_two_cell_types_give_a_kernel_each(space):
    """A piece per cell type and local facet, so the cell type split has nothing left."""
    triangles = LagrangeElement(PlacedTriangle(), 1, (3,))
    tetrahedra = LagrangeElement(PlacedTetrahedron(), 1, (3,))
    mixed = parametrized_domain((triangles, tetrahedra))
    mixed_space = function_space(
        mixed, (LagrangeElement(PlacedTriangle(), 1), LagrangeElement(PlacedTetrahedron(), 1))
    )
    v = TestFunction(mixed_space)
    facets = sub_entities_of(mixed, 1)

    pulled = pull_back_to_entity(restrict(v, facets) * dx(facets))

    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == 3 + 4
    reference_facets = [term.domain for term in pulled.terms]
    assert reference_facets.count(EntityDomain(Interval())) == 3
    assert reference_facets.count(EntityDomain(Triangle())) == 4


# A triangle laid on the paraboloid z = u^2 + v^2, as the worked example of
# an exterior facet integral over a curved domain.
#
# Why the reference triangle is the thing laid on it, rather than a region
# of the plane or a mesh:
#
# - Facets come from entities. sub_entities_of asks each cell type for the
#   inclusion of its sub-entities, and only an AbstractEntity knows it has
#   three edges and which vertices each has. A region of the plane is a
#   parametrized domain the paraboloid would compose with just as happily, and
#   test_regions.py integrates over one, but it has no entities, so there
#   is nothing for sub_entities_of to enumerate.
# - A reference cell's own parametrization is the identity, so composing it with the
#   paraboloid gives a domain whose geometry is closed form end to end and
#   every density below comes out as a number that can be checked. A mesh
#   sums a basis against coordinate dofs that live outside UFLx, so the
#   same code builds and pulls back to the same three kernels but leaves
#   the arithmetic unevaluated.
# - Its three facets differ: two legs of length one and a hypotenuse of
#   root two, each with its own inclusion, so the fan out into one kernel
#   per facet is visible rather than symmetric.
#
# What it is, then, is one curved triangular patch and not a mesh of them --
# which is also all UFLx ever holds, one cell type and a map out of it.


def paraboloid_area_element(u: float, v: float) -> float:
    """sqrt(1 + |grad f|^2) for the graph of f(u, v) = u^2 + v^2."""
    return math.sqrt(1.0 + 4.0 * (u * u + v * v))


# The arc length density of facet k of the triangle, lifted onto the
# paraboloid: the norm of d/ds of s -> phi(inclusion(s)). Facets 0 and 1 are
# the legs from the origin and facet 2 is the hypotenuse.
facet_densities = (
    lambda s: math.sqrt(1.0 + 4.0 * s * s),
    lambda s: math.sqrt(1.0 + 4.0 * s * s),
    lambda s: math.sqrt(2.0 + (4.0 * s - 2.0) ** 2),
)

# The length of the lifted hypotenuse, from the closed form of
# the integral of sqrt(2 + (4s - 2)^2) over the unit interval.
lifted_hypotenuse_length = (
    2.0 * math.sqrt(6.0) + math.log((2.0 + math.sqrt(6.0)) / (math.sqrt(6.0) - 2.0))
) / 4.0


@pytest.fixture
def patch():
    """The reference triangle laid on the paraboloid z = u^2 + v^2."""
    return composed_domain(entity_domain(PlacedTriangle()), Paraboloidal())


def test_a_triangle_laid_on_a_paraboloid_is_a_curved_patch(patch):
    """Two dimensional in three coordinates, with the graph's own geometry.

    Nothing here is paraboloid specific: the area element is sqrt(det g)
    of the composed parametrization and the normal is the cross product of its
    Jacobian's columns, which for a graph come out as the formulae a
    textbook gives.
    """
    point = Point([RealScalar(0.2), RealScalar(0.3)], entity_domain(PlacedTriangle()))

    assert patch.topological_dimension == 2
    assert patch.geometric_dimension == 3
    assert measure_of(patch, point) == pytest.approx(paraboloid_area_element(0.2, 0.3))

    expected = [-0.4, -0.6, 1.0]
    scale = paraboloid_area_element(0.2, 0.3)
    assert normal_of(patch, point) == pytest.approx([c / scale for c in expected])


@pytest.mark.parametrize("k", range(3))
@pytest.mark.parametrize("s", [0.3, 0.75])
def test_a_facet_of_a_curved_patch_is_a_curve_in_space(patch, k, s):
    """Its measure is arc length, and it varies along a facet that curves.

    The flat hypotenuse has length root two everywhere; lifted onto the
    paraboloid its density is sqrt(2 + (4s - 2)^2), which is smallest in
    the middle. That is the one measure over another domain, with no
    arc length formula written anywhere.
    """
    facet = sub_entities_of(patch, 1).pieces[k]
    point = Point([RealScalar(s)], entity_domain(Interval()))

    assert facet.topological_dimension == 1
    assert facet.geometric_dimension == 3
    assert measure_of(facet, point) == pytest.approx(facet_densities[k](s))


def test_a_curved_facet_in_space_has_no_single_normal(patch):
    """A curve in space has a normal plane, not a normal.

    What a surface form wants along such a facet is the conormal in the
    surface, which needs the patch the facet is a facet of. The patch's
    own normal is codimension one and is in the test above.
    """
    facet = sub_entities_of(patch, 1).pieces[2]
    point = Point([RealScalar(0.3)], entity_domain(Interval()))

    with pytest.raises(ValueError, match="single direction only where"):
        UnitNormal(facet, point).value_shape


def test_the_measure_of_a_curved_facet_integrates_to_its_arc_length(patch):
    """The density UFLx gives, integrated, is the length of the lifted curve.

    Simpson's rule on the reference facet, refined, which is what a
    quadrature in a generated kernel would be doing. The convergence is
    the point: a fixed rule is off by its own error and not by a geometry
    mistake, the integrand being a square root rather than a polynomial.
    """
    facet = sub_entities_of(patch, 1).pieces[2]

    def density(s: float) -> float:
        return measure_of(facet, Point([RealScalar(s)], entity_domain(Interval())))

    def simpson(intervals: int) -> float:
        h = 1.0 / intervals
        total = density(0.0) + density(1.0)
        total += sum((4.0 if i % 2 else 2.0) * density(i * h) for i in range(1, intervals))
        return total * h / 3.0

    errors = [abs(simpson(n) - lifted_hypotenuse_length) for n in (4, 16, 64)]

    assert errors[0] < 1.0e-3
    assert errors[1] < errors[0] / 100.0
    assert errors[2] < errors[1] / 100.0
    assert simpson(64) == pytest.approx(lifted_hypotenuse_length, abs=1.0e-8)


def test_a_form_on_the_facets_of_a_curved_patch_pulls_back_per_facet(patch):
    """The whole example, as a form: three kernels over the reference facet.

    Each carries the facet's arc length, read at the facet's own point,
    and the patch's inverse Jacobian, read at the point the inclusion
    gives inside the triangle. The first is what the measure contributes
    and the second is why a restriction belongs in the integrand: folded
    into the measure, the gradient would have become the facet's
    tangential one and the derivative across the facet would be gone.
    """
    space = function_space(patch, LagrangeElement(PlacedTriangle(), 1))
    u, v = Coefficient(space), TestFunction(space)
    facets = sub_entities_of(patch, 1)

    pulled = pull_back_to_entity(
        inner(restrict(grad(u), facets), restrict(grad(v), facets)) * dx(facets)
    )

    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == 3
    for piece, term in zip(facets.pieces, pulled.terms, strict=True):
        assert isinstance(term, Integral)
        assert isinstance(piece, ComposedDomain)
        assert term.domain == EntityDomain(Interval())
        quantities = quantities_in(term.integrand)
        assert quantities["VolumeElement"] == VolumeElement(piece, term.variable)
        assert quantities["JacobianInverse"] == JacobianInverse(
            patch,
            ImagePoint(term.variable.to_entity_coordinates(Interval()), piece.inner_domain),
        )


def test_a_sum_of_a_cell_integral_and_a_facet_one_pulls_back(mesh, space):
    """A cell integral and a facet integral added, which a - L is.

    The two fan out differently -- one over cell types, one over pieces --
    so a sum of them is where both have to reach the terms. It pulls back
    to one kernel over the reference triangle and one per facet over the
    reference interval.
    """
    u, v = Coefficient(space), TestFunction(space)
    facets = sub_entities_of(mesh, 1)

    form = inner(u, v) * dx(mesh) + inner(restrict(u, facets), restrict(v, facets)) * dx(facets)

    pulled = pull_back_to_entity(form)

    assert isinstance(pulled, IntegralSum)
    domains = [term.domain for term in pulled.terms]
    assert domains.count(EntityDomain(PlacedTriangle())) == 1
    assert domains.count(EntityDomain(Interval())) == 3
