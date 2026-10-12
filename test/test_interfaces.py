"""Test an interface, the domain two others meet along, presented from each side.

UFL's ``dS`` is a measure kind. INTEGRATION.md section 5 says it is the
one measure over a different domain, an interface, with the two sides in
the integrand as pullbacks along each side's embedding. Nothing in that
asks for a cell, so the first tests glue two reference triangles in the
plane with closed form maps and the last glue two boxes; between them is
the case UFL knows, a mesh with itself, which pulls back to one kernel
per pair of facets.

The chart transition is never an object. The ``-`` side's point is the
preimage, through its own map, of the image of the interface's point,
and the tests play the library by substituting the point that preimage
stands for.
"""

import math

import pytest
from conftest import Interval, LagrangeElement, Quadrilateral, Triangle
from conftest import Point as PointEntity
from test_facet_domains import (
    Placed,
    PlacedInterval,
    PlacedTetrahedron,
    PlacedTriangle,
    quantities_in,
)
from test_regions import Box

from uflx import (
    Coefficient,
    TestFunction,
    avg,
    composed_domain,
    dx,
    function_space,
    grad,
    inner,
    interface,
    jump,
    parametrized_domain,
    restrict,
)
from uflx.algorithms import pull_back_to_entity, replace, simplify
from uflx.domains import (
    RD,
    AbstractCoordinateDomain,
    AbstractParametrization,
    EmbedsAnInterface,
    EntityDomain,
    LabelledDomain,
    LabelledParametrization,
    MadeOfPieces,
    entity_domain,
    union_domain,
)
from uflx.expressions import AbstractExpression, RealScalar, expression_sum
from uflx.functions import extract_domain
from uflx.geometry import (
    ImagePoint,
    JacobianInverse,
    OutwardNormal,
    PreimagePoint,
    SpatialCoordinate,
    UnitNormal,
    VolumeElement,
    expand_geometry,
)
from uflx.graphs import GraphNode, as_graph
from uflx.integrals import Integral, IntegralSum
from uflx.interfaces import CoordinateRepresentation, Embedding, Interface
from uflx.parametrizations import (
    AffineParametrization,
    ComposedDomain,
    ComposedParametrization,
    included_sub_entity,
    sub_entities_of,
    sub_entity_domain,
)
from uflx.points import Point
from uflx.restrictions import Restriction
from uflx.tensors import Matrix, Vector
from uflx.variables import AbstractVariable


def as_numbers(expression, length: int) -> list[float]:
    """A vector valued expression, expanded, as numbers."""
    expanded = expand_geometry(expression)
    assert isinstance(expanded, AbstractExpression)
    return [as_number(expanded.component(i)) for i in range(length)]


def as_number(expression) -> float:
    """A scalar expression, simplified, as a number."""
    simplified = simplify(expression)
    assert isinstance(simplified, AbstractExpression)
    return simplified.as_float()


def preimage_in(expression):
    """The one preimage point an expression reads, which stands for a library's answer."""
    (point,) = [node for node in as_graph(expression) if isinstance(node, PreimagePoint)]
    return point


@pytest.fixture
def plus():
    """The reference triangle placed in the plane by the identity, as a domain of R^2."""
    return composed_domain(
        entity_domain(PlacedTriangle()),
        AffineParametrization(RD(2), RD(2), ((1.0, 0.0), (0.0, 1.0)), (0.0, 0.0)),
    )


@pytest.fixture
def minus():
    """The reference triangle reflected through (1/2, 1/2), sharing the hypotenuse."""
    return composed_domain(
        entity_domain(PlacedTriangle()),
        AffineParametrization(RD(2), RD(2), ((-1.0, 0.0), (0.0, -1.0)), (1.0, 1.0)),
    )


def hypotenuse_of(domain):
    """Facet 2 of a placed triangle mapped on, as a piece of its facets."""
    return composed_domain(
        sub_entity_domain(PlacedTriangle(), 1, 2), domain.parametrization(PlacedTriangle())
    )


@pytest.fixture
def glued(plus, minus):
    """The interface the two triangles meet along, both as their facet 2.

    The reflection reverses the facet, so the chart transition is
    ``s -> 1 - s``.
    """
    return interface(hypotenuse_of(plus), (("+", hypotenuse_of(plus)), ("-", hypotenuse_of(minus))))


@pytest.fixture
def on_the_facet():
    """A point of the reference interval, away from either end."""
    return Point([RealScalar(0.3)], entity_domain(Interval()))


def interior_facets_of(mesh, label=None):
    """UFL's dS as a domain: the interface of a mesh's facets with themselves."""
    facets = sub_entities_of(mesh, 1, label)
    return interface(facets, (("+", facets), ("-", facets)))


@pytest.fixture
def mesh():
    """A mesh of placed triangles in the plane."""
    return parametrized_domain(LagrangeElement(PlacedTriangle(), 1, (2,)))


@pytest.fixture
def space(mesh):
    """A P1 space on it."""
    return function_space(mesh, LagrangeElement(PlacedTriangle(), 1))


def test_an_interface_is_the_domain_it_is_over_with_a_presentation_per_side(glued, plus):
    """Its geometry is the domain's, and each side is a domain of its own."""
    assert glued.domain == hypotenuse_of(plus)
    assert glued.labels == ("+", "-")
    assert glued.parametrizations == hypotenuse_of(plus).parametrizations
    assert glued.geometric_dimension == 2
    assert glued.topological_dimension == 1

    for side in glued.embeddings:
        assert isinstance(side, Embedding)
        assert isinstance(side, EmbedsAnInterface)
        assert side.interface == glued
        assert side == glued.embedding(side.label)
        assert side.image == glued.image(side.label)
    assert glued.embedding("+") != glued.embedding("-")

    assert isinstance(glued, MadeOfPieces)
    assert glued.pieces == (glued,)


def test_an_embeddings_parametrization_is_its_images_named_by_the_embedding(glued, minus):
    """The map is the image's. What the side adds is identity."""
    (parametrization,) = glued.embedding("-").parametrizations
    assert isinstance(parametrization, LabelledParametrization)
    assert parametrization.parametrization == hypotenuse_of(minus).parametrizations[0]
    assert parametrization.label == glued.embedding("-")
    assert parametrization != glued.embedding("+").parametrizations[0]
    assert parametrization.is_affine


def test_the_measure_of_an_interface_is_its_domains(glued, on_the_facet):
    """The hypotenuse has length root two, read off the domain's own map."""
    density = VolumeElement(glued, on_the_facet).expand_geometry()
    assert density.as_float() == pytest.approx(math.sqrt(2.0))

    integral = RealScalar(1.0) * dx(glued)
    assert isinstance(integral, Integral)
    assert integral.measure.density == VolumeElement(glued)


def test_the_jump_is_side_symmetric_in_the_measure(plus, minus, on_the_facet):
    """Either side induces the same density, differing by an isometry of the facet."""
    from_plus = interface(
        hypotenuse_of(plus), (("+", hypotenuse_of(plus)), ("-", hypotenuse_of(minus)))
    )
    from_minus = interface(
        hypotenuse_of(minus), (("-", hypotenuse_of(minus)), ("+", hypotenuse_of(plus)))
    )

    assert from_plus != from_minus
    assert VolumeElement(from_plus, on_the_facet).expand_geometry().as_float() == pytest.approx(
        VolumeElement(from_minus, on_the_facet).expand_geometry().as_float()
    )


def test_the_primary_embeddings_transition_is_the_interfaces_own(glued, on_the_facet):
    """The side the domain is the image of is read at the integral's point."""
    assert glued.primary == "+"
    assert glued.embedding("+").transition(on_the_facet) is on_the_facet


def test_the_other_embeddings_transition_is_written_as_a_preimage(glued, on_the_facet):
    """``psi_-^-1(psi_+(X))``: an image through the domain's map, pulled back through the side's."""
    transition = glued.embedding("-").transition(on_the_facet)

    assert isinstance(transition, PreimagePoint)
    assert transition.domain == entity_domain(Interval())
    assert transition.parametrization == glued.embedding("-").parametrizations[0]
    assert transition.ambient_point == ImagePoint(on_the_facet, glued.domain)


def test_a_point_of_another_region_is_none_of_a_sides_business(glued):
    """A point of a mesh's cell is not a point of the interface's region."""
    elsewhere = Point([RealScalar(0.2), RealScalar(0.3)], entity_domain(Triangle()))
    assert glued.embedding("-").transition(elsewhere) is elsewhere


def test_the_transition_is_the_reflection_of_the_facet(glued, minus, on_the_facet):
    """Substituting what the preimage stands for lands on the same point of the plane.

    The library's answer here is ``1 - s``. Carried through the ``-``
    side's map it gives the point the ``+`` side's map gives at ``s``,
    which is what makes the two charts one set.
    """
    transition = glued.embedding("-").transition(on_the_facet)
    answer = Point([RealScalar(1.0) - RealScalar(0.3)], entity_domain(Interval()))

    through_minus = SpatialCoordinate(hypotenuse_of(minus), transition)
    substituted = replace(expand_geometry(through_minus), {transition: answer})
    assert isinstance(substituted, AbstractExpression)
    through_plus = SpatialCoordinate(glued.domain, on_the_facet)

    assert [as_number(substituted.component(i)) for i in range(2)] == pytest.approx(
        as_numbers(through_plus, 2)
    )


def test_the_outward_normals_of_the_two_sides_are_opposite(glued, on_the_facet):
    """Out of each triangle across the shared edge, and the two point against each other.

    Each is the reference facet's outward normal carried by the side's
    map. The reflection does not flip it, where it would flip a cross
    product of the mapped facet's own Jacobian: outwardness is about the
    interior, which INTEGRATION.md section 4 says.
    """
    plus = OutwardNormal(glued.embedding("+")).reconstruct_with_variable(on_the_facet)
    minus = OutwardNormal(glued.embedding("-")).reconstruct_with_variable(on_the_facet)
    assert isinstance(minus.variable, PreimagePoint)

    root_half = 1.0 / math.sqrt(2.0)
    assert as_numbers(plus, 2) == pytest.approx([root_half, root_half])
    assert as_numbers(minus, 2) == pytest.approx([-root_half, -root_half])

    direction = as_numbers(UnitNormal(glued, on_the_facet), 2)
    assert abs(sum(a * b for a, b in zip(direction, as_numbers(plus, 2)))) == pytest.approx(1.0)


def test_an_outward_normal_of_a_reference_facet_is_the_entitys_own():
    """With the cell map the identity, the transport does nothing."""
    from test_boundaries import outward_normals

    for k in range(3):
        facet = sub_entity_domain(PlacedTriangle(), 1, k)
        point = Point([RealScalar(0.3)], entity_domain(Interval()))
        assert as_numbers(OutwardNormal(facet, point), 2) == pytest.approx(outward_normals[k])
        assert PlacedTriangle().facet_outward_normal(k) == pytest.approx(outward_normals[k])


def test_an_outward_normal_of_a_tetrahedrons_face_points_away_from_the_other_vertex():
    """Out of a tetrahedron across each face, which the fourth vertex is inside of."""
    cell = PlacedTetrahedron()
    for k in range(4):
        normal = cell.facet_outward_normal(k)
        on_face = cell.sub_entity_vertices(2)[k]
        (inside,) = [v for v in range(4) if v not in on_face]
        origin = cell.vertices[on_face[0]]
        towards = sum(n * (a - b) for n, a, b in zip(normal, cell.vertices[inside], origin))
        assert towards < 0
        assert sum(c * c for c in normal) == pytest.approx(1.0)


def test_an_outward_normal_wants_a_facet(glued, mesh):
    """Out of a cell across an edge of a tetrahedron is no direction."""
    tetrahedra = parametrized_domain(LagrangeElement(PlacedTetrahedron(), 1, (3,)))
    edge = sub_entities_of(tetrahedra, 2).pieces[0]
    point = Point([RealScalar(0.3)], entity_domain(Interval()))

    with pytest.raises(ValueError, match="is not a facet"):
        OutwardNormal(edge, point).expand_geometry()
    with pytest.raises(ValueError, match="is not a sub-entity of a cell"):
        OutwardNormal(
            mesh, Point([RealScalar(0.2), RealScalar(0.3)], entity_domain(Triangle()))
        ).expand_geometry()


def test_an_interface_refuses_what_is_not_two_presentations_of_one_domain(plus, minus):
    """Fewer than two sides, two of one label, or an image of another dimension."""
    gamma = hypotenuse_of(plus)
    with pytest.raises(ValueError, match="two embeddings"):
        Interface(gamma, (("+", gamma),))
    with pytest.raises(ValueError, match="told apart by label"):
        Interface(gamma, (("+", gamma), ("+", hypotenuse_of(minus))))
    with pytest.raises(ValueError, match="topological dimension"):
        Interface(gamma, (("+", gamma), ("-", minus)))
    with pytest.raises(ValueError, match="has no side"):
        gamma_interface = Interface(gamma, (("+", gamma), ("-", hypotenuse_of(minus))))
        gamma_interface.embedding("top")


def test_a_field_is_restricted_to_a_side_and_not_to_the_interface(space, mesh):
    """Which map a field is read through is the whole question, so it is named."""
    gamma = interior_facets_of(mesh)
    u = Coefficient(space)

    restricted = restrict(u, gamma.embedding("-"))
    assert isinstance(restricted, Restriction)
    assert restricted.domain == gamma.embedding("-")
    assert restricted.reads_on == gamma

    with pytest.raises(ValueError, match="Name it"):
        restrict(u, gamma)


def test_the_jump_is_a_difference_of_restrictions_on_one_domain(space, mesh):
    """``iota_+^* u - iota_-^* u`` typechecks because both are fields on the interface."""
    gamma = interior_facets_of(mesh)
    u = Coefficient(space)

    assert jump(u, gamma) == restrict(u, gamma.embedding("+")) - restrict(u, gamma.embedding("-"))
    assert avg(u, gamma) == RealScalar(0.5) * (
        restrict(u, gamma.embedding("+")) + restrict(u, gamma.embedding("-"))
    )
    assert extract_domain(jump(u, gamma)) == gamma

    integral = jump(u, gamma) * restrict(TestFunction(space), gamma.embedding("+")) * dx(gamma)
    assert isinstance(integral, Integral)
    assert integral.domain == gamma


def test_a_field_reaches_an_interface_only_through_a_side(space, mesh):
    """An unrestricted field in an interface integral is refused, as in a facet one."""
    gamma = interior_facets_of(mesh)
    with pytest.raises(ValueError, match="through a map"):
        Coefficient(space) * dx(gamma)


def test_a_side_of_an_unsplit_interface_names_no_one_map(space, mesh):
    """Three facets a side, so no one transition until the integral is split."""
    gamma = interior_facets_of(mesh)
    point = Point([RealScalar(0.3)], entity_domain(Interval()))
    with pytest.raises(ValueError, match="Split the integral by piece"):
        gamma.embedding("-").transition(point)


def test_the_interface_of_a_mesh_with_itself_is_every_pair_of_facets(mesh):
    """Three facets a side gives nine pairs, each an interface of one facet a side."""
    gamma = interior_facets_of(mesh)
    facets = sub_entities_of(mesh, 1)

    assert gamma.domain == facets
    assert gamma.image("+") == facets
    assert gamma.image("-") == facets

    pairs = gamma.pieces
    assert len(pairs) == 9
    assert len(set(pairs)) == 9
    facet_pairs = []
    for pair in pairs:
        assert isinstance(pair, Interface)
        assert pair.domain == pair.image("+")
        assert pair.primary == "+"
        assert pair.pieces == (pair,)
        facet_pairs.append(
            (
                included_sub_entity(pair.image("+")),
                included_sub_entity(pair.image("-")),
            )
        )
    assert facet_pairs == [((1, k), (1, m)) for k in range(3) for m in range(3)]


def test_pulling_an_interface_integral_back_gives_one_kernel_per_pair(space, mesh):
    """The whole of an interior facet integral, from the one measure.

    Each term is over the reference facet and carries the pair's volume
    element and the ``+`` side's outward normal. The ``+`` restriction
    reads the cell's inverse Jacobian at the included point, as a facet
    integral does; the ``-`` restriction reads it at the point of the
    ``-`` facet that the chart transition gives, carried through the
    ``-`` facet's inclusion. Nine kernels, telling their pair.
    """
    gamma = interior_facets_of(mesh)
    u, v = Coefficient(space), TestFunction(space)

    pulled = pull_back_to_entity(
        inner(jump(grad(u), gamma), OutwardNormal(gamma.embedding("+")))
        * restrict(v, gamma.embedding("+"))
        * dx(gamma)
    )

    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == 9
    for pair, term in zip(gamma.pieces, pulled.terms, strict=True):
        assert isinstance(term, Integral)
        assert isinstance(pair, Interface)
        assert term.domain == EntityDomain(Interval())
        variable = term.variable.to_entity_coordinates(Interval())

        quantities = quantities_in(term.integrand)
        assert quantities["VolumeElement"] == VolumeElement(pair, term.variable)
        assert quantities["OutwardNormal"] == OutwardNormal(pair.embedding("+"), term.variable)

        plus_piece, minus_piece = pair.image("+"), pair.image("-")
        assert isinstance(plus_piece, ComposedDomain)
        assert isinstance(minus_piece, ComposedDomain)
        at_plus = ImagePoint(variable, plus_piece.inner_domain)
        (side_map,) = pair.embedding("-").parametrizations
        at_minus = ImagePoint(
            PreimagePoint(ImagePoint(variable, pair.domain), side_map), minus_piece.inner_domain
        )
        inverses = {node for node in as_graph(term.integrand) if isinstance(node, JacobianInverse)}
        assert inverses == {JacobianInverse(mesh, at_plus), JacobianInverse(mesh, at_minus)}


def test_a_label_survives_to_the_pair(space, mesh):
    """A marked interface, as ``dS(3)``, keeps its mark on every piece of every pair."""
    gamma = interior_facets_of(mesh, "interior-3")
    pulled = pull_back_to_entity(jump(Coefficient(space), gamma) * dx(gamma))

    assert isinstance(pulled, IntegralSum)
    for term in pulled.terms:
        density = quantities_in(term.integrand)["VolumeElement"]
        pair = density.domain
        assert isinstance(pair, Interface)
        for label in pair.labels:
            image = pair.image(label)
            assert isinstance(image, LabelledDomain)
            assert image.label == "interior-3"


def test_the_pairs_of_a_mixed_mesh_match_facet_types():
    """A triangle's edge and a tetrahedron's face are no interface.

    Nine pairs of edges and sixteen of faces, and nothing between.
    """
    mesh = parametrized_domain(
        (LagrangeElement(PlacedTriangle(), 1, (3,)), LagrangeElement(PlacedTetrahedron(), 1, (3,)))
    )
    pairs = interior_facets_of(mesh).pieces
    regions = [{p.source for p in pair.parametrizations} for pair in pairs]

    assert len(pairs) == 9 + 16
    assert all(len(r) == 1 for r in regions)
    assert sum(r == {entity_domain(Interval())} for r in regions) == 9
    assert sum(r == {entity_domain(Triangle())} for r in regions) == 16


def test_two_boxes_meet_along_a_face_with_no_cells_anywhere(on_the_facet):
    """The same object with no reference cell in it.

    Each side is a face of its box, parametrized by that box's own
    inclusion, and each box says which way is out of it. The measure is
    one, as a region's faces are measured by their region.
    """
    left = Box(((0.0, 1.0), (0.0, 1.0)))
    right = Box(((1.0, 2.0), (0.0, 1.0)))
    shared = left.face(axis=0, upper=True)
    gamma = interface(shared, (("left", shared), ("right", right.face(axis=0, upper=False))))
    point = Point([RealScalar(0.5)], shared.parameter_region)

    assert gamma.pieces == (gamma,)
    assert gamma.topological_dimension == 1
    integral = RealScalar(1.0) * dx(gamma)
    assert isinstance(integral, Integral)
    assert VolumeElement(gamma, point).expand_geometry().as_float() == pytest.approx(1.0)

    out_of_left = OutwardNormal(gamma.embedding("left")).reconstruct_with_variable(point)
    out_of_right = OutwardNormal(gamma.embedding("right")).reconstruct_with_variable(point)
    assert as_numbers(out_of_left, 2) == pytest.approx([1.0, 0.0])
    assert as_numbers(out_of_right, 2) == pytest.approx([-1.0, 0.0])

    with pytest.raises(ValueError, match="is not made of cells"):
        pull_back_to_entity(integral)


class PlacedQuadrilateral(Placed, Quadrilateral):
    """The reference square, with its vertices where Basix puts them."""

    vertices = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0))


@pytest.fixture
def along_the_edge():
    """The reference interval laid along the hypotenuse, as a domain of its own.

    Traversed the other way from either triangle's facet, ``s -> (s, 1 - s)``,
    so that neither side's chart is the interface's.
    """
    return composed_domain(
        entity_domain(Interval()),
        AffineParametrization(RD(1), RD(2), ((1.0,), (-1.0,)), (0.0, 1.0)),
    )


def test_an_interval_laid_along_the_edge_is_an_interface_of_its_own(
    plus, minus, along_the_edge, on_the_facet
):
    """Section 5's general case: the interface is a facet of neither side.

    No side is primary, so both are read through a preimage. The library's
    answers here are ``1 - s`` for the ``+`` side and ``s`` for the ``-``,
    and substituting them lands each side on the interface's own point.
    """
    gamma = interface(along_the_edge, (("+", hypotenuse_of(plus)), ("-", hypotenuse_of(minus))))
    assert gamma.primary is None

    density = VolumeElement(gamma, on_the_facet).expand_geometry()
    assert density.as_float() == pytest.approx(math.sqrt(2.0))

    here = as_numbers(SpatialCoordinate(gamma, on_the_facet), 2)
    assert here == pytest.approx([0.3, 0.7])

    answers = {"+": 1.0 - 0.3, "-": 0.3}
    for side in gamma.embeddings:
        transition = side.transition(on_the_facet)
        assert isinstance(transition, PreimagePoint)
        through_side = expand_geometry(SpatialCoordinate(side.image, transition))
        answer = Point([RealScalar(answers[side.label])], entity_domain(Interval()))
        substituted = replace(through_side, {transition: answer})
        assert isinstance(substituted, AbstractExpression)
        assert [as_number(substituted.component(i)) for i in range(2)] == pytest.approx(here)

    root_half = 1.0 / math.sqrt(2.0)
    out_of_plus = OutwardNormal(gamma.embedding("+")).reconstruct_with_variable(on_the_facet)
    out_of_minus = OutwardNormal(gamma.embedding("-")).reconstruct_with_variable(on_the_facet)
    assert as_numbers(out_of_plus, 2) == pytest.approx([root_half, root_half])
    assert as_numbers(out_of_minus, 2) == pytest.approx([-root_half, -root_half])


def test_a_multiplier_on_the_interface_couples_the_two_sides(space, mesh):
    """A field living on the interface, with the two sides read on it through preimages.

    The interface is a mesh of intervals of its own, so a function space
    on it is an ordinary one, and a mortar coupling is a jump against a
    multiplier. Pulling back gives one kernel per pair of facets, the
    multiplier read at the interval's point and both sides at the point
    their own chart gives, which is the library's to compute.
    """
    along = parametrized_domain(LagrangeElement(Interval(), 1, (2,)))
    multipliers = function_space(along, LagrangeElement(Interval(), 1))
    facets = sub_entities_of(mesh, 1)
    gamma = interface(along, (("+", facets), ("-", facets)))
    assert gamma.primary is None

    u, mu = Coefficient(space), TestFunction(multipliers)
    coupling = mu * jump(u, gamma) * dx(gamma)
    assert isinstance(coupling, Integral)
    with pytest.raises(ValueError, match="through a map"):
        mu * u * dx(gamma)

    flux = mu * inner(jump(grad(u), gamma), OutwardNormal(gamma.embedding("+"))) * dx(gamma)
    pulled = pull_back_to_entity(flux)
    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == 9
    for pair, term in zip(gamma.pieces, pulled.terms, strict=True):
        assert isinstance(pair, Interface)
        assert pair.domain == along
        assert term.domain == EntityDomain(Interval())
        variable = term.variable.to_entity_coordinates(Interval())
        read_at = [
            node.variable.entity_point
            for node in as_graph(term.integrand)
            if isinstance(node, JacobianInverse) and isinstance(node.variable, ImagePoint)
        ]
        sides = set()
        for point in read_at:
            assert isinstance(point, PreimagePoint)
            assert point.ambient_point == ImagePoint(variable, along)
            assert isinstance(point.parametrization, LabelledParametrization)
            sides.add(point.parametrization.label.label)
        assert sides == {"+", "-"}


def test_a_reflected_neighbour_has_the_same_facet_map_and_the_opposite_outward_normal(
    plus, on_the_facet
):
    """INTEGRATION.md section 4's claim, as sharp as it gets.

    Reflecting the triangle across its hypotenuse gives a neighbour whose
    facet is parametrized by the very same points in the very same order,
    so the cross product of the facet's own Jacobian is identical on the
    two sides and cannot tell them apart. The reference outward normal
    carried by ``J^-T`` can: the reflection flips it, as outwardness
    requires.
    """
    reflected = composed_domain(
        entity_domain(PlacedTriangle()),
        AffineParametrization(RD(2), RD(2), ((0.0, -1.0), (-1.0, 0.0)), (1.0, 1.0)),
    )
    gamma = interface(
        hypotenuse_of(plus), (("+", hypotenuse_of(plus)), ("-", hypotenuse_of(reflected)))
    )

    assert as_numbers(
        SpatialCoordinate(hypotenuse_of(reflected), on_the_facet), 2
    ) == pytest.approx(as_numbers(SpatialCoordinate(hypotenuse_of(plus), on_the_facet), 2))
    assert as_numbers(UnitNormal(hypotenuse_of(reflected), on_the_facet), 2) == pytest.approx(
        as_numbers(UnitNormal(hypotenuse_of(plus), on_the_facet), 2)
    )

    out_of_plus = OutwardNormal(gamma.embedding("+")).reconstruct_with_variable(on_the_facet)
    out_of_minus = OutwardNormal(gamma.embedding("-")).reconstruct_with_variable(on_the_facet)
    assert as_numbers(out_of_minus, 2) == pytest.approx([-c for c in as_numbers(out_of_plus, 2)])


def test_two_intervals_meet_at_a_point():
    """An interface in one dimension: a point, counted once, with a normal that is a sign.

    The normal to a point is a convention, which is why ``UnitNormal``
    refuses; the entity supplies the convention as which way is out of
    the interval, and the outward normal carries it.
    """
    left = composed_domain(
        entity_domain(PlacedInterval()),
        AffineParametrization(RD(1), RD(1), ((1.0,),), (0.0,)),
    )
    right = composed_domain(
        entity_domain(PlacedInterval()),
        AffineParametrization(RD(1), RD(1), ((1.0,),), (1.0,)),
    )
    end_of_left = composed_domain(
        sub_entity_domain(PlacedInterval(), 0, 1), left.parametrization(PlacedInterval())
    )
    start_of_right = composed_domain(
        sub_entity_domain(PlacedInterval(), 0, 0), right.parametrization(PlacedInterval())
    )
    gamma = interface(end_of_left, (("left", end_of_left), ("right", start_of_right)))
    at_the_point = Point([], entity_domain(PointEntity()))

    assert gamma.topological_dimension == 0
    assert VolumeElement(gamma, at_the_point).expand_geometry().as_float() == pytest.approx(1.0)
    with pytest.raises(NotImplementedError, match="is a sign"):
        UnitNormal(gamma, at_the_point).value_shape

    out_of_left = OutwardNormal(gamma.embedding("left")).reconstruct_with_variable(at_the_point)
    out_of_right = OutwardNormal(gamma.embedding("right")).reconstruct_with_variable(at_the_point)
    assert as_numbers(out_of_left, 1) == pytest.approx([1.0])
    assert as_numbers(out_of_right, 1) == pytest.approx([-1.0])


def test_a_one_dimensional_mesh_has_an_interior_facet_integral():
    """``dS`` on a mesh of intervals: four pairs of end points, each a kernel over a point."""
    mesh = parametrized_domain(LagrangeElement(PlacedInterval(), 1, (1,)))
    space = function_space(mesh, LagrangeElement(PlacedInterval(), 1))
    u, v = Coefficient(space), TestFunction(space)
    gamma = interior_facets_of(mesh)

    pulled = pull_back_to_entity(jump(u, gamma) * restrict(v, gamma.embedding("+")) * dx(gamma))
    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == 4
    assert all(term.domain == EntityDomain(PointEntity()) for term in pulled.terms)


def test_a_triangle_meets_a_quadrilateral():
    """Two meshes of different cell types, sharing interval facets: three by four pairs."""
    triangles = parametrized_domain(LagrangeElement(PlacedTriangle(), 1, (2,)))
    quadrilaterals = parametrized_domain(LagrangeElement(PlacedQuadrilateral(), 1, (2,)))
    gamma = interface(
        sub_entities_of(triangles, 1),
        (("+", sub_entities_of(triangles, 1)), ("-", sub_entities_of(quadrilaterals, 1))),
    )
    assert len(gamma.pieces) == 12

    u = Coefficient(function_space(triangles, LagrangeElement(PlacedTriangle(), 1)))
    w = Coefficient(function_space(quadrilaterals, LagrangeElement(PlacedQuadrilateral(), 1)))
    v = TestFunction(function_space(triangles, LagrangeElement(PlacedTriangle(), 1)))
    across = restrict(u, gamma.embedding("+")) - restrict(w, gamma.embedding("-"))
    with pytest.raises(ValueError, match="is not a sub-entity of a cell of"):
        restrict(w, gamma.embedding("+"))

    pulled = pull_back_to_entity(across * restrict(v, gamma.embedding("+")) * dx(gamma))
    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == 12
    pairs = []
    for term in pulled.terms:
        pair = quantities_in(term.integrand)["VolumeElement"].domain
        assert isinstance(pair, Interface)
        pairs.append(
            (
                included_sub_entity(pair.image("+")),
                included_sub_entity(pair.image("-")),
            )
        )
    assert pairs == [((1, k), (1, m)) for k in range(3) for m in range(4)]


def dot(a, b) -> float:
    """The dot product of two lists of numbers."""
    return sum(x * y for x, y in zip(a, b, strict=True))


class Sheared(AbstractParametrization):
    """The map of the plane ``(x, y) -> (x, y + x^2)``, which bends every line not vertical."""

    @property
    def source(self) -> AbstractCoordinateDomain:
        """The plane."""
        return RD(2)

    @property
    def target(self) -> AbstractCoordinateDomain:
        """The plane."""
        return RD(2)

    def value(self, point: AbstractVariable) -> AbstractExpression:
        """Add the square of the first coordinate to the second."""
        x, y = point.component(0), point.component(1)
        return Vector([x, y + x * x])

    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """A lower triangular matrix that depends on the point."""
        x = point.component(0)
        return Matrix([[RealScalar(1.0), RealScalar(0.0)], [RealScalar(2.0) * x, RealScalar(1.0)]])

    def __repr__(self) -> str:
        """Representation."""
        return "Sheared()"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, Sheared)

    def __hash__(self) -> int:
        """Hash."""
        return hash("Sheared")


def test_a_curved_interface_has_arc_length_and_opposite_outward_normals(plus, minus, on_the_facet):
    """Two triangles bent by one map of the plane, sharing an edge that is now a curve.

    The density is the curve's arc length and not the flat facet's, and
    either side induces it. The outward normals are unit, orthogonal to
    the curve, opposite, and really outward: the ``+`` side's points
    away from a point inside the ``+`` triangle. The ``-`` side's
    Jacobian depends on its point, so its preimage has to be answered
    before anything is a number; the answer is ``1 - s``.
    """
    bent_plus = composed_domain(plus, Sheared())
    bent_minus = composed_domain(minus, Sheared())
    gamma = interface(
        hypotenuse_of(bent_plus),
        (("+", hypotenuse_of(bent_plus)), ("-", hypotenuse_of(bent_minus))),
    )

    s = 0.3
    # The curve is (1 - s, s + (1 - s)^2), whose speed is this.
    arc_speed = math.sqrt(1.0 + (2.0 * s - 1.0) ** 2)
    assert arc_speed != pytest.approx(math.sqrt(2.0))
    assert VolumeElement(gamma, on_the_facet).expand_geometry().as_float() == pytest.approx(
        arc_speed
    )
    answer = Point([RealScalar(1.0 - s)], entity_domain(Interval()))
    other_way = VolumeElement(hypotenuse_of(bent_minus), answer).expand_geometry()
    assert other_way.as_float() == pytest.approx(arc_speed)

    out_of_plus = OutwardNormal(gamma.embedding("+")).reconstruct_with_variable(on_the_facet)
    out_of_minus = OutwardNormal(gamma.embedding("-")).reconstruct_with_variable(on_the_facet)
    nu_plus = as_numbers(out_of_plus, 2)
    expanded = expand_geometry(out_of_minus)
    substituted = replace(expanded, {preimage_in(out_of_minus.variable): answer})
    assert isinstance(substituted, AbstractExpression)
    nu_minus = [as_number(substituted.component(i)) for i in range(2)]

    tangent = [-1.0, 2.0 * s - 1.0]
    assert dot(nu_plus, nu_plus) == pytest.approx(1.0)
    assert dot(nu_plus, tangent) == pytest.approx(0.0)
    assert nu_minus == pytest.approx([-c for c in nu_plus])

    here = [1.0 - s, s + (1.0 - s) ** 2]
    inside_plus = [0.2, 0.2 + 0.2**2]
    assert dot(nu_plus, [a - b for a, b in zip(inside_plus, here, strict=True)]) < 0


@pytest.fixture
def flat_in_space():
    """The reference triangle laid in the plane z = 0 of space."""
    return composed_domain(
        entity_domain(PlacedTriangle()),
        AffineParametrization(RD(2), RD(3), ((1.0, 0.0), (0.0, 1.0), (0.0, 0.0)), (0.0, 0.0, 0.0)),
    )


def folded_onto(apex):
    """A triangle in space sharing the flat one's hypotenuse as facet 2, with vertex 0 at the apex.

    An affine map sending vertex 1 to (1, 0, 0), vertex 2 to (0, 1, 0) and
    vertex 0 to the apex, so the triangle shares the edge and leans out
    of the plane.
    """
    matrix = [[(1.0, 0.0, 0.0)[i] - apex[i], (0.0, 1.0, 0.0)[i] - apex[i]] for i in range(3)]
    return composed_domain(
        entity_domain(PlacedTriangle()), AffineParametrization(RD(2), RD(3), matrix, apex)
    )


def test_two_triangles_meeting_at_an_angle_have_conormals_and_not_a_normal(
    flat_in_space, on_the_facet
):
    """A manifold interface: the outward normals are conormals, each in its own plane.

    The second triangle leans out of the plane, so the two tangent planes
    differ along the shared edge and the two conormals are not opposite.
    Each is unit, orthogonal to the edge and in its own triangle's plane,
    which is what out of a surface across its edge means; an ambient
    normal to the edge is not a single direction at all.
    """
    leaning = folded_onto((1.0, 1.0, 1.0))
    gamma = interface(
        hypotenuse_of(flat_in_space),
        (("flat", hypotenuse_of(flat_in_space)), ("leaning", hypotenuse_of(leaning))),
    )
    assert as_numbers(SpatialCoordinate(hypotenuse_of(leaning), on_the_facet), 3) == pytest.approx(
        as_numbers(SpatialCoordinate(gamma, on_the_facet), 3)
    )
    with pytest.raises(ValueError, match="codimension one"):
        UnitNormal(gamma, on_the_facet).value_shape

    nu_flat = as_numbers(
        OutwardNormal(gamma.embedding("flat")).reconstruct_with_variable(on_the_facet), 3
    )
    nu_leaning = as_numbers(
        OutwardNormal(gamma.embedding("leaning")).reconstruct_with_variable(on_the_facet), 3
    )
    tangent = [-1.0, 1.0, 0.0]

    assert nu_flat == pytest.approx([1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0), 0.0])
    assert dot(nu_leaning, nu_leaning) == pytest.approx(1.0)
    assert dot(nu_leaning, tangent) == pytest.approx(0.0)
    # Away from the apex, within the leaning plane: minus (1/2, 1/2, 1) normalised.
    assert nu_leaning == pytest.approx([-c / math.sqrt(1.5) for c in (0.5, 0.5, 1.0)])
    assert dot(nu_flat, nu_leaning) == pytest.approx(-1.0 / math.sqrt(3.0))


def test_three_surfaces_meeting_along_an_edge(flat_in_space, on_the_facet):
    """A T-junction: an interface of three sides, which a jump refuses and a sum writes.

    Three triangles in space share one edge. Each has its outward conormal
    along it, a field on each is read on the edge through its own side,
    and the integral pulls back to one kernel reading three cells, one
    directly and two through preimages. A jump is of two sides, so it
    says so; a Kirchhoff-type sum of the three is written as the three
    restrictions it is.
    """
    up = folded_onto((1.0, 1.0, 1.0))
    down = folded_onto((1.0, 1.0, -1.0))
    gamma = interface(
        hypotenuse_of(flat_in_space),
        (
            ("flat", hypotenuse_of(flat_in_space)),
            ("up", hypotenuse_of(up)),
            ("down", hypotenuse_of(down)),
        ),
    )
    assert gamma.labels == ("flat", "up", "down")
    assert gamma.pieces == (gamma,)

    normals = {
        label: as_numbers(
            OutwardNormal(gamma.embedding(label)).reconstruct_with_variable(on_the_facet), 3
        )
        for label in gamma.labels
    }
    for normal in normals.values():
        assert dot(normal, normal) == pytest.approx(1.0)
        assert dot(normal, [-1.0, 1.0, 0.0]) == pytest.approx(0.0)
    assert normals["up"][2] == pytest.approx(-normals["down"][2])

    fields = {
        label: Coefficient(function_space(domain, LagrangeElement(PlacedTriangle(), 1)))
        for label, domain in (("flat", flat_in_space), ("up", up), ("down", down))
    }
    with pytest.raises(ValueError, match="has 3 sides"):
        jump(fields["flat"], gamma)
    with pytest.raises(ValueError, match="has 3 sides"):
        avg(fields["flat"], gamma)

    total = expression_sum(
        restrict(fields[label], gamma.embedding(label)) for label in gamma.labels
    )
    assert extract_domain(total) == gamma
    pulled = pull_back_to_entity(total * dx(gamma))
    assert isinstance(pulled, Integral)
    assert pulled.domain == EntityDomain(Interval())
    read_through = set()
    for node in as_graph(pulled.integrand):
        if isinstance(node, Coefficient):
            assert isinstance(node.variable, ImagePoint)
            point = node.variable.entity_point
            if isinstance(point, PreimagePoint):
                assert isinstance(point.parametrization, LabelledParametrization)
                read_through.add(point.parametrization.label.label)
            else:
                read_through.add("flat")
    assert read_through == {"flat", "up", "down"}


def read_through(term: Integral) -> dict[str, str | None]:
    """Which side each function in a pulled back term is read through, or None if directly."""
    sides: dict[str, str | None] = {}
    for node in as_graph(term.integrand):
        if isinstance(node, Coefficient | TestFunction):
            assert isinstance(node.variable, ImagePoint)
            point = node.variable.entity_point
            if isinstance(point, PreimagePoint):
                assert isinstance(point.parametrization, LabelledParametrization)
                sides[type(node).__name__] = point.parametrization.label.label
            else:
                assert isinstance(point.domain, EntityDomain)
                (entity,) = point.domain.cell_types
                assert point == term.variable.to_entity_coordinates(entity)
                sides[type(node).__name__] = None
    return sides


def test_a_mesh_and_a_submesh_of_its_cells_meet_as_an_interface_of_codimension_zero():
    """A mixed domain form: a field on a mesh against a field on a submesh of its cells.

    Two meshes whose cells coincide are two images of one set, so
    they are an interface like any other, of codimension zero. A cell is
    its own sub-entity there, so the pieces are the cells mapped on, and
    the submesh's field is read at the point its own cell gives, which is
    the entity map between the two meshes, a preimage the library answers.
    """
    tri = PlacedTriangle()
    parent = parametrized_domain(LagrangeElement(tri, 1, (2,)))
    submesh = parametrized_domain(LagrangeElement(tri, 2, (2,)))
    u = Coefficient(function_space(parent, LagrangeElement(tri, 1)))
    v = TestFunction(function_space(submesh, LagrangeElement(tri, 1)))

    cells = sub_entities_of(parent, 0, "marked")
    assert cells == union_domain(
        [LabelledDomain(composed_domain(entity_domain(tri), parent.parametrization(tri)), "marked")]
    )
    gamma = interface(cells, (("mesh", cells), ("submesh", sub_entities_of(submesh, 0, "marked"))))
    assert gamma.primary == "mesh"
    assert gamma.topological_dimension == 2
    assert len(gamma.pieces) == 1

    with pytest.raises(ValueError, match="through a map"):
        u * restrict(v, gamma.embedding("submesh")) * dx(gamma)

    pulled = pull_back_to_entity(
        restrict(u, gamma.embedding("mesh")) * restrict(v, gamma.embedding("submesh")) * dx(gamma)
    )
    assert isinstance(pulled, IntegralSum)
    (term,) = pulled.terms
    assert isinstance(term, Integral)
    assert term.domain == EntityDomain(tri)
    assert read_through(term) == {"Coefficient": None, "TestFunction": "submesh"}


def test_a_mesh_and_a_submesh_of_its_facets_meet_along_those_facets():
    """A field on the boundary against a field in the bulk, the facets being the submesh's cells.

    The mesh presents the facets as its sub-entities of codimension one
    and the submesh as its cells, so there is one pair per facet of the
    mesh. The bulk field is read through the facet's inclusion into its
    cell, as under ``ds``, and the boundary field at the point its own
    cell gives, through the map between the facet and the submesh cell.
    """
    tri = PlacedTriangle()
    mesh = parametrized_domain(LagrangeElement(tri, 1, (2,)))
    boundary = parametrized_domain(LagrangeElement(Interval(), 1, (2,)))
    u = Coefficient(function_space(mesh, LagrangeElement(tri, 1)))
    mu = TestFunction(function_space(boundary, LagrangeElement(Interval(), 1)))

    facets = sub_entities_of(mesh, 1)
    gamma = interface(facets, (("mesh", facets), ("submesh", sub_entities_of(boundary, 0))))
    assert gamma.topological_dimension == 1
    assert len(gamma.pieces) == 3

    coupling = inner(
        restrict(grad(u), gamma.embedding("mesh")), OutwardNormal(gamma.embedding("mesh"))
    )
    pulled = pull_back_to_entity(coupling * restrict(mu, gamma.embedding("submesh")) * dx(gamma))
    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == 3
    for pair, term in zip(gamma.pieces, pulled.terms, strict=True):
        assert isinstance(pair, Interface)
        assert isinstance(term, Integral)
        assert term.domain == EntityDomain(Interval())
        assert included_sub_entity(pair.image("mesh")) == included_sub_entity(pair.domain)
        assert read_through(term) == {"Coefficient": None, "TestFunction": "submesh"}
        inverses = [node for node in as_graph(term.integrand) if isinstance(node, JacobianInverse)]
        assert len(inverses) == 1
        assert inverses[0].domain == mesh


def test_a_submesh_is_of_entities_of_its_mesh():
    """Cells of a higher dimension than the entities they are said to be are refused."""
    tri, tet = PlacedTriangle(), PlacedTetrahedron()
    tetrahedra = parametrized_domain(LagrangeElement(tet, 1, (3,)))
    faces = parametrized_domain(LagrangeElement(tri, 1, (3,)))
    edges = parametrized_domain(LagrangeElement(Interval(), 1, (3,)))

    with pytest.raises(ValueError, match="topological dimension"):
        interface(
            sub_entities_of(tetrahedra, 2),
            (("mesh", sub_entities_of(tetrahedra, 2)), ("submesh", sub_entities_of(faces, 0))),
        )

    ridges = sub_entities_of(tetrahedra, 2)
    along_the_edges = interface(ridges, (("mesh", ridges), ("submesh", sub_entities_of(edges, 0))))
    assert len(along_the_edges.pieces) == 6


def embedded(gamma, omega, by=None):
    """A domain presented from a mesh: carried into each cell type's coordinates and mapped on.

    Written out, as the language has no spelling for it: the one
    construction is ``interface(domain, sides)``, and this is what a side
    is when the interface is a facet of no cell of the mesh.
    """
    (region,) = {parametrization.source for parametrization in gamma.parametrizations}
    return union_domain(
        [
            composed_domain(
                composed_domain(
                    region, CoordinateRepresentation(gamma, omega.parametrization(cell), by)
                ),
                omega.parametrization(cell),
            )
            for cell in omega.cell_types
        ]
    )


def answered(expression, answer):
    """An expression with its preimage point answered, as the library would.

    An affine map's Jacobian reads no point, so the preimage may have
    dropped out already, and then there is nothing to answer.
    """
    expanded = expand_geometry(expression)
    replacements: dict[GraphNode, GraphNode] = {
        node: answer for node in as_graph(expanded) if isinstance(node, PreimagePoint)
    }
    assert len(replacements) <= 1
    substituted = replace(expanded, replacements) if replacements else expanded
    assert isinstance(substituted, AbstractExpression)
    return substituted


def test_a_periodic_boundary_is_an_interface_whose_other_side_is_a_translation_away(
    plus, on_the_facet
):
    """The left edge of one triangle meets the left edge of its translate, through the shift.

    The right side's map is ``phi_B^-1 . T . psi_Gamma``: the edge's point
    shifted by one and pulled back into the translated triangle, where it
    lands on that triangle's own left edge. The preimage is the library's
    to answer, and the answer here is ``(0, s)``. The Jacobian is the chain
    rule and comes out as the edge's own tangent, the shift being rigid.
    """
    shift = AffineParametrization(RD(2), RD(2), ((1.0, 0.0), (0.0, 1.0)), (1.0, 0.0))
    translate = composed_domain(entity_domain(PlacedTriangle()), shift)
    left_edge = composed_domain(
        sub_entity_domain(PlacedTriangle(), 1, 1), plus.parametrization(PlacedTriangle())
    )
    gamma = interface(
        left_edge, (("left", left_edge), ("right", embedded(left_edge, translate, shift)))
    )

    assert gamma.primary == "left"
    (pair,) = gamma.pieces
    assert isinstance(pair, Interface)
    right = pair.embedding("right")
    assert right.transition(on_the_facet) is on_the_facet
    (carried,) = right.image.parametrizations
    assert isinstance(carried, ComposedParametrization)
    assert isinstance(carried.inner, CoordinateRepresentation)
    assert carried.inner.carries == left_edge

    u = Coefficient(function_space(translate, LagrangeElement(PlacedTriangle(), 1)))
    image = Restriction(u, right, on_the_facet).image_point
    assert image.domain == entity_domain(PlacedTriangle())
    expanded = expand_geometry(image)
    preimage = preimage_in(expanded)
    assert preimage.parametrization == translate.parametrization(PlacedTriangle())
    shifted = expand_geometry(preimage.ambient_point)
    assert isinstance(shifted, AbstractExpression)
    assert [as_number(shifted.component(i)) for i in range(2)] == pytest.approx([1.0, 0.3])

    answer = Point([RealScalar(0.0), RealScalar(0.3)], entity_domain(PlacedTriangle()))
    assert as_numbers(SpatialCoordinate(translate, answer), 2) == pytest.approx([1.0, 0.3])

    tangent = answered(carried.inner.jacobian(on_the_facet), answer)
    assert [as_number(tangent.component(i, 0)) for i in range(2)] == pytest.approx([0.0, 1.0])
    assert answered(VolumeElement(right.image, on_the_facet), answer).as_float() == (
        pytest.approx(1.0)
    )

    with pytest.raises(ValueError, match="is not a facet"):
        OutwardNormal(right).reconstruct_with_variable(on_the_facet).expand_geometry()


def test_an_interface_that_is_a_facet_of_neither_side_is_read_through_preimages(
    plus, minus, along_the_edge, on_the_facet
):
    """Both sides carried, so no side is primary and both are read at the interface's point.

    The answers are ``1 - s`` for the ``+`` side and ``s`` for the ``-``,
    as in the matching case, but now they sit inside the sides' maps
    rather than in a chart transition.
    """
    gamma = interface(
        along_the_edge,
        (("+", embedded(along_the_edge, plus)), ("-", embedded(along_the_edge, minus))),
    )
    assert gamma.primary is None
    (pair,) = gamma.pieces
    assert isinstance(pair, Interface)
    here = as_numbers(SpatialCoordinate(pair, on_the_facet), 2)

    # The preimage lands in each triangle's own coordinates, where the
    # point (0.3, 0.7) of the plane sits at itself in the first and at its
    # reflection in the second.
    answers = {"+": (0.3, 0.7), "-": (0.7, 0.3)}
    for side in pair.embeddings:
        assert side.transition(on_the_facet) is on_the_facet
        answer = Point(
            [RealScalar(c) for c in answers[side.label]], entity_domain(PlacedTriangle())
        )
        through_side = answered(SpatialCoordinate(side.image, on_the_facet), answer)
        assert [as_number(through_side.component(i)) for i in range(2)] == pytest.approx(here)


def test_a_mortar_coupling_of_two_meshes_through_a_multiplier(on_the_facet):
    """Two meshes of different geometry, an interface mesh of their own, a multiplier on it.

    One kernel over the reference interval: the multiplier read at the
    interval's point, each mesh's field at the preimage through its own
    map, which is point location the library does.
    """
    tri = PlacedTriangle()
    coarse = parametrized_domain(LagrangeElement(tri, 1, (2,)))
    fine = parametrized_domain(LagrangeElement(tri, 2, (2,)))
    along = parametrized_domain(LagrangeElement(Interval(), 1, (2,)))
    u = Coefficient(function_space(coarse, LagrangeElement(tri, 1)))
    w = Coefficient(function_space(fine, LagrangeElement(tri, 1)))
    mu = TestFunction(function_space(along, LagrangeElement(Interval(), 1)))

    gamma = interface(along, (("coarse", embedded(along, coarse)), ("fine", embedded(along, fine))))
    coupling = (
        mu
        * (restrict(u, gamma.embedding("coarse")) - restrict(w, gamma.embedding("fine")))
        * dx(gamma)
    )

    pulled = pull_back_to_entity(coupling)
    assert isinstance(pulled, IntegralSum)
    (term,) = pulled.terms
    assert isinstance(term, Integral)
    assert term.domain == EntityDomain(Interval())
    # Each mesh's field is read at the interval's point carried into that
    # mesh's cell, and the carrying map is the preimage through that mesh.
    through = []
    for node in as_graph(term.integrand):
        if isinstance(node, Coefficient):
            assert isinstance(node.variable, ImagePoint)
            assert node.variable.entity_point == term.variable.to_entity_coordinates(Interval())
            carried = node.variable.parametrization
            assert isinstance(carried, CoordinateRepresentation)
            assert carried.carries == along
            through.append(carried.through)
        elif isinstance(node, TestFunction):
            assert node.variable == term.variable.to_entity_coordinates(Interval())
    assert sorted(map(repr, through)) == sorted(
        map(repr, [coarse.parametrization(tri), fine.parametrization(tri)])
    )


def test_a_well_in_three_dimensions_is_an_interface_of_codimension_two():
    """A line of intervals in a mesh of tetrahedra, which is a facet of nothing.

    The rock's field is read at the point of the tetrahedron the well
    passes through, which is where the kernel lands; there is no normal
    to a line in space and no out of a cell across it.
    """
    well = parametrized_domain(LagrangeElement(Interval(), 1, (3,)))
    rock = parametrized_domain(LagrangeElement(PlacedTetrahedron(), 1, (3,)))
    p = Coefficient(function_space(rock, LagrangeElement(PlacedTetrahedron(), 1)))
    lam = TestFunction(function_space(well, LagrangeElement(Interval(), 1)))
    s = Point([RealScalar(0.3)], entity_domain(Interval()))

    gamma = interface(well, (("well", well), ("rock", embedded(well, rock))))
    assert gamma.topological_dimension == 1
    assert gamma.embedding("rock").geometric_dimension == 3

    pulled = pull_back_to_entity(restrict(p, gamma.embedding("rock")) * lam * dx(gamma))
    assert isinstance(pulled, IntegralSum)
    (term,) = pulled.terms
    assert isinstance(term, Integral)
    (field,) = [node for node in as_graph(term.integrand) if isinstance(node, Coefficient)]
    assert isinstance(field.variable, ImagePoint)
    assert field.variable.domain == entity_domain(PlacedTetrahedron())

    with pytest.raises(ValueError, match="codimension one"):
        UnitNormal(gamma, s).value_shape
    (pair,) = gamma.pieces
    assert isinstance(pair, Interface)
    with pytest.raises(ValueError, match="is not a facet"):
        OutwardNormal(pair.embedding("rock")).reconstruct_with_variable(s).expand_geometry()


def test_a_union_presented_piece_by_piece_pairs_each_piece_with_its_own():
    """The facets of one mesh carried into another: three pairs and not nine.

    A carried piece says which piece it carries, so the fan out pairs
    it with that one alone.
    """
    tri = PlacedTriangle()
    one = parametrized_domain(LagrangeElement(tri, 1, (2,)))
    other = parametrized_domain(LagrangeElement(tri, 2, (2,)))
    facets = sub_entities_of(one, 1)
    carried = union_domain(
        [piece for facet in facets.pieces for piece in embedded(facet, other).pieces]
    )

    gamma = interface(facets, (("one", facets), ("other", carried)))
    pairs = gamma.pieces
    assert len(pairs) == 3
    for facet, pair in zip(facets.pieces, pairs, strict=True):
        assert isinstance(pair, Interface)
        assert pair.domain == facet
        (map_of_other,) = pair.image("other").parametrizations
        assert isinstance(map_of_other, ComposedParametrization)
        assert isinstance(map_of_other.inner, CoordinateRepresentation)
        assert map_of_other.inner.carries == facet


def test_a_domain_is_carried_from_one_cells_coordinates_by_a_map_that_joins_them(plus):
    """What a coordinate representation refuses: a union, no map between; and a box is fine."""
    tri = PlacedTriangle()
    box_face = Box(((0.0, 1.0), (0.0, 1.0))).face(axis=0, upper=True)
    mesh = parametrized_domain(LagrangeElement(tri, 1, (2,)))
    in_space = parametrized_domain(LagrangeElement(tri, 1, (3,)))
    edge = composed_domain(sub_entity_domain(tri, 1, 1), plus.parametrization(tri))

    # The map needs no cell on either end; presenting its image as a domain does.
    of_a_box = CoordinateRepresentation(box_face, mesh.parametrization(tri))
    assert of_a_box.source == box_face.parameter_region
    assert of_a_box.target == entity_domain(tri)
    with pytest.raises(ValueError, match="piece by piece"):
        CoordinateRepresentation(sub_entities_of(mesh, 1), mesh.parametrization(tri))
    with pytest.raises(ValueError, match="without a map between them"):
        CoordinateRepresentation(edge, in_space.parametrization(tri))
    with pytest.raises(ValueError, match="does not take the coordinates"):
        CoordinateRepresentation(
            edge,
            in_space.parametrization(tri),
            AffineParametrization(RD(3), RD(3), ((1, 0, 0), (0, 1, 0), (0, 0, 1)), (0, 0, 0)),
        )
    lift = AffineParametrization(RD(2), RD(3), ((1, 0), (0, 1), (0, 0)), (0, 0, 0))
    assert CoordinateRepresentation(
        edge, in_space.parametrization(tri), lift
    ).target == entity_domain(tri)


def test_an_interface_immersed_in_one_mesh_has_no_jump_to_see():
    """Both sides carried into the same mesh by the same map are the same side.

    Point location names one cell for a point of the interface, so the
    two restrictions of a field on that mesh are read at one point and
    the pulled back jump is a difference of one node with itself. That
    is right: without an enriched space there is no jump across a line
    through a cell, and the language says so by making the two sides
    one rather than by pretending otherwise.
    """
    tri = PlacedTriangle()
    mesh = parametrized_domain(LagrangeElement(tri, 1, (2,)))
    along = parametrized_domain(LagrangeElement(Interval(), 1, (2,)))
    u = Coefficient(function_space(mesh, LagrangeElement(tri, 1)))
    mu = TestFunction(function_space(along, LagrangeElement(Interval(), 1)))

    gamma = interface(along, (("+", embedded(along, mesh)), ("-", embedded(along, mesh))))
    assert gamma.image("+") == gamma.image("-")
    assert gamma.embedding("+") != gamma.embedding("-")

    pulled = pull_back_to_entity(mu * jump(u, gamma) * dx(gamma))
    assert isinstance(pulled, IntegralSum)
    (term,) = pulled.terms
    assert isinstance(term, Integral)
    fields = [node for node in as_graph(term.integrand) if isinstance(node, Coefficient)]
    assert len(fields) == 1
