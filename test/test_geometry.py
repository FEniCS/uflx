"""Test geometry."""

import pytest

from uflx import Coefficient, TestFunction, dx, function_space, inner, parametrized_domain
from uflx.algorithms import pull_back_to_entity, simplify
from uflx.basis_functions import EvaluatedBasisFunction
from uflx.domains import RD, EntityDomain
from uflx.expressions import RealScalar, Sqrt
from uflx.functions import create_variable
from uflx.geometry import (
    AbstractGeometricQuantity,
    Jacobian,
    JacobianDeterminant,
    JacobianInverse,
    JacobianInverseTranspose,
    JacobianTranspose,
    MetricTensor,
    PulledBackPoint,
    PushedForwardPoint,
    SingleSpatialCoordinate,
    SpatialCoordinate,
    TangentialProjector,
    UnitNormal,
    VolumeElement,
    expand_geometry,
)
from uflx.graphs import as_graph
from uflx.integrals import Integral
from uflx.points import Point
from uflx.tensors import FlattenedTensorMap, Identity, Matrix

cells_and_gdims = [
    ("interval", 1),
    ("interval", 2),
    ("interval", 3),
    ("interval", 4),
    ("triangle", 2),
    ("triangle", 3),
    ("triangle", 4),
    ("quadrilateral", 2),
    ("quadrilateral", 3),
    ("quadrilateral", 4),
    ("tetrahedron", 3),
    ("tetrahedron", 4),
    ("hexahedron", 3),
    ("hexahedron", 4),
]


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_expand_geometry(cell, gdim, lagrange_element):
    """Test expansion of Jacobian."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cell_types[0].topological_dimension

    point = Point([RealScalar(1.0)] * tdim, EntityDomain(domain.cell_types[0]))
    j = Jacobian(domain, point)
    assert j.value_shape == (gdim, tdim)

    mat = j.expand_geometry()
    assert mat.value_shape == (gdim, tdim)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_inverse_expand_geometry(cell, gdim, lagrange_element):
    """Test expansion of Jacobian inverse."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cell_types[0].topological_dimension

    point = Point([RealScalar(1.0)] * tdim, EntityDomain(domain.cell_types[0]))
    j = JacobianInverse(domain, point)
    assert j.value_shape == (tdim, gdim)

    mat = j.expand_geometry()
    assert mat.value_shape == (tdim, gdim)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_tranpose_expand_geometry(cell, gdim, lagrange_element):
    """Test expansion of Jacobian inverse transpose."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cell_types[0].topological_dimension

    point = Point([RealScalar(1.0)] * tdim, EntityDomain(domain.cell_types[0]))
    j = JacobianTranspose(domain, point)
    assert j.value_shape == (tdim, gdim)

    mat = j.expand_geometry()
    assert mat.value_shape == (tdim, gdim)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_inverse_transpose_expand_geometry(cell, gdim, lagrange_element):
    """Test expansion of Jacobian inverse transpose."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cell_types[0].topological_dimension

    point = Point([RealScalar(1.0)] * tdim, EntityDomain(domain.cell_types[0]))
    j = JacobianInverseTranspose(domain, point)
    assert j.value_shape == (gdim, tdim)

    mat = j.expand_geometry()
    assert mat.value_shape == (gdim, tdim)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_pushed_forward_point_expands_to_ambient_coordinates(cell, gdim, lagrange_element):
    """Expanding a pushed forward point gives explicit ambient coordinates."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    (entity,) = domain.cell_types
    tdim = entity.topological_dimension

    entity_point = Point([RealScalar(0.25)] * tdim, EntityDomain(entity))
    pushed = PushedForwardPoint(entity_point, domain)
    assert pushed.domain == RD(gdim)

    expanded = pushed.expand_geometry()
    assert isinstance(expanded, Point)
    assert expanded.domain == RD(gdim)
    assert not expanded.in_entity_coordinates


def test_pulled_back_point_lands_in_entity_coordinates(lagrange_element):
    """A pulled back point lies in the entity's coordinates, not on the mesh."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (3,)))
    (entity,) = domain.cell_types

    parametrization = domain.parametrization(entity)
    pulled = PulledBackPoint(Point([RealScalar(1.0)] * 3, RD(3)), parametrization)

    assert pulled.domain == EntityDomain(entity)
    assert pulled.in_entity_coordinates
    assert pulled.parametrization == parametrization


def test_the_two_mapped_points_do_not_collide(lagrange_element):
    """A pushed forward point is never equal to a pulled back one."""
    domain = parametrized_domain(lagrange_element("interval", 1, (1,)))
    (entity,) = domain.cell_types

    entity_point = Point([RealScalar(0.25)], EntityDomain(entity))
    ambient_point = Point([RealScalar(0.25)], RD(1))

    pushed = PushedForwardPoint(entity_point, domain)
    pulled = PulledBackPoint(ambient_point, domain.parametrization(entity))

    assert pushed != pulled
    assert pulled != pushed
    assert len({pushed, pulled}) == 2


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_derivatives_are_entity_derivatives(cell, gdim, lagrange_element):
    """Test that the Jacobian differentiates in the tdim directions of the cell's coordinates."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cell_types[0].topological_dimension

    point = Point([RealScalar(1.0)] * tdim, EntityDomain(domain.cell_types[0]))
    mat = Jacobian(domain, point).expand_geometry()
    derivatives = {
        node.derivative
        for node in as_graph(mat).ordered_nodes()
        if isinstance(node, EvaluatedBasisFunction)
    }
    assert derivatives == {tuple(int(d == i) for d in range(tdim)) for i in range(tdim)}


def coordinate_dof_entries(expression):
    """The (node, component) of each coordinate DOF an expression refers to."""
    return {
        node.init_args[0]
        for node in as_graph(expression).ordered_nodes()
        if isinstance(node, FlattenedTensorMap)
    }


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_coordinate_dofs(cell, gdim, lagrange_element):
    """Test that the Jacobian uses the coordinate DOFs of the map it differentiates."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cell_types[0].topological_dimension
    point = Point([RealScalar(0.1)] * tdim, EntityDomain(domain.cell_types[0]))

    x_dofs = coordinate_dof_entries(PushedForwardPoint(point, domain).expand_geometry())
    assert len(x_dofs) > 0
    assert coordinate_dof_entries(Jacobian(domain, point).expand_geometry()) == x_dofs


def mass_form(cell, gdim, lagrange_element):
    """A form with geometry in it once pulled back, but no gradients."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    space = function_space(domain, lagrange_element(cell, 1))
    return inner(Coefficient(space), TestFunction(space)) * dx


def test_a_jacobian_is_told_where_it_is_evaluated(lagrange_element):
    """A pulled back integral hands its geometry the variable standing for the point."""
    pulled = pull_back_to_entity(mass_form("triangle", 2, lagrange_element))

    assert isinstance(pulled, Integral)
    jacobians = [n for n in as_graph(pulled) if isinstance(n, AbstractGeometricQuantity)]
    assert len(jacobians) > 0
    for j in jacobians:
        assert j.point is not None
        assert j.point == pulled.variable
        assert j.point.domain == j.parametrization.source


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_a_pulled_back_integral_expands_its_geometry(cell, gdim, lagrange_element):
    """Expanding a pulled back form leaves no geometry behind.

    A Jacobian built during a pull back used to keep point=None, so
    expanding one asserted instead of giving an expression.
    """
    expanded = expand_geometry(pull_back_to_entity(mass_form(cell, gdim, lagrange_element)))

    assert not any(isinstance(n, AbstractGeometricQuantity) for n in as_graph(expanded))


def test_geometry_refuses_a_variable_in_the_wrong_coordinates(lagrange_element):
    """A Jacobian is evaluated at a point of its map's source, not anywhere else."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    (cell,) = domain.cell_types
    (other_cell,) = parametrized_domain(lagrange_element("interval", 1, (1,))).cell_types

    ambient = create_variable(domain)
    assert not ambient.in_entity_coordinates
    assert Jacobian(domain).reconstruct_with_variable(ambient).point is None

    foreign = (
        ambient.to_entity_coordinates(cell)
        .to_ambient_coordinates(parametrized_domain(lagrange_element("interval", 1, (1,))))
        .to_entity_coordinates(other_cell)
    )
    assert Jacobian(domain).reconstruct_with_variable(foreign).point is None

    entity = ambient.to_entity_coordinates(cell)
    assert Jacobian(domain).reconstruct_with_variable(entity).point == entity


@pytest.fixture
def mixed_mesh(lagrange_element):
    """A mesh of triangles and quadrilaterals, so no one cell type is the cell type."""
    return parametrized_domain(
        [lagrange_element("triangle", 1, (2,)), lagrange_element("quadrilateral", 1, (2,))]
    )


def test_a_jacobian_on_a_mixed_mesh_has_a_shape(mixed_mesh):
    """A domain with several cell types still has one dimension pair."""
    assert len(mixed_mesh.cell_types) == 2
    assert Jacobian(mixed_mesh).value_shape == (2, 2)


def test_a_jacobian_resolves_its_cell_from_its_point(mixed_mesh):
    """A point lies in a cell's coordinate domain, so it names the cell."""
    for cell in mixed_mesh.cell_types:
        tdim = cell.topological_dimension
        point = Point([RealScalar(0.25)] * tdim, EntityDomain(cell))

        j = Jacobian(mixed_mesh, point)

        assert j.parametrization == mixed_mesh.parametrization(cell)
        assert j.expand_geometry().value_shape == (2, 2)


def test_a_mixed_mesh_accepts_a_variable_on_any_of_its_cells(mixed_mesh):
    """Which cell a quantity is on is settled by the variable it is given."""
    ambient = create_variable(mixed_mesh)
    for cell in mixed_mesh.cell_types:
        variable = ambient.to_entity_coordinates(cell)

        told = Jacobian(mixed_mesh).reconstruct_with_variable(variable)

        assert told.point == variable
        assert told.parametrization == mixed_mesh.parametrization(cell)


def test_a_jacobian_without_a_point_knows_no_cell(mixed_mesh):
    """Being generic over cell types is useful, but it cannot be expanded."""
    with pytest.raises(ValueError, match="not been told where"):
        Jacobian(mixed_mesh).expand_geometry()


geometric_quantities = [
    Jacobian,
    JacobianDeterminant,
    VolumeElement,
    MetricTensor,
    TangentialProjector,
    UnitNormal,
    JacobianInverse,
    JacobianTranspose,
    JacobianInverseTranspose,
]


@pytest.mark.parametrize("quantity", geometric_quantities)
def test_a_geometric_quantity_shows_its_arguments(quantity, lagrange_element):
    """All of them say what they are of, not just their class name."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))

    assert repr(quantity(domain)).startswith(f"{quantity.__name__}(")
    assert repr(domain) in repr(quantity(domain))


@pytest.mark.parametrize("quantity", geometric_quantities)
def test_a_geometric_quantity_cannot_be_mutated(quantity, lagrange_element):
    """Its hash comes from its arguments, so moving them would lose it in a dict.

    These nodes are dictionary keys while an expression is being rewritten.
    """
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    (cell,) = domain.cell_types
    node = quantity(domain)
    found_by = {node: "here"}

    with pytest.raises(AttributeError):
        node.point = Point([RealScalar(0.25)] * 2, EntityDomain(cell))

    assert node in found_by


def test_pushing_forward_needs_a_point_in_a_cells_coordinates(lagrange_element):
    """There is nothing to carry a point of the ambient coordinates forward."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    ambient = Point([RealScalar(0.25)] * 2, RD(2))

    with pytest.raises(ValueError, match="pushed forward from a cell's coordinates"):
        PushedForwardPoint(ambient, domain)


def test_pulling_back_needs_a_point_in_ambient_coordinates(lagrange_element):
    """A point already in a cell's coordinates has nothing to pull back."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    (cell,) = domain.cell_types
    entity_point = Point([RealScalar(0.25)] * 2, EntityDomain(cell))

    with pytest.raises(ValueError, match="nothing to pull back"):
        PulledBackPoint(entity_point, domain.parametrization(cell))


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_the_metric_is_square_in_the_cells_own_dimension(cell, gdim, lagrange_element):
    """The metric measures in the cell's coordinates, so it is tdim by tdim."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    (entity,) = domain.cell_types
    tdim = entity.topological_dimension

    assert MetricTensor(domain).value_shape == (tdim, tdim)
    assert Jacobian(domain).value_shape == (gdim, tdim)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_the_metric_is_symmetric(cell, gdim, lagrange_element):
    """G = J^T J, so g[i, j] and g[j, i] are the same sum."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    (entity,) = domain.cell_types
    tdim = entity.topological_dimension
    point = Point([RealScalar(0.25)] * tdim, EntityDomain(entity))

    g = MetricTensor(domain, point).expand_geometry()

    assert isinstance(g, Matrix)
    for i in range(tdim):
        for j in range(tdim):
            assert simplify(g.component(i, j)) == simplify(g.component(j, i))


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_the_volume_element_is_the_metrics_gram_determinant(cell, gdim, lagrange_element):
    """On a manifold the factor an integral picks up is sqrt(det g).

    For a square map the same identity holds, but the volume element
    takes the direct and much cheaper abs(det J) instead, so the two are
    equal as numbers without being equal as expressions.
    """
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    (entity,) = domain.cell_types
    tdim = entity.topological_dimension
    if gdim == tdim:
        pytest.skip("Not a manifold, so the direct determinant is used.")
    point = Point([RealScalar(0.25)] * tdim, EntityDomain(entity))

    g = MetricTensor(domain, point).expand_geometry()
    assert isinstance(g, Matrix)

    assert simplify(Sqrt(g.compute_determinant())) == simplify(
        VolumeElement(domain, point).expand_geometry()
    )


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_the_volume_element_is_the_determinants_magnitude_on_a_cell(cell, gdim, lagrange_element):
    """Where the map is square the density is abs(det J), and the sign is kept."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    (entity,) = domain.cell_types
    tdim = entity.topological_dimension
    if gdim != tdim:
        pytest.skip("A manifold, so the determinant does not exist.")
    point = Point([RealScalar(0.25)] * tdim, EntityDomain(entity))

    assert simplify(abs(JacobianDeterminant(domain, point).expand_geometry())) == simplify(
        VolumeElement(domain, point).expand_geometry()
    )


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_a_manifold_map_has_no_determinant(cell, gdim, lagrange_element):
    """A non-square matrix has none, and the density is what was wanted anyway."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    (entity,) = domain.cell_types
    tdim = entity.topological_dimension
    if gdim == tdim:
        pytest.skip("Not a manifold, so the determinant exists.")
    point = Point([RealScalar(0.25)] * tdim, EntityDomain(entity))

    with pytest.raises(ValueError, match="has no determinant"):
        JacobianDeterminant(domain, point).value_shape
    with pytest.raises(ValueError, match="has no determinant"):
        JacobianDeterminant(domain, point).expand_geometry()

    assert VolumeElement(domain, point).value_shape == ()


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_the_projector_is_square_in_the_ambient_dimension(cell, gdim, lagrange_element):
    """It acts on the ambient coordinates, keeping the tangent directions."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))

    assert TangentialProjector(domain).value_shape == (gdim, gdim)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_the_projector_is_the_identity_only_on_a_square_map(cell, gdim, lagrange_element):
    """Where a cell can move in every ambient direction there is nothing to project."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    (entity,) = domain.cell_types
    tdim = entity.topological_dimension
    point = Point([RealScalar(0.25)] * tdim, EntityDomain(entity))

    expanded = TangentialProjector(domain, point).expand_geometry()

    if gdim == tdim:
        assert expanded == Identity(gdim)
    else:
        assert not isinstance(expanded, Identity)
        assert expanded.value_shape == (gdim, gdim)


def test_the_projector_is_idempotent(lagrange_element):
    """Projecting an already projected vector changes nothing."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (3,)))
    projector = TangentialProjector(domain)

    assert simplify(projector @ projector) == projector


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_the_normal_exists_exactly_at_codimension_one(cell, gdim, lagrange_element):
    """One direction is orthogonal to the tangent space only when one is left over."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    (entity,) = domain.cell_types
    tdim = entity.topological_dimension

    if gdim == tdim + 1:
        assert UnitNormal(domain).value_shape == (gdim,)
    else:
        with pytest.raises(ValueError, match="codimension"):
            UnitNormal(domain).value_shape


def test_a_domain_of_no_topological_dimension_has_no_normal_to_compute(lagrange_element):
    """The normal to a point is a sign, which is a convention rather than a value."""
    domain = parametrized_domain(lagrange_element("point", 1, (1,)))

    with pytest.raises(NotImplementedError, match="sign"):
        UnitNormal(domain).value_shape


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_the_spatial_coordinate_is_the_maps_value(cell, gdim, lagrange_element):
    """X is phi(X), so expanding it is asking the map what it gives."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    (entity,) = domain.cell_types
    tdim = entity.topological_dimension
    point = Point([RealScalar(0.25)] * tdim, EntityDomain(entity))

    coordinates = SpatialCoordinate(domain, point)

    assert coordinates.value_shape == (gdim,)
    assert coordinates.expand_geometry() == domain.parametrization(entity).value(point)


def test_both_ways_of_indexing_a_coordinate_check_their_range(lagrange_element):
    """`x[i]` and `x.component(i)` are the same path, so they agree.

    They used not to: one raised and the other gave a nonsense node.
    """
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    coordinates = SpatialCoordinate(domain)

    assert coordinates[1] == coordinates.component(1)
    for out_of_range in [-1, 2]:
        with pytest.raises(IndexError, match="out of range"):
            coordinates[out_of_range]
        with pytest.raises(IndexError, match="out of range"):
            coordinates.component(out_of_range)


def test_one_coordinate_stays_a_node_of_its_own(lagrange_element):
    """`sin(x[0])` should stay readable until geometry is expanded."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    (entity,) = domain.cell_types
    point = Point([RealScalar(0.25)] * 2, EntityDomain(entity))

    single = SpatialCoordinate(domain, point)[1]

    assert isinstance(single, SingleSpatialCoordinate)
    assert single.value_shape == ()
    assert single.expand_geometry() == SpatialCoordinate(domain, point).expand_geometry().component(
        1
    )
    with pytest.raises(ValueError, match="scalar"):
        single.component(0)


def test_one_coordinate_keeps_which_coordinate_it_is_when_told_a_point(lagrange_element):
    """It carries more than a domain and a point, so it rebuilds itself."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    (cell,) = domain.cell_types
    variable = create_variable(domain).to_entity_coordinates(cell)

    told = SingleSpatialCoordinate(domain, 1).reconstruct_with_variable(variable)

    assert told.point == variable
    assert told == SingleSpatialCoordinate(domain, 1, variable)


def test_a_pushed_forward_point_is_the_spatial_coordinate_as_a_point(lagrange_element):
    """The two are one computation; the difference is only the type."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (3,)))
    (entity,) = domain.cell_types
    point = Point([RealScalar(0.25)] * 2, EntityDomain(entity))

    pushed = PushedForwardPoint(point, domain).expand_geometry()
    coordinates = SpatialCoordinate(domain, point).expand_geometry()

    assert isinstance(pushed, Point)
    assert [pushed.component(i) for i in range(3)] == [coordinates.component(i) for i in range(3)]
