"""Test geometry."""

import pytest

from uflx import Coefficient, TestFunction, dx, function_space, inner, parametrized_domain
from uflx.algorithms import pull_back_to_entity
from uflx.basis_functions import EvaluatedBasisFunction
from uflx.domains import RD, EntityDomain
from uflx.expressions import RealScalar
from uflx.functions import create_variable
from uflx.geometry import (
    AbstractJacobian,
    Jacobian,
    JacobianInverse,
    JacobianInverseTranspose,
    JacobianTranspose,
    PulledBackPoint,
    PushedForwardPoint,
    expand_geometry,
)
from uflx.graphs import as_graph
from uflx.integrals import Integral
from uflx.points import Point
from uflx.tensors import FlattenedTensorMap

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
    jacobians = [n for n in as_graph(pulled) if isinstance(n, AbstractJacobian)]
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

    assert not any(isinstance(n, AbstractJacobian) for n in as_graph(expanded))


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
