"""Test geometry."""

import pytest

from uflx import parametrized_domain
from uflx.basis_functions import EvaluatedBasisFunction
from uflx.domains import RD, EntityDomain
from uflx.expressions import RealScalar
from uflx.geometry import (
    Jacobian,
    JacobianInverse,
    JacobianInverseTranspose,
    JacobianTranspose,
    PulledBackPoint,
    PushedForwardPoint,
)
from uflx.graphs import as_graph
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
    assert pushed.domain == domain

    expanded = pushed.expand_geometry()
    assert isinstance(expanded, Point)
    assert expanded.domain == RD(gdim)
    assert not expanded.in_entity_coordinates


def test_pulled_back_point_lands_in_entity_coordinates(lagrange_element):
    """A pulled back point lies in the entity's coordinates, not on the mesh."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (3,)))
    (entity,) = domain.cell_types

    pulled = PulledBackPoint(Point([RealScalar(1.0)] * 3, RD(3)), domain)

    assert pulled.domain == EntityDomain(entity)
    assert pulled.in_entity_coordinates
    assert pulled.parametrized_domain == domain


def test_the_two_mapped_points_do_not_collide(lagrange_element):
    """A pushed forward point is never equal to a pulled back one."""
    domain = parametrized_domain(lagrange_element("interval", 1, (1,)))
    (entity,) = domain.cell_types

    entity_point = Point([RealScalar(0.25)], EntityDomain(entity))
    ambient_point = Point([RealScalar(0.25)], RD(1))

    pushed = PushedForwardPoint(entity_point, domain)
    pulled = PulledBackPoint(ambient_point, domain)

    assert pushed != pulled
    assert pulled != pushed
    assert len({pushed, pulled}) == 2


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_derivatives_are_entity_derivatives(cell, gdim, lagrange_element):
    """Test that the Jacobian differentiates in the tdim directions of the cell's coordinates."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cells[0].topological_dimension

    point = Point([RealScalar(1.0)] * tdim, EntityDomain(domain.cells[0]))
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
    tdim = domain.cells[0].topological_dimension
    point = Point([RealScalar(0.1)] * tdim, EntityDomain(domain.cells[0]))

    x_dofs = coordinate_dof_entries(PushedForwardPoint(point, domain).expand_geometry())
    assert len(x_dofs) > 0
    assert coordinate_dof_entries(Jacobian(domain, point).expand_geometry()) == x_dofs
