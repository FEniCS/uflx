"""Test geometry."""

import pytest

from uflx import coordinate_element
from uflx.basis_functions import EvaluatedBasisFunction
from uflx.expressions import RealScalar
from uflx.geometry import (
    CoordinateDofs,
    Jacobian,
    JacobianInverse,
    JacobianInverseTranspose,
    JacobianTranspose,
    ReferenceToPhysical,
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
    domain = coordinate_element(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cells[0].topological_dimension

    point = Point([RealScalar(1.0)] * tdim, is_reference=True)
    j = Jacobian(domain, point)
    assert j.value_shape == (gdim, tdim)

    mat = j.expand_geometry()
    assert mat.value_shape == (gdim, tdim)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_inverse_expand_geometry(cell, gdim, lagrange_element):
    """Test expansion of Jacobian inverse."""
    domain = coordinate_element(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cells[0].topological_dimension

    point = Point([RealScalar(1.0)] * tdim, is_reference=True)
    j = JacobianInverse(domain, point)
    assert j.value_shape == (tdim, gdim)

    mat = j.expand_geometry()
    assert mat.value_shape == (tdim, gdim)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_tranpose_expand_geometry(cell, gdim, lagrange_element):
    """Test expansion of Jacobian inverse transpose."""
    domain = coordinate_element(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cells[0].topological_dimension

    point = Point([RealScalar(1.0)] * tdim, is_reference=True)
    j = JacobianTranspose(domain, point)
    assert j.value_shape == (tdim, gdim)

    mat = j.expand_geometry()
    assert mat.value_shape == (tdim, gdim)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_inverse_transpose_expand_geometry(cell, gdim, lagrange_element):
    """Test expansion of Jacobian inverse transpose."""
    domain = coordinate_element(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cells[0].topological_dimension

    point = Point([RealScalar(1.0)] * tdim, is_reference=True)
    j = JacobianInverseTranspose(domain, point)
    assert j.value_shape == (gdim, tdim)

    mat = j.expand_geometry()
    assert mat.value_shape == (gdim, tdim)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_derivatives_are_reference_derivatives(cell, gdim, lagrange_element):
    """Test that the Jacobian differentiates in the tdim reference directions."""
    domain = coordinate_element(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cells[0].topological_dimension

    point = Point([RealScalar(1.0)] * tdim, is_reference=True)
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
        node.init_args[1]
        for node in as_graph(expression).ordered_nodes()
        if isinstance(node, FlattenedTensorMap)
    }


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_jacobian_coordinate_dofs(cell, gdim, lagrange_element):
    """Test that the Jacobian uses the coordinate DOFs of the map it differentiates."""
    domain = coordinate_element(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cells[0].topological_dimension
    point = Point([RealScalar(0.1)] * tdim, is_reference=True)

    x_dofs = coordinate_dof_entries(ReferenceToPhysical(point, domain).expand_geometry())
    assert len(x_dofs) > 0
    assert coordinate_dof_entries(Jacobian(domain, point).expand_geometry()) == x_dofs


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_geometry_reads_the_coordinate_dofs(cell, gdim, lagrange_element):
    """Test that the entries of the geometry are read from the coordinate DOFs of the domain."""
    domain = coordinate_element(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cells[0].topological_dimension
    point = Point([RealScalar(0.1)] * tdim, is_reference=True)

    for expression in [
        ReferenceToPhysical(point, domain).expand_geometry(),
        Jacobian(domain, point).expand_geometry(),
    ]:
        arrays = {
            node.array
            for node in as_graph(expression).ordered_nodes()
            if isinstance(node, FlattenedTensorMap)
        }
        assert arrays == {CoordinateDofs(domain)}
