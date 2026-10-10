"""Test the value shapes of points mapped between reference and physical cells."""

import pytest

from uflx import parametrized_domain
from uflx.domains import RD, EntityDomain
from uflx.expressions import RealScalar
from uflx.geometry import ImagePoint, PreimagePoint
from uflx.points import Point

cells_and_gdims = [
    ("interval", 1),
    ("interval", 2),
    ("interval", 3),
    ("triangle", 2),
    ("triangle", 3),
    ("quadrilateral", 3),
    ("tetrahedron", 3),
]


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_reference_to_physical_shape(cell, gdim, lagrange_element):
    """A reference point mapped to a physical cell has gdim components."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cell_types[0].topological_dimension

    point = ImagePoint(
        Point([RealScalar(0.1)] * tdim, EntityDomain(domain.cell_types[0])),
        domain,
    )
    assert point.value_shape == (gdim,)
    assert point.dim == gdim
    assert point.expand_geometry().value_shape == (gdim,)


@pytest.mark.parametrize(("cell", "gdim"), cells_and_gdims)
def test_physical_to_reference_shape(cell, gdim, lagrange_element):
    """A physical point mapped to the reference cell has tdim components."""
    domain = parametrized_domain(lagrange_element(cell, 1, (gdim,)))
    tdim = domain.cell_types[0].topological_dimension

    (cell_type,) = domain.cell_types
    point = PreimagePoint(
        Point([RealScalar(0.1)] * gdim, RD(gdim)), domain.parametrization(cell_type)
    )
    assert point.value_shape == (tdim,)
    assert point.dim == tdim
