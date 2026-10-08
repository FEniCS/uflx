"""Test equality of points mapped between reference and physical cells."""

from uflx import coordinate_element
from uflx.expressions import RealScalar
from uflx.geometry import PhysicalToReference, ReferenceToPhysical
from uflx.points import Point


def test_physical_to_reference(lagrange_element):
    """A point mapped to the reference cell equals itself, and is on the reference cell."""
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    x = Point([RealScalar(0.1), RealScalar(0.2)])

    point = PhysicalToReference(x, domain)
    assert point == point
    assert point == PhysicalToReference(x, domain)
    assert hash(point) == hash(PhysicalToReference(x, domain))
    assert point.is_reference


def test_point_maps_differ(lagrange_element):
    """The two maps of the same point are different points."""
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    x = Point([RealScalar(0.1), RealScalar(0.2)], is_reference=True)

    assert PhysicalToReference(x, domain) != ReferenceToPhysical(x, domain)
    assert ReferenceToPhysical(x, domain) != PhysicalToReference(x, domain)
    assert len({PhysicalToReference(x, domain), ReferenceToPhysical(x, domain)}) == 2
