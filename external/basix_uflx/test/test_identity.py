"""Test that Basix cells and elements are identified by their init_args."""

import basix
import numpy as np
import pytest

import basix_uflx
from basix_uflx import BasixCell
from basix_uflx.finite_element import BlockedElement, BlockedOrdering


def _quadrature():
    return basix_uflx.quadrature_element("triangle", degree=2)


def _values():
    return [
        BasixCell(basix.CellType.triangle),
        basix_uflx.element("Lagrange", "triangle", 2),
        basix_uflx.element("Lagrange", "triangle", 1, shape=(2,)),
        basix_uflx.element("Lagrange", "triangle", 1, shape=(2, 2)),
        basix_uflx.element("Lagrange", "triangle", 1, shape=(2, 2), symmetry=True),
        basix_uflx.mixed_element(
            [basix_uflx.element("Lagrange", "triangle", 2), basix_uflx.real_element("triangle", ())]
        ),
        _quadrature(),
        basix_uflx.real_element("triangle", (2,)),
        BlockedElement(
            basix_uflx.element("Lagrange", "triangle", 1), (2,), None, BlockedOrdering.xxyyzz
        ),
    ]


def test_distinct_values_are_unequal():
    """Values built differently are unequal."""
    values = _values()
    for i, a in enumerate(values):
        for b in values[i + 1 :]:
            assert a != b


def test_quadrature_elements_are_equal_exactly():
    """Quadrature elements are compared exactly and not to a tolerance."""
    points, weights = basix.make_quadrature(basix.CellType.triangle, 2)
    a = basix_uflx.quadrature_element("triangle", points=points, weights=weights)
    b = basix_uflx.quadrature_element("triangle", points=points + 1e-15, weights=weights)
    assert np.allclose(points, points + 1e-15)
    assert a != b
    assert a == basix_uflx.quadrature_element("triangle", points=points.copy(), weights=weights)


def test_blocked_symmetry_false_is_no_symmetry():
    """Symmetry False builds the same element as no symmetry."""
    sub = basix_uflx.element("Lagrange", "triangle", 1)
    a = basix_uflx.blocked_element(sub, (2, 2), symmetry=False)
    b = basix_uflx.blocked_element(sub, (2, 2))
    assert a is b


@pytest.mark.parametrize(
    "kwargs",
    [{}, {"dof_ordering": [0, 2, 1, 3, 4, 5]}, {"dtype": np.float32}],
    ids=["default", "dof_ordering", "float32"],
)
def test_a_wrapped_element_is_the_element_built_alike(kwargs):
    """Guards against wrap_element reading back a description that differs from element's."""
    wrapped = basix_uflx.wrap_element(
        basix.create_element(
            basix.ElementFamily.P,
            basix.CellType.triangle,
            2,
            basix.LagrangeVariant.gll_warped,
            **kwargs,
        )
    )
    assert wrapped is basix_uflx.element("Lagrange", "triangle", 2, **kwargs)


def test_a_custom_element_keeps_its_data_exactly():
    """Guards against a custom element losing the shape of an empty array or an entry."""
    p1 = basix.create_element(basix.ElementFamily.P, basix.CellType.triangle, 1)
    custom = basix.create_custom_element(
        p1.cell_type,
        p1.value_shape,
        p1.wcoeffs,
        p1.x,
        p1.M,
        p1.interpolation_nderivs,
        p1.map_type,
        p1.sobolev_space,
        p1.discontinuous,
        p1.embedded_subdegree,
        p1.embedded_superdegree,
        p1.polyset_type,
    )
    wrapped = basix_uflx.wrap_element(custom)
    points = np.array([[0.2, 0.3], [0.5, 0.1]])
    assert np.array_equal(wrapped.tabulate(1, points), custom.tabulate(1, points))
    assert basix_uflx.wrap_element(custom) is wrapped
