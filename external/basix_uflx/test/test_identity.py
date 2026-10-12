"""Test that Basix cells and elements are identified by their init_args."""

import os
import subprocess
import sys
import textwrap

import basix
import numpy as np
import pytest
import uflx
from uflx.graphs.graphs import fingerprint

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
        basix_uflx.element("Lagrange", "triangle", 1, shape=(2, 2), symmetry=True),
        basix_uflx.mixed_element(
            [basix_uflx.element("Lagrange", "triangle", 2), basix_uflx.real_element("triangle", ())]
        ),
        _quadrature(),
        basix_uflx.real_element("triangle", (2,)),
    ]


@pytest.mark.parametrize("index", range(len(_values())))
def test_built_twice_is_one_value(index):
    """Values built from equal arguments are one object."""
    assert _values()[index] is _values()[index]


@pytest.mark.parametrize("index", range(len(_values())))
def test_rebuilt_from_its_arguments(index):
    """A value can be rebuilt from its init_args."""
    value = _values()[index]
    assert type(value)(*value.init_args) is value


def test_distinct_values_are_unequal():
    """Values built differently are unequal."""
    values = _values()
    for i, a in enumerate(values):
        for b in values[i + 1 :]:
            assert a != b


def test_quadrature_elements_are_equal_exactly():
    """Quadrature elements are compared exactly and not to a tolerance."""
    points, weights = basix.make_quadrature(basix.CellType.triangle, 2)
    identity = uflx.maps.IdentityReferenceMap()
    a = basix_uflx.QuadratureElement(basix.CellType.triangle, points, weights, identity)
    b = basix_uflx.QuadratureElement(basix.CellType.triangle, points + 1e-15, weights, identity)
    assert np.allclose(points, points + 1e-15)
    assert a != b
    assert a == basix_uflx.QuadratureElement(
        basix.CellType.triangle, points.copy(), weights, identity
    )


def test_quadrature_points_cannot_change():
    """Changing the input array does not change the element."""
    points, weights = basix.make_quadrature(basix.CellType.triangle, 2)
    identity = uflx.maps.IdentityReferenceMap()
    e = basix_uflx.QuadratureElement(basix.CellType.triangle, points, weights, identity)
    before = hash(e)
    points[0, 0] = 1.0
    assert hash(e) == before
    assert isinstance(e.init_args[1], tuple)


def test_blocked_ordering_is_part_of_the_element():
    """Blocked elements with different orderings are unequal."""
    sub = basix_uflx.element("Lagrange", "triangle", 1)
    a = BlockedElement(sub, (2,), None, BlockedOrdering.xyzxyz)
    b = BlockedElement(sub, (2,), None, BlockedOrdering.xxyyzz)
    assert a != b


def test_blocked_symmetry_false_is_no_symmetry():
    """Symmetry False builds the same element as no symmetry."""
    sub = basix_uflx.element("Lagrange", "triangle", 1)
    a = basix_uflx.blocked_element(sub, (2, 2), symmetry=False)
    b = basix_uflx.blocked_element(sub, (2, 2))
    assert a is b


def test_every_process_agrees_on_an_elements_hash_whatever_its_hash_seed():
    """Guards against a wrapped Basix object hashing differently on two MPI ranks."""
    script = textwrap.dedent("""
        import basix_uflx

        elements = (
            basix_uflx.element("Lagrange", "triangle", 2, shape=(2,)),
            basix_uflx.quadrature_element("triangle", degree=2),
            basix_uflx.real_element("triangle", ()),
        )
        print([hash(e) for e in elements])
    """)
    printed = {
        seed: subprocess.run(
            [sys.executable, "-c", script],
            env={**os.environ, "PYTHONHASHSEED": seed},
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        for seed in ("1", "2")
    }
    assert printed["1"] == printed["2"]


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


def test_an_element_prints_as_its_description():
    """Guards against an element printing the address of the Basix object it wraps."""
    printed = repr(basix_uflx.element("Lagrange", "triangle", 2))
    assert "ElementFamily.P" in printed
    assert " at 0x" not in printed


def test_fingerprints_tell_apart_what_equality_tells_apart():
    """Blocked elements differing only in ordering have different digests."""
    sub = basix_uflx.element("Lagrange", "triangle", 1)
    a = BlockedElement(sub, (2,), None, BlockedOrdering.xyzxyz)
    b = BlockedElement(sub, (2,), None, BlockedOrdering.xxyyzz)
    assert fingerprint(a) != fingerprint(b)
