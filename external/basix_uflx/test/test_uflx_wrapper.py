# Copyright (c) 2024-2027 Matthew Scroggs
# FEniCS Project
# SPDX-License-Identifier: MIT
"""Tests for Basix's UFLx wrapper."""

import basix
import numpy as np
import pytest

import basix_uflx


@pytest.mark.parametrize(
    ("input_args", "input_kwargs"),
    [
        (("Lagrange", "triangle", 2), {}),
        (("Lagrange", basix.CellType.triangle, 2), {}),
        ((basix.ElementFamily.P, basix.CellType.triangle, 2), {}),
        ((basix.ElementFamily.P, "triangle", 2), {}),
        (("Lagrange", "triangle", 1), {"shape": (2,)}),
        (("Lagrange", "triangle", 2), {"shape": (2,)}),
        (("Lagrange", basix.CellType.triangle, 2), {"shape": (2,)}),
        ((basix.ElementFamily.P, basix.CellType.triangle, 2), {"shape": (2,)}),
        ((basix.ElementFamily.P, "triangle", 2), {"shape": (2,)}),
        (("Lagrange", "triangle", 2), {"shape": (2, 2)}),
        (("Lagrange", basix.CellType.triangle, 2), {"shape": (2, 2)}),
        ((basix.ElementFamily.P, basix.CellType.triangle, 2), {"shape": (2, 2)}),
        ((basix.ElementFamily.P, "triangle", 2), {"shape": (2, 2)}),
    ],
)
def test_element(input_args, input_kwargs):
    """Test element initialisation and tabulate shape."""
    e = basix_uflx.element(*input_args, **input_kwargs)
    table = e.tabulate(0, np.array([[0, 0]]))
    assert table.shape == (1, 1, e.dim, e.reference_value_size)


@pytest.mark.parametrize(
    "cell",
    [
        basix.CellType.triangle,
        basix.CellType.quadrilateral,
        basix.CellType.tetrahedron,
        basix.CellType.prism,
    ],
)
@pytest.mark.parametrize("degree", [1, 3, 6])
@pytest.mark.parametrize("shape", [(), (1,), (2,), (3,), (5,), (2, 2), (3, 3), (4, 1), (5, 1, 7)])
def test_quadrature_element(cell, degree, shape):
    """Test quadrature elements."""
    scalar_e = basix_uflx.quadrature_element(cell, (), degree=degree)
    e = basix_uflx.quadrature_element(cell, shape, degree=degree)

    size = 1
    for i in shape:
        size *= i

    assert e.reference_value_size == scalar_e.reference_value_size * size
    assert e.dim == scalar_e.dim * size


def test_wrap_element():
    """Test direct wrapping of a Basix element."""
    e = basix.create_element(basix.ElementFamily.P, basix.CellType.triangle, 1)
    basix_uflx.wrap_element(e)


def test_dof_ordering():
    """Test non-standard DOF ordering."""
    e = basix_uflx.element(basix.ElementFamily.P, basix.CellType.triangle, 1)

    e_reordered = basix_uflx.element(
        basix.ElementFamily.P, basix.CellType.triangle, 1, dof_ordering=[1, 2, 0]
    )

    points = np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [0.5, 0.5]])

    table = e.tabulate(0, points)[0]
    table2 = e_reordered.tabulate(0, points)[0]

    for i, j in enumerate([1, 2, 0]):
        assert np.allclose(table[:, i], table2[:, j])
