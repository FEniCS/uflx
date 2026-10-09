# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Test finite elements."""


def test_lagrange_element(entity, lagrange_element):
    """Test Lagrange element properties."""
    element = lagrange_element(entity.name, 2)

    assert element.cell == entity
    assert element.entity_value_shape == ()
    assert element.lagrange_superdegree == 2


def test_describes_affine_map_on_simplex(lagrange_element, simplex):
    """A degree 1 element on a simplex describes an affine map, degree 2 does not."""
    assert lagrange_element(simplex.name, 1).describes_affine_map
    assert not lagrange_element(simplex.name, 2).describes_affine_map


def test_describes_affine_map_on_tensor_product_cell(lagrange_element):
    """Even at degree 1 a tensor-product cell's map is multilinear, not affine."""
    assert not lagrange_element("quadrilateral", 1).describes_affine_map
    assert not lagrange_element("hexahedron", 1).describes_affine_map
