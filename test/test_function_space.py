# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Test function spaces."""

import pytest

from uflx import function_space, parametrized_domain
from uflx.domains import AbstractParametrizedDomain


@pytest.mark.parametrize(
    ("cell_name", "gdim"),
    [
        ("triangle", 2),
        ("triangle", 3),
        ("quadrilateral", 2),
        ("quadrilateral", 3),
        ("tetrahedron", 3),
    ],
)
def test_function_space(cell_name, gdim, lagrange_element):
    """Test function space with single cell."""
    domain = parametrized_domain(lagrange_element(cell_name, 1, (3,)))
    space = function_space(domain, lagrange_element(cell_name, 1))
    assert isinstance(space.domain, AbstractParametrizedDomain)
    assert len(space.elements) == len(space.domain.cells) == 1


@pytest.mark.parametrize("gdim", [2, 3])
def test_function_space_multiple_cells(gdim, lagrange_element):
    """Test function space with multiple cells."""
    domain = parametrized_domain(
        [
            lagrange_element("triangle", 1, (3,)),
            lagrange_element("quadrilateral", 1, (3,)),
        ]
    )
    space = function_space(
        domain,
        [lagrange_element("triangle", 2), lagrange_element("quadrilateral", 2)],
    )
    assert isinstance(space.domain, AbstractParametrizedDomain)
    assert len(space.elements) == len(space.domain.cells) == 2
