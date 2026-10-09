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


def test_element_must_live_on_a_cell_of_the_domain(lagrange_element):
    """An element on a cell the domain does not have is rejected."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (3,)))

    with pytest.raises(ValueError, match="not defined on a cell of its domain"):
        function_space(domain, lagrange_element("quadrilateral", 1))
