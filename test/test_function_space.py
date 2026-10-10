# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Test function spaces."""

import pytest

from uflx import function_space, parametrized_domain
from uflx.domains import AbstractParametrizedDomain
from uflx.functions import Coefficient, TestFunction


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
    assert len(space.elements) == len(space.domain.cell_types) == 1


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
    assert len(space.elements) == len(space.domain.cell_types) == 2


def test_element_must_live_on_a_cell_of_the_domain(lagrange_element):
    """An element on a cell the domain does not have is rejected."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (3,)))

    with pytest.raises(ValueError, match="not defined on a cell of its domain"):
        function_space(domain, lagrange_element("quadrilateral", 1))


def test_separately_built_spaces_are_equal(lagrange_element):
    """Two spaces of the same elements on the same domain are the same space."""

    def build(degree: int):
        domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
        return function_space(domain, lagrange_element("triangle", degree))

    V, W = build(1), build(1)

    assert V is not W
    assert V == W
    assert hash(V) == hash(W)
    assert len({V, W}) == 1
    assert V != build(2)


def test_spaces_on_different_domains_differ(lagrange_element):
    """A space is not equal to the same elements on a different domain."""
    flat = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    surface = parametrized_domain(lagrange_element("triangle", 1, (3,)))

    assert function_space(flat, lagrange_element("triangle", 1)) != function_space(
        surface, lagrange_element("triangle", 1)
    )


def test_functions_on_equal_spaces_are_equal(lagrange_element):
    """Space equality reaches the functions, so the graph deduplicates them."""

    def build():
        domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
        return function_space(domain, lagrange_element("triangle", 1))

    V, W = build(), build()

    assert Coefficient(V, "f") == Coefficient(W, "f")
    assert TestFunction(V) == TestFunction(W)
    # An unlabelled coefficient still gets its own generated label.
    assert Coefficient(V) != Coefficient(W)
