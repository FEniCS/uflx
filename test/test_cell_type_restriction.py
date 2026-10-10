"""Restricting a domain, a space and a function to one cell type."""

import pytest

from uflx import Coefficient, TestFunction, function_space, parametrized_domain
from uflx.domains import entity_domain


@pytest.fixture
def mixed_mesh(lagrange_element):
    """A mesh of triangles and quadrilaterals, so no one cell type is the cell type."""
    return parametrized_domain(
        [lagrange_element("triangle", 1, (2,)), lagrange_element("quadrilateral", 1, (2,))]
    )


def test_a_mixed_domain_restricts_to_each_of_its_cell_types(mixed_mesh):
    """Each restriction has one cell type and keeps that type's map."""
    for cell in mixed_mesh.cell_types:
        restricted = mixed_mesh.restricted_to(cell)

        assert restricted.cell_types == (cell,)
        assert restricted.geometric_dimension == mixed_mesh.geometric_dimension
        # The map keeps its element. The map object itself differs, because a
        # map is indexed against the domain whose coordinate dofs it sums.
        assert restricted.parametrization_element(cell) == mixed_mesh.parametrization_element(cell)


def test_restricting_a_domain_of_one_cell_type_changes_nothing(lagrange_element):
    """The restriction of a domain that is already restricted is itself."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    (cell,) = domain.cell_types

    assert domain.restricted_to(cell) == domain


def test_restricting_a_domain_to_a_cell_it_lacks_is_rejected(mixed_mesh, lagrange_element):
    """Asking for a cell type the domain does not have is a caller mistake."""
    (interval,) = parametrized_domain(lagrange_element("interval", 1, (1,))).cell_types

    with pytest.raises(ValueError, match="not a cell type"):
        mixed_mesh.restricted_to(interval)


def test_an_entity_domain_restricts_to_itself(lagrange_element):
    """An entity domain has one cell type already, which is its entity."""
    (triangle,) = parametrized_domain(lagrange_element("triangle", 1, (2,))).cell_types
    (interval,) = parametrized_domain(lagrange_element("interval", 1, (1,))).cell_types
    domain = entity_domain(triangle)

    assert domain.restricted_to(triangle) is domain

    with pytest.raises(ValueError, match="not the entity"):
        domain.restricted_to(interval)


def test_a_space_restricts_to_the_elements_on_one_cell(mixed_mesh, lagrange_element):
    """A space holding one element per cell type keeps the one for that cell."""
    elements = {cell: lagrange_element(cell.name, 1) for cell in mixed_mesh.cell_types}
    space = function_space(mixed_mesh, list(elements.values()))

    for cell, element in elements.items():
        restricted = space.restricted_to(cell)

        assert restricted.elements == (element,)
        assert restricted.domain == mixed_mesh.restricted_to(cell)


def test_a_space_keeps_every_element_on_the_cell_it_restricts_to(lagrange_element):
    """Several elements on one cell all survive, since the cell is what selects."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    (cell,) = domain.cell_types
    elements = (lagrange_element("triangle", 0), lagrange_element("triangle", 1))
    space = function_space(domain, elements)

    assert space.restricted_to(cell).elements == elements


def test_a_space_with_no_element_on_a_cell_is_rejected(mixed_mesh, lagrange_element):
    """A space cannot be restricted to a cell it has no basis on."""
    space = function_space(mixed_mesh, lagrange_element("triangle", 1))
    (quad,) = (c for c in mixed_mesh.cell_types if c != space.elements[0].cell)

    with pytest.raises(ValueError, match="no element on"):
        space.restricted_to(quad)


def test_a_function_restricts_onto_the_restricted_space(mixed_mesh, lagrange_element):
    """A function follows its space, keeping what identifies it."""
    space = function_space(
        mixed_mesh, [lagrange_element("triangle", 1), lagrange_element("quadrilateral", 1)]
    )
    coefficient = Coefficient(space)
    argument = TestFunction(space)

    for cell in mixed_mesh.cell_types:
        restricted = coefficient.restricted_to(cell)
        assert isinstance(restricted, Coefficient)
        assert restricted.function_space == space.restricted_to(cell)
        assert restricted.label == coefficient.label

        assert isinstance(argument.restricted_to(cell), TestFunction)


def test_restricting_a_function_that_does_not_move_gives_it_back(lagrange_element):
    """A function on a space of one cell type is its own restriction."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    (cell,) = domain.cell_types
    coefficient = Coefficient(function_space(domain, lagrange_element("triangle", 1)))

    assert coefficient.restricted_to(cell) is coefficient
