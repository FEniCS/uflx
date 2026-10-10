"""Restricting a domain, a space and a function to one cell type."""

import pytest

from uflx import Coefficient, TestFunction, dx, function_space, inner, parametrized_domain
from uflx.algorithms import pull_back_to_entity
from uflx.domains import EntityDomain, entity_domain
from uflx.geometry import AbstractGeometricQuantity, expand_geometry
from uflx.graphs import as_graph
from uflx.integrals import Integral, IntegralSum


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


@pytest.fixture
def mixed_space(mixed_mesh, lagrange_element):
    """A space with one element per cell type of the mixed mesh."""
    return function_space(
        mixed_mesh, [lagrange_element(cell.name, 1) for cell in mixed_mesh.cell_types]
    )


@pytest.fixture
def mixed_form(mixed_space, mixed_mesh):
    """A mass form over the mixed mesh, with no gradients in it."""
    return inner(Coefficient(mixed_space), TestFunction(mixed_space)) * dx(mixed_mesh)


def test_an_integral_knows_the_domain_it_is_over(mixed_form, mixed_mesh):
    """The domain is read off the functions in the integrand."""
    assert mixed_form.domain == mixed_mesh


def test_an_integral_restricts_to_one_cell_type(mixed_form, mixed_mesh):
    """Restricting takes the integrand, the dummy variable and the measure with it."""
    for cell in mixed_mesh.cell_types:
        restricted = mixed_form.restricted_to(cell)

        assert restricted.domain == mixed_mesh.restricted_to(cell)
        assert restricted.variable.domain == mixed_mesh.restricted_to(cell)
        assert restricted.measure == mixed_form.measure.with_domain(restricted.domain)


def test_splitting_by_cell_type(mixed_form, mixed_mesh, lagrange_element):
    """A domain of one cell type has nothing to split."""
    assert len(mixed_form.split_by_cell_type().terms) == len(mixed_mesh.cell_types)

    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, lagrange_element("triangle", 1))
    single = inner(Coefficient(space), TestFunction(space)) * dx(domain)
    assert isinstance(single, Integral)
    assert single.split_by_cell_type() is None


def test_integrals_add_to_a_sum_that_flattens(mixed_form, mixed_mesh):
    """A sum of integrals is flat, however it is built up."""
    first, second = (mixed_form.restricted_to(c) for c in mixed_mesh.cell_types)

    total = first + second
    assert isinstance(total, IntegralSum)
    assert total.terms == (first, second)
    assert len((total + first).terms) == 3
    assert len((first + total).terms) == 3
    assert total == first + second

    with pytest.raises(ValueError, match="empty sum"):
        IntegralSum(())


def test_a_mixed_mesh_form_pulls_back_to_one_integral_per_cell_type(mixed_form, mixed_mesh):
    """Each cell type gets its own integral, its own coordinates and its own map."""
    pulled = pull_back_to_entity(mixed_form)

    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == len(mixed_mesh.cell_types)

    for cell, term in zip(mixed_mesh.cell_types, pulled.terms, strict=True):
        assert term.variable.domain == EntityDomain(cell)
        jacobians = [n for n in as_graph(term) if isinstance(n, AbstractGeometricQuantity)]
        assert len(jacobians) > 0
        for jacobian in jacobians:
            assert jacobian.parametrization.source == EntityDomain(cell)


def test_a_pulled_back_mixed_form_expands_its_geometry(mixed_form):
    """Every term's geometry expands, since every term knows its cell."""
    expanded = expand_geometry(pull_back_to_entity(mixed_form))

    assert isinstance(expanded, IntegralSum)
    assert not any(isinstance(n, AbstractGeometricQuantity) for n in as_graph(expanded))


def test_a_form_over_one_cell_type_is_still_a_plain_integral(lagrange_element):
    """The path for a single cell type is untouched by the fan out."""
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, lagrange_element("triangle", 1))
    form = inner(Coefficient(space), TestFunction(space)) * dx(domain)

    assert isinstance(pull_back_to_entity(form), Integral)


def test_pulling_a_mixed_integral_back_directly_is_rejected(mixed_form):
    """There is no one set of coordinates to pull several cell types back to."""
    with pytest.raises(ValueError, match="Split it by cell type"):
        mixed_form.pull_back_to_entity({})
