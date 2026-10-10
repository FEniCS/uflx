"""Test map algorithms."""

import pytest

from uflx import TestFunction, TrialFunction, dx, function_space, grad, inner, parametrized_domain
from uflx.algorithms import pull_back_to_entity
from uflx.domains import EntityDomain
from uflx.functions import (
    AbstractFunction,
    Coefficient,
    FiniteElementVariable,
)
from uflx.graphs import as_graph
from uflx.integrals import Integral


def test_mass_matrix(lagrange_element):
    """Test a mass matrix."""
    element = lagrange_element("triangle", 2)
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(u, v) * dx
    assert isinstance(form, Integral)

    pulled_form = pull_back_to_entity(form)
    assert isinstance(pulled_form, Integral)

    functions = [node for node in as_graph(form) if isinstance(node, AbstractFunction)]
    pulled_functions = [
        node for node in as_graph(pulled_form) if isinstance(node, AbstractFunction)
    ]

    assert len(functions) == 2
    assert len(pulled_functions) == 2

    for f in functions:
        assert isinstance(f, AbstractFunction) and not f.in_entity_coordinates
    for f in pulled_functions:
        assert isinstance(f, AbstractFunction) and f.in_entity_coordinates


def test_stuffness_matrix(lagrange_element):
    """Test a stiffness matrix."""
    element = lagrange_element("triangle", 2)
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    u = TrialFunction(space)
    v = TestFunction(space)
    form = inner(grad(u), grad(v)) * dx
    assert isinstance(form, Integral)

    pulled_form = pull_back_to_entity(form)
    assert isinstance(pulled_form, Integral)

    functions = [node for node in as_graph(form) if isinstance(node, AbstractFunction)]
    pulled_functions = [
        node for node in as_graph(pulled_form) if isinstance(node, AbstractFunction)
    ]

    assert len(functions) == 2
    assert len(pulled_functions) == 2

    for f in functions:
        assert isinstance(f, AbstractFunction) and not f.in_entity_coordinates
    for f in pulled_functions:
        assert isinstance(f, AbstractFunction) and f.in_entity_coordinates


def test_linear_form(lagrange_element):
    """Test a linear form."""
    element = lagrange_element("triangle", 2)
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    v = TestFunction(space)
    form = v * dx
    assert isinstance(form, Integral)

    pulled_form = pull_back_to_entity(form)
    assert isinstance(pulled_form, Integral)

    functions = [node for node in as_graph(form) if isinstance(node, AbstractFunction)]
    pulled_functions = [
        node for node in as_graph(pulled_form) if isinstance(node, AbstractFunction)
    ]

    assert len(functions) == 1
    assert len(pulled_functions) == 1

    for f in functions:
        assert isinstance(f, AbstractFunction) and not f.in_entity_coordinates
    for f in pulled_functions:
        assert isinstance(f, AbstractFunction) and f.in_entity_coordinates


def test_coefficient_mass_matrix_like_form(lagrange_element):
    """Test that a Coefficient pulls back like an Argument does."""
    element = lagrange_element("triangle", 2)
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    w = Coefficient(space)
    v = TestFunction(space)
    form = inner(w, v) * dx
    assert isinstance(form, Integral)

    pulled_form = pull_back_to_entity(form)
    assert isinstance(pulled_form, Integral)

    functions = [node for node in as_graph(form) if isinstance(node, AbstractFunction)]
    pulled_functions = [
        node for node in as_graph(pulled_form) if isinstance(node, AbstractFunction)
    ]

    assert len(functions) == 2
    assert len(pulled_functions) == 2

    for f in functions:
        assert isinstance(f, AbstractFunction) and not f.in_entity_coordinates
    for f in pulled_functions:
        assert isinstance(f, AbstractFunction) and f.in_entity_coordinates

    reference_coefficients = [
        f for f in pulled_functions if isinstance(f, Coefficient) and f.in_entity_coordinates
    ]
    assert len(reference_coefficients) == 1
    assert reference_coefficients[0].label == w.label


def test_coefficient_gradient_pulls_back(lagrange_element):
    """Test that grad(Coefficient) pulls back the same way grad(Argument) does."""
    element = lagrange_element("triangle", 2)
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    w = Coefficient(space)
    v = TestFunction(space)
    form = inner(grad(w), grad(v)) * dx
    assert isinstance(form, Integral)

    pulled_form = pull_back_to_entity(form)
    assert isinstance(pulled_form, Integral)

    pulled_functions = [
        node for node in as_graph(pulled_form) if isinstance(node, AbstractFunction)
    ]
    reference_coefficients = [
        f for f in pulled_functions if isinstance(f, Coefficient) and f.in_entity_coordinates
    ]
    assert len(reference_coefficients) == 1
    assert reference_coefficients[0].label == w.label


def test_distinct_coefficients_stay_distinguishable_after_pull_back(lagrange_element):
    """Test that two distinct Coefficients on the same space don't collapse together."""
    element = lagrange_element("triangle", 2)
    domain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    w1 = Coefficient(space)
    w2 = Coefficient(space)
    v = TestFunction(space)
    assert w1.label != w2.label

    form = inner(w1 + w2, v) * dx
    pulled_form = pull_back_to_entity(form)

    reference_coefficients = [
        node
        for node in as_graph(pulled_form)
        if isinstance(node, Coefficient) and node.in_entity_coordinates
    ]
    assert len(reference_coefficients) == 2
    assert {f.label for f in reference_coefficients} == {w1.label, w2.label}


def test_pull_back_rebases_the_space(lagrange_element):
    """Pulling a form back puts its functions on the entity's domain."""
    element = lagrange_element("triangle", 2)
    domain = parametrized_domain(lagrange_element("triangle", 1, (3,)))
    space = function_space(domain, element)
    form = inner(TrialFunction(space), TestFunction(space)) * dx

    pulled = pull_back_to_entity(form)
    functions = [node for node in as_graph(pulled) if isinstance(node, AbstractFunction)]

    assert len(functions) == 2
    for f in functions:
        assert f.function_space.domain == EntityDomain(element.cell)
        # Nothing is stored: the domain is what makes this true.
        assert f.in_entity_coordinates
        assert f.variable is not None and f.variable.domain == EntityDomain(element.cell)


def test_a_function_on_an_entity_domain_needs_no_flag(lagrange_element):
    """Building a function on an entity's domain is enough to put it in those coordinates."""
    element = lagrange_element("triangle", 2)
    mesh = parametrized_domain(lagrange_element("triangle", 1, (3,)))

    assert not Coefficient(function_space(mesh, element)).in_entity_coordinates
    assert Coefficient(function_space(EntityDomain(element.cell), element)).in_entity_coordinates


def test_pulling_a_variable_back_changes_only_where_it_lives(lagrange_element):
    """A pulled back variable keeps its label and takes the cell's domain."""
    element = lagrange_element("triangle", 2)
    mesh = parametrized_domain(lagrange_element("triangle", 1, (3,)))
    v = FiniteElementVariable(mesh, "x")

    pulled = v.to_entity_coordinates(element.cell)

    assert pulled.label == v.label
    assert pulled.domain == EntityDomain(element.cell)
    assert pulled != v
    assert len({pulled, v}) == 2


def test_a_variable_pulls_back_only_to_a_cell_of_its_domain(lagrange_element):
    """A facet is a sub-entity of a cell, not a cell of the mesh."""
    mesh = parametrized_domain(lagrange_element("triangle", 1, (3,)))
    (triangle,) = mesh.cell_types
    (facet,) = {f for f in triangle.sub_entities(1)}
    v = FiniteElementVariable(mesh, "x")

    assert v.to_entity_coordinates(triangle).domain == EntityDomain(triangle)
    with pytest.raises(AssertionError):
        v.to_entity_coordinates(facet)
