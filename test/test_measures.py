"""Test measures."""

import pytest

from uflx import Coefficient, TestFunction, dx, function_space, inner, parametrized_domain
from uflx.algorithms import pull_back_to_entity
from uflx.domains import EntityDomain
from uflx.geometry import VolumeElement
from uflx.integrals import AbstractMeasure, Integral, Measure


@pytest.fixture
def mesh(lagrange_element):
    """A triangle mesh in two dimensions."""
    return parametrized_domain(lagrange_element("triangle", 1, (2,)))


@pytest.fixture
def space(mesh, lagrange_element):
    """A P1 space on that mesh."""
    return function_space(mesh, lagrange_element("triangle", 1))


def test_a_measure_names_the_domain_it_integrates_over(mesh):
    """Which domain is meant is said, never inferred."""
    assert dx(mesh).domain == mesh
    assert dx(mesh) == Measure(mesh)


def test_a_measures_density_is_its_domains_volume_element(mesh):
    """What a function is integrated against, and what a pull back picks up."""
    assert dx(mesh).density == VolumeElement(mesh)


def test_a_measure_is_retargeted_onto_another_domain(mesh):
    """Retargeting leaves the measure it came from alone."""
    (cell,) = mesh.cell_types
    entity = EntityDomain(cell)

    retargeted = dx(mesh).with_domain(entity)

    assert retargeted.domain == entity
    assert retargeted.density == VolumeElement(entity)
    assert dx(mesh).domain == mesh


def test_measures_compare_and_hash_by_value(mesh, lagrange_element):
    """A measure is a description, so two of the same description are one.

    Integral equality reaches for the measure, and a measure is a
    dictionary key while an expression is being rewritten.
    """
    other = parametrized_domain(lagrange_element("triangle", 2, (2,)))

    assert dx(mesh) == dx(mesh)
    assert hash(dx(mesh)) == hash(dx(mesh))
    assert dx(mesh) != dx(other)
    assert len({dx(mesh), dx(mesh), dx(other)}) == 2


def test_a_measure_rebuilds_from_its_arguments(mesh):
    """Its init_args are what reconstructing a graph node hands back."""
    measure = dx(mesh)

    assert Measure(*measure.init_args) == measure


def test_a_measure_has_no_successors(mesh):
    """Its density is not part of an integral's expression.

    The density is told where it is evaluated only when a pull back
    multiplies it into the integrand, so walking an integral's graph must
    not reach for it.
    """
    assert dx(mesh).successors == set()


def test_a_measure_shows_its_domain(mesh):
    """It says what it is of, not just its class name."""
    assert repr(dx(mesh)).startswith("Measure(")
    assert repr(mesh) in repr(dx(mesh))


def test_a_measure_cannot_be_mutated(mesh):
    """Its hash comes from its arguments, so moving them would lose it in a dict."""
    measure = dx(mesh)
    found_by = {measure: "here"}

    with pytest.raises(AttributeError):
        setattr(measure, "domain", EntityDomain(mesh.cell_types[0]))

    assert measure in found_by


def test_every_measure_is_a_measure():
    """The domain, the density and the retargeting are the whole interface."""
    for name in ["domain", "density", "with_domain", "init_args"]:
        assert name in AbstractMeasure.__abstractmethods__


def test_an_integral_takes_its_domain_from_its_measure(mesh, space):
    """The integral asks the measure rather than reading the integrand."""
    integral = inner(Coefficient(space), TestFunction(space)) * dx(mesh)

    assert isinstance(integral, Integral)
    assert integral.domain == mesh
    assert integral.variable.domain == mesh


def test_integrating_a_function_over_a_domain_it_is_not_on_is_rejected(space, lagrange_element):
    """A function reaches another domain through a map, not by being integrated there."""
    elsewhere = parametrized_domain(lagrange_element("triangle", 2, (2,)))

    with pytest.raises(ValueError, match="Cannot integrate a function on"):
        inner(Coefficient(space), TestFunction(space)) * dx(elsewhere)


def test_an_integrand_spanning_two_domains_is_rejected(mesh, space, lagrange_element):
    """Two subdomains meeting at an interface need a map between them, not one dx."""
    other_mesh = parametrized_domain(lagrange_element("triangle", 2, (2,)))
    other_space = function_space(other_mesh, lagrange_element("triangle", 1))

    with pytest.raises(ValueError, match="Cannot integrate a function on"):
        inner(Coefficient(space), TestFunction(other_space)) * dx(mesh)


def test_pulling_back_retargets_the_measure(mesh, space):
    """Change of variables moves the measure, not only the integrand."""
    (cell,) = mesh.cell_types
    integral = inner(Coefficient(space), TestFunction(space)) * dx(mesh)
    assert isinstance(integral, Integral)

    pulled = pull_back_to_entity(integral)

    assert isinstance(pulled, Integral)
    assert integral.measure.domain == mesh
    assert pulled.measure.domain == EntityDomain(cell)
    assert pulled.domain == EntityDomain(cell)
