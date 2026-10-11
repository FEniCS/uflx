"""Test the legacy UFL spellings, which build the native ones.

Each test says what UFL writes and what UFLx makes of it. The two that
refuse are the point of the module as much as the ones that work: a
legacy name that cannot be given its legacy meaning says so rather than
giving a form whose sign or whose sides are wrong.
"""

from typing import cast

import pytest
from conftest import Interval, LagrangeElement, Triangle
from test_facet_domains import PlacedTriangle

from uflx import (
    Coefficient,
    TestFunction,
    facets_of,
    function_space,
    grad,
    inner,
    parametrized_domain,
)
from uflx.algorithms import pull_back_to_entity
from uflx.domains import EntityDomain, LabelledDomain, labelled_domain
from uflx.geometry import AbstractGeometricQuantity, UnitNormal
from uflx.graphs import as_graph
from uflx.integrals import Integral, IntegralSum
from uflx.legacy import Cells, ExteriorFacets, dS, ds, dx, facet_normal
from uflx.restrictions import Restriction


@pytest.fixture
def mesh():
    """A mesh of placed triangles, whose facets can be included."""
    return parametrized_domain(LagrangeElement(PlacedTriangle(), 1, (2,)))


@pytest.fixture
def space(mesh):
    """A P1 space on it."""
    return function_space(mesh, LagrangeElement(PlacedTriangle(), 1))


def kernels(form) -> tuple[Integral, ...]:
    """The integrals a form pulls back to."""
    pulled = pull_back_to_entity(form)
    terms = pulled.terms if isinstance(pulled, IntegralSum) else (pulled,)
    assert all(isinstance(term, Integral) for term in terms)
    return cast(tuple[Integral, ...], terms)


def test_a_cell_measure_reads_its_domain_off_the_integrand(mesh, space):
    """``u * v * dx``, with the domain inferred rather than named."""
    u, v = Coefficient(space), TestFunction(space)

    form = u * v * dx

    assert isinstance(form, Integral)
    assert form.domain == mesh
    assert [term.domain for term in kernels(form)] == [EntityDomain(PlacedTriangle())]


def test_an_exterior_facet_measure_integrates_over_the_facets(mesh, space):
    """``u * v * ds`` is the one measure over the facets of the inferred domain."""
    u, v = Coefficient(space), TestFunction(space)

    form = u * v * ds

    assert isinstance(form, Integral)
    assert form.domain == facets_of(mesh, ExteriorFacets())
    assert [term.domain for term in kernels(form)] == [EntityDomain(Interval())] * 3


def test_a_field_in_a_facet_integral_is_restricted_for_you(mesh, space):
    """What the legacy spelling leaves out, written in: the field is read there."""
    u, v = Coefficient(space), TestFunction(space)

    form = u * v * ds

    assert isinstance(form, Integral)
    restrictions = [
        n for n in as_graph(form.integrand).ordered_nodes() if isinstance(n, Restriction)
    ]
    assert len(restrictions) == 1
    assert restrictions[0].domain == facets_of(mesh, ExteriorFacets())


def test_a_gradient_in_a_facet_integral_is_restricted_whole(mesh, space):
    """``inner(grad(u), n) * v * ds``, where wrapping the wrong subtree would be wrong.

    The restriction has to take the gradient whole. Wrapping ``u``
    instead would give the facet's tangential gradient and lose the
    derivative across the facet, and nothing would complain.
    """
    u, v = Coefficient(space), TestFunction(space)
    n = UnitNormal(facets_of(mesh, ExteriorFacets()))

    form = inner(grad(u), n) * v * ds
    assert isinstance(form, Integral)

    restricted = [
        type(node.field).__name__
        for node in as_graph(form.integrand).ordered_nodes()
        if isinstance(node, Restriction)
    ]
    assert sorted(restricted) == ["Grad", "TestFunction"]
    assert len(kernels(form)) == 3


def test_a_cell_integral_and_a_facet_one_add(mesh, space):
    """``u * v * dx + u * v * ds``, which is the first thing anyone writes."""
    u, v = Coefficient(space), TestFunction(space)

    form = u * v * dx + u * v * ds

    assert isinstance(form, IntegralSum)
    domains = [term.domain for term in kernels(form)]
    assert domains.count(EntityDomain(PlacedTriangle())) == 1
    assert domains.count(EntityDomain(Interval())) == 3


def test_a_marked_part_of_the_facets_is_another_domain(mesh, space):
    """``ds(3)``: a subdomain is a domain, labelled, and not an id on a measure."""
    u, v = Coefficient(space), TestFunction(space)

    marked = u * v * ds(3)
    unmarked = u * v * ds
    other = u * v * ds(4)
    assert isinstance(marked, Integral)
    assert isinstance(unmarked, Integral)
    assert isinstance(other, Integral)

    assert marked.domain == facets_of(mesh, ExteriorFacets(3))
    assert marked.domain != unmarked.domain
    assert marked.domain != other.domain
    assert len(kernels(marked)) == 3


def test_a_marked_tag_is_still_there_in_the_kernels(mesh, space):
    """What a label is for: a consumer can tell which facets a kernel is for.

    The label is on each piece, so splitting and pulling back keep it,
    where a label on the domain of facets alone would have been dropped
    at the split.
    """
    u, v = Coefficient(space), TestFunction(space)

    for term in kernels(u * v * ds(3)):
        read_on = [
            node.domain
            for node in as_graph(term.integrand).ordered_nodes()
            if isinstance(node, AbstractGeometricQuantity)
        ]
        assert all(isinstance(domain, LabelledDomain) for domain in read_on)
        labels = {domain.label for domain in read_on if isinstance(domain, LabelledDomain)}
        assert labels == {ExteriorFacets(3)}


def test_a_marked_part_of_the_cells_is_another_domain_too(mesh, space):
    """``dx(2)``: the same mechanism, so marked cells need no special case."""
    u, v = Coefficient(space), TestFunction(space)

    marked = u * v * dx(2)
    unmarked = u * v * dx
    assert isinstance(marked, Integral)
    assert isinstance(unmarked, Integral)

    assert marked.domain == labelled_domain(mesh, Cells(2))
    assert marked.domain != unmarked.domain
    assert len(kernels(marked)) == 1


def test_a_named_domain_is_used_rather_than_inferred(mesh, space):
    """``ds(domain=mesh)``, for an integrand that does not say."""
    u, v = Coefficient(space), TestFunction(space)

    form = u * v * ds(domain=mesh)
    assert isinstance(form, Integral)
    assert form.domain == facets_of(mesh, ExteriorFacets())


def test_an_integrand_that_does_not_say_which_domain_is_refused(space):
    """Inference can be wrong, so it refuses rather than guesses."""
    with pytest.raises(ValueError, match="says which domain"):
        UnitNormal(facets_of(parametrized_domain(LagrangeElement(PlacedTriangle(), 1, (2,))))) * ds

    other = parametrized_domain(LagrangeElement(PlacedTriangle(), 2, (2,)))
    elsewhere = function_space(other, LagrangeElement(PlacedTriangle(), 1))
    with pytest.raises(ValueError, match="on 2 domains"):
        inner(Coefficient(space), TestFunction(elsewhere)) * ds


def test_an_interior_facet_measure_says_what_it_would_need(mesh, space):
    """``dS`` cannot be built: an interface wants a parametrization per side."""
    u, v = Coefficient(space), TestFunction(space)

    with pytest.raises(NotImplementedError, match="parametrization per side"):
        u * v * dS


def test_a_facet_normal_says_why_it_has_no_sign(mesh):
    """``FacetNormal`` is the outward normal, and the sign is not there yet."""
    with pytest.raises(NotImplementedError, match="needs an orientation"):
        facet_normal(mesh)


def test_the_legacy_dx_says_it_is_not_the_native_one(mesh):
    """Importing it shadows ``uflx.dx``, whose argument is a domain."""
    with pytest.raises(ValueError, match="legacy dx"):
        dx(mesh)
    with pytest.raises(ValueError, match="legacy ds"):
        ds(mesh)


def test_a_domain_with_no_facet_geometry_still_refuses(lagrange_element):
    """The inclusion is the entity's to supply, and an unplaced cell has none."""
    plain = parametrized_domain(lagrange_element("triangle", 1, (2,)))
    space = function_space(plain, lagrange_element("triangle", 1))
    assert Triangle() in plain.cell_types

    with pytest.raises(NotImplementedError, match="does not say where its sub-entities are"):
        Coefficient(space) * TestFunction(space) * ds
