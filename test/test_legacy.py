"""Test the legacy UFL spellings, which build the native ones.

Each test says what UFL writes and what UFLx makes of it. The two that
refuse are the point of the module as much as the ones that work: a
legacy name that cannot be given its legacy meaning says so rather than
giving a form whose sign or whose sides are wrong.
"""

import subprocess
import sys
from typing import cast

import pytest
from conftest import Interval, LagrangeElement, Triangle
from test_facet_domains import PlacedTriangle
from test_interfaces import interior_facets_of

from uflx import (
    Coefficient,
    TestFunction,
    TrialFunction,
    function_space,
    grad,
    inner,
    parametrized_domain,
)
from uflx.algorithms import pull_back_to_entity
from uflx.domains import EntityDomain, LabelledDomain, labelled_domain
from uflx.expressions import AbstractExpression, RealScalar
from uflx.geometry import AbstractGeometricQuantity, OutwardNormal, UnitNormal, VolumeElement
from uflx.graphs import as_graph
from uflx.integrals import Integral, IntegralSum
from uflx.integrals import dx as measure_of
from uflx.interfaces import Interface
from uflx.legacy import (
    Cells,
    ExteriorFacets,
    InteriorFacets,
    SideMarker,
    add_plus_minus_restriction,
    avg,
    dS,
    ds,
    dx,
    facet_normal,
    jump,
    minus,
    plus,
    remove_plus_minus_restriction,
)
from uflx.parametrizations import sub_entities_of
from uflx.restrictions import Restriction
from uflx.restrictions import avg as native_avg
from uflx.restrictions import jump as native_jump


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
    assert form.domain == sub_entities_of(mesh, 1, ExteriorFacets())
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
    assert restrictions[0].domain == sub_entities_of(mesh, 1, ExteriorFacets())


def test_a_gradient_in_a_facet_integral_is_restricted_whole(mesh, space):
    """``inner(grad(u), n) * v * ds``, where wrapping the wrong subtree would be wrong.

    The restriction has to take the gradient whole. Wrapping ``u``
    instead would give the facet's tangential gradient and lose the
    derivative across the facet, and nothing would complain.
    """
    u, v = Coefficient(space), TestFunction(space)
    n = UnitNormal(sub_entities_of(mesh, 1, ExteriorFacets()))

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

    assert marked.domain == sub_entities_of(mesh, 1, ExteriorFacets(3))
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
    assert form.domain == sub_entities_of(mesh, 1, ExteriorFacets())


def test_an_integrand_that_does_not_say_which_domain_is_refused(space):
    """Inference can be wrong, so it refuses rather than guesses."""
    with pytest.raises(ValueError, match="says which domain"):
        UnitNormal(
            sub_entities_of(parametrized_domain(LagrangeElement(PlacedTriangle(), 1, (2,))), 1)
        ) * ds

    other = parametrized_domain(LagrangeElement(PlacedTriangle(), 2, (2,)))
    elsewhere = function_space(other, LagrangeElement(PlacedTriangle(), 1))
    with pytest.raises(ValueError, match="on 2 domains"):
        inner(Coefficient(space), TestFunction(elsewhere)) * ds


def restrictions_in(integral: Integral) -> set:
    """The fields an integral reads through a map, and the domain each is read through."""
    return {
        (node.field, node.domain)
        for node in as_graph(integral.integrand)
        if isinstance(node, Restriction)
    }


def test_an_interior_facet_measure_integrates_over_the_interface(mesh, space):
    """``dS`` is the one measure over the interface of the mesh with itself.

    ``jump`` and ``avg`` resolve to the native ones once the measure
    knows the interface, so the two spellings build the same integral.
    """
    u, v = Coefficient(space), TestFunction(space)
    gamma = interior_facets_of(mesh, InteriorFacets(None))

    integral = jump(u) * avg(v) * dS
    assert isinstance(integral, Integral)
    assert integral.domain == gamma
    native = native_jump(u, gamma) * native_avg(v, gamma) * measure_of(gamma)
    assert isinstance(native, Integral)
    assert integral.measure == native.measure
    assert restrictions_in(integral) == restrictions_in(native)
    assert restrictions_in(integral) == {
        (u, gamma.side("+")),
        (u, gamma.side("-")),
        (v, gamma.side("+")),
        (v, gamma.side("-")),
    }


def test_a_side_is_named_by_plus_and_minus(mesh, space):
    """``plus(u)`` stands where a restriction to the ``+`` side will."""
    u = Coefficient(space)
    gamma = interior_facets_of(mesh, InteriorFacets(None))

    marked = plus(u)
    assert isinstance(marked, SideMarker)
    assert marked.side == "+"
    assert marked.value_shape == u.value_shape

    integral = plus(u) * minus(u) * dS
    assert isinstance(integral, Integral)
    assert restrictions_in(integral) == {(u, gamma.side("+")), (u, gamma.side("-"))}

    with pytest.raises(ValueError, match="marked with a side already"):
        plus(minus(u))
    with pytest.raises(ValueError, match="A side is"):
        SideMarker(u, "left")


def test_adding_plus_minus_restriction_makes_an_expression_callable_with_a_side(mesh, space):
    """``u("+")``, UFL's own spelling, on request and for a class the user names.

    A patch of a core class, made at a line the user wrote naming the
    class, undone by another. Before the request, and in a fresh
    interpreter that only imported ``uflx``, there is no such call.
    """
    u = Coefficient(space)
    gamma = interior_facets_of(mesh, InteriorFacets(None))

    assert "__call__" not in vars(AbstractExpression)
    with pytest.raises(TypeError):
        u("+")  # pyrefly: ignore

    add_plus_minus_restriction(AbstractExpression)
    try:
        add_plus_minus_restriction(AbstractExpression)
        marked = u("+")  # pyrefly: ignore
        assert marked == plus(u)
        assert grad(u)("-") == minus(grad(u))  # pyrefly: ignore

        integral = u("+") * u("-") * dS  # pyrefly: ignore
        assert isinstance(integral, Integral)
        assert restrictions_in(integral) == {(u, gamma.side("+")), (u, gamma.side("-"))}

        with pytest.raises(ValueError, match="not with"):
            u(0.5)  # pyrefly: ignore
    finally:
        remove_plus_minus_restriction(AbstractExpression)
    assert "__call__" not in vars(AbstractExpression)
    remove_plus_minus_restriction(AbstractExpression)

    untouched = (
        "import uflx\n"
        "from uflx.expressions import AbstractExpression\n"
        "assert not any('__call__' in vars(c) for c in AbstractExpression.__mro__ "
        "if c is not object)"
    )
    fresh = subprocess.run([sys.executable, "-c", untouched], check=False)
    assert fresh.returncode == 0


def test_plus_minus_restriction_is_as_wide_as_the_class_named(space):
    """Naming ``Coefficient`` makes coefficients callable and leaves a gradient alone."""
    u = Coefficient(space)

    add_plus_minus_restriction(Coefficient)
    try:
        assert u("+") == plus(u)  # pyrefly: ignore
        with pytest.raises(TypeError):
            grad(u)("+")  # pyrefly: ignore
        with pytest.raises(TypeError):
            TestFunction(space)("+")  # pyrefly: ignore
    finally:
        remove_plus_minus_restriction(Coefficient)
    assert "__call__" not in vars(Coefficient)


def test_plus_minus_restriction_does_not_override_a_call_that_is_already_there():
    """Something else making expressions callable is a conflict, not a thing to replace.

    And a class that is not an expression has no side to be read of.
    """
    AbstractExpression.__call__ = lambda self, *args: args  # pyrefly: ignore
    try:
        with pytest.raises(ValueError, match="callable already"):
            add_plus_minus_restriction(AbstractExpression)
        with pytest.raises(ValueError, match="callable already"):
            add_plus_minus_restriction(Coefficient)
        remove_plus_minus_restriction(AbstractExpression)
        assert "__call__" in vars(AbstractExpression)
    finally:
        del AbstractExpression.__call__

    with pytest.raises(ValueError, match="not a kind of expression"):
        add_plus_minus_restriction(Integral)  # pyrefly: ignore


def test_an_unmarked_field_under_dS_is_refused(mesh, space):
    """A field reaches the interface only through a side, as UFL also insists."""
    u, v = Coefficient(space), TestFunction(space)
    with pytest.raises(ValueError, match="through a map"):
        u * v * dS
    with pytest.raises(ValueError, match="through a map"):
        jump(u) * v * dS


def test_a_facet_normal_is_read_off_the_facets_under_ds(mesh, space):
    """``FacetNormal(mesh)`` under ``ds`` is the outward normal of the facets."""
    u, v = Coefficient(space), TestFunction(space)
    n = facet_normal(mesh)
    assert isinstance(n, OutwardNormal)

    integral = inner(grad(u), n) * v * ds
    assert isinstance(integral, Integral)
    facets = sub_entities_of(mesh, 1, ExteriorFacets(None))
    assert restrictions_in(integral) == {(grad(u), facets), (v, facets)}
    normals = {node for node in as_graph(integral.integrand) if isinstance(node, OutwardNormal)}
    assert normals == {OutwardNormal(facets)}

    pulled = pull_back_to_entity(integral)
    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == 3


def test_a_facet_normal_is_read_off_a_side_under_dS(mesh, space):
    """``n("+")`` is the normal out of the ``+`` side, and ``jump(v, n)`` uses both."""
    u, v = Coefficient(space), TestFunction(space)
    n = facet_normal(mesh)
    gamma = interior_facets_of(mesh, InteriorFacets(None))

    assert jump(v, n) == plus(v) * plus(n) + minus(v) * minus(n)
    assert jump(grad(u), n) == inner(plus(grad(u)), plus(n)) + inner(minus(grad(u)), minus(n))

    integral = inner(avg(grad(u)), plus(n)) * jump(v) * dS
    assert isinstance(integral, Integral)
    normals = {node for node in as_graph(integral.integrand) if isinstance(node, OutwardNormal)}
    assert normals == {OutwardNormal(gamma.side("+"))}

    pulled = pull_back_to_entity(jump(grad(u), n) * jump(v) * dS)
    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == 9
    for pair, term in zip(gamma.pieces, pulled.terms, strict=True):
        assert isinstance(pair, Interface)
        normals = {node for node in as_graph(term.integrand) if isinstance(node, OutwardNormal)}
        assert {node.domain for node in normals} == {pair.side("+"), pair.side("-")}


def test_a_marked_part_of_the_interface_is_another_domain(mesh, space):
    """``dS(3)`` and ``dS(4)`` label both sides' facets, and the label reaches the pairs."""
    u, v = Coefficient(space), TestFunction(space)

    assert jump(u) * jump(v) * dS(3) != jump(u) * jump(v) * dS(4)

    pulled = pull_back_to_entity(jump(u) * jump(v) * dS(3))
    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == 9
    for term in pulled.terms:
        density = next(q for q in as_graph(term.integrand) if isinstance(q, VolumeElement))
        pair = density.domain
        assert isinstance(pair, Interface)
        for label in pair.labels:
            presentation = pair.presentation(label)
            assert isinstance(presentation, LabelledDomain)
            assert presentation.label == InteriorFacets(3)


def test_a_discontinuous_galerkin_form_builds_and_pulls_back(mesh, space):
    """Interior penalty for Poisson in UFL's spelling, with ``dx``, ``ds`` and ``dS``.

    One cell kernel, three exterior facet kernels per ``ds`` term and nine
    pairs per ``dS`` term.
    """
    u, v = TrialFunction(space), TestFunction(space)
    n = facet_normal(mesh)
    alpha = RealScalar(10.0)

    a = (
        inner(grad(u), grad(v)) * dx
        - inner(avg(grad(u)), jump(v, n)) * dS
        - inner(jump(u, n), avg(grad(v))) * dS
        + alpha * jump(u) * jump(v) * dS
        - inner(grad(u), n) * v * ds
        - u * inner(grad(v), n) * ds
        + alpha * u * v * ds
    )
    assert isinstance(a, IntegralSum)
    assert len(a.terms) == 7

    pulled = pull_back_to_entity(a)
    assert isinstance(pulled, IntegralSum)
    assert len(pulled.terms) == 1 + 3 * 9 + 3 * 3


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
