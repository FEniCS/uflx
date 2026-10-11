# Copyright (C) 2026 Jack S. Hale
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Legacy UFL spellings, for code being moved onto UFLx.

Opt in by importing from here; nothing in `uflx` imports this module.
UFL writes a measure with no domain and reads the domain off the
integrand, and writes a field on a mesh directly in a facet integral:

    a = inner(grad(u), n) * v * ds

UFLx writes the same form by naming the domain and reading the field
there through its parametrization:

    facets = facets_of(mesh)
    a = inner(restrict(grad(u), facets), n) * restrict(v, facets) * dx(facets)

The second is what this module builds from the first, so old code runs
while it is being moved. The native spelling is the one to write.

What inference costs is the reason this is not the core language, and
INTEGRATION.md note 2 is the argument: an integral that guesses its
domain cannot be told that it guessed wrong, and where domain equality
is structural the guess can be wrong silently, two meshes built from the
same coordinate element comparing equal. So the measures here refuse
rather than choose when the integrand does not say, and the form that
names its domain can always be written instead.

The module is named for what it is rather than for a plural of what it
contains, bending the convention in DESIGN_CHOICES.md, since what it
contains is a vocabulary and not a kind of object.

A subdomain id is a domain and not an id on a measure, as INTEGRATION.md
note 3 says, so `ds(3)` integrates over the facets labelled 3 --
`LabelledDomain` per piece -- and the label is still there to be read off
a kernel afterwards. Which entities a consumer marked 3 is the
consumer's, as its cells are; UFLx records only that a selection was
meant and which one.
"""

from __future__ import annotations

from typing import Any, NamedTuple

from uflx.algorithms import replace
from uflx.domains import (
    AbstractParametrizedDomain,
    AbstractCellularDomain,
    UnionDomain,
    labelled_domain,
)
from uflx.expressions import AbstractExpression
from uflx.functions import AbstractFunction
from uflx.geometry import AbstractGeometricQuantity
from uflx.graphs import GraphNode, NodeOrder, as_graph
from uflx.integrals import Integral
from uflx.integrals import dx as measure_of
from uflx.parametrizations import facets_of
from uflx.restrictions import restrict


class Cells(NamedTuple):
    """What a legacy cell measure means, as a domain's label."""

    subdomain_id: Any = None


class ExteriorFacets(NamedTuple):
    """What a legacy exterior facet measure means, as a domain's label.

    UFLx does not interpret a label: which facets a domain has is
    external, as a mesh's cells are. This records what the consumer was
    told to pick -- the exterior facets, of one marked subdomain or of
    all of them -- so that `ds(1)` and `ds(2)` are different domains
    rather than one, which is what the label is for.
    """

    subdomain_id: Any = None


class InteriorFacets(NamedTuple):
    """What a legacy interior facet measure would mean, as a domain's label."""

    subdomain_id: Any = None


def _inferred_domain(integrand: AbstractExpression) -> AbstractParametrizedDomain:
    """Read the domain off the functions in an integrand.

    Args:
        integrand: The expression being integrated

    Returns:
        The domain its functions agree on

    Raises:
        ValueError: If no function says which domain, or they disagree, or
            it is not a domain made of cells
    """
    domains = [
        node.function_space.domain
        for node in as_graph(integrand).ordered_nodes()
        if isinstance(node, AbstractFunction) and not node.in_entity_coordinates
    ]
    distinct = list(dict.fromkeys(domains))
    if len(distinct) == 0:
        raise ValueError(
            "No function in this integrand says which domain it is on, so a measure "
            "with no domain has nothing to infer from. Pass domain=, or write the "
            "native spelling, which names the domain."
        )
    if len(distinct) > 1:
        raise ValueError(
            f"The functions in this integrand are on {len(distinct)} different domains, "
            f"so there is no one domain to infer. An integrand spanning two domains is "
            f"what the native spelling is for: name the domain and read each field "
            f"there through its own parametrization."
        )
    (domain,) = distinct
    if not isinstance(domain, AbstractParametrizedDomain):
        raise ValueError(f"{domain!r} is not made of cells, so it has no facets.")
    return domain


def _is_on(node: GraphNode, graph, domain: AbstractParametrizedDomain) -> bool:
    """Whether a subtree is entirely a field on the given domain.

    A subtree qualifies when it holds a function on that domain and no
    geometry of anywhere else, so that a facet normal stops the wrapping
    and a gradient does not.
    """
    nodes = graph.descendants(node)
    holds_field = any(
        isinstance(n, AbstractFunction) and n.function_space.domain == domain for n in nodes
    )
    holds_other_geometry = any(
        isinstance(n, AbstractGeometricQuantity) and n.domain != domain for n in nodes
    )
    return holds_field and not holds_other_geometry


def _restrict_fields(
    integrand: AbstractExpression,
    domain: AbstractParametrizedDomain,
    onto: AbstractParametrizedDomain,
) -> AbstractExpression:
    """Read every field on the domain at a point of the domain it is restricted to.

    The largest subtrees, walked from the root: a gradient is wrapped
    whole, so it keeps the cell's inverse Jacobian. Wrapping the function
    inside it instead would give the tangential gradient and lose the
    derivative across the facet.

    Args:
        integrand: The expression being integrated
        domain: The domain its fields are on
        onto: The domain they are read at a point of

    Returns:
        The integrand with its fields restricted
    """
    graph = as_graph(integrand)
    restrictions: dict[GraphNode, GraphNode] = {}
    covered: set[GraphNode] = set()
    for node in graph.ordered_nodes(NodeOrder.roots_first):
        if node in covered:
            continue
        if _is_on(node, graph, domain):
            assert isinstance(node, AbstractExpression)
            restrictions[node] = restrict(node, onto)
            covered |= graph.descendants(node)
    if len(restrictions) == 0:
        return integrand
    restricted = replace(integrand, restrictions)
    assert isinstance(restricted, AbstractExpression)
    return restricted


class CellMeasure:
    """`dx` as UFL writes it: a cell measure that reads its domain off the integrand.

    Call it to mark a subdomain, as `dx(1)`, which labels the domain, or
    to name the domain, as `dx(domain=mesh)`. `uflx.dx` is the native
    measure and takes the domain itself; importing this one shadows it,
    which is the point of importing it.
    """

    def __init__(self, subdomain_id: Any = None, domain: AbstractParametrizedDomain | None = None):
        """Initialise.

        Args:
            subdomain_id: What the consumer marked the cells with
            domain: The domain to integrate over, if it is not to be inferred
        """
        self._subdomain_id = subdomain_id
        self._domain = domain

    def __call__(
        self, subdomain_id: Any = None, domain: AbstractParametrizedDomain | None = None
    ) -> CellMeasure:
        """Get this measure over a marked subdomain, or over a named domain.

        Raises:
            ValueError: If the subdomain id is a domain, which is the
                native spelling of `dx`
        """
        if isinstance(subdomain_id, AbstractParametrizedDomain):
            raise ValueError(
                "This is the legacy dx, whose argument is a subdomain id. Write "
                "dx(domain=...) for a domain, or import dx from uflx for the native "
                "measure of a domain."
            )
        return CellMeasure(subdomain_id, domain or self._domain)

    @property
    def subdomain_id(self) -> Any:
        """What the consumer marked the cells with."""
        return self._subdomain_id

    def domain_of(self, integrand: AbstractExpression) -> AbstractParametrizedDomain:
        """The domain an integrand of this measure is over."""
        return self._domain or _inferred_domain(integrand)

    def __rmul__(self, integrand: AbstractExpression) -> Integral:
        """Right multiply by an expression to form an integral."""
        if not isinstance(integrand, AbstractExpression):
            return NotImplemented
        domain: AbstractCellularDomain = self.domain_of(integrand)
        if self._subdomain_id is not None:
            domain = labelled_domain(domain, Cells(self._subdomain_id))
        integral = integrand * measure_of(domain)
        assert isinstance(integral, Integral)
        return integral

    def __repr__(self) -> str:
        """Representation."""
        return f"CellMeasure({self._subdomain_id!r}, {self._domain!r})"


class ExteriorFacetMeasure(CellMeasure):
    """`ds` as UFL writes it: an exterior facet measure over the facets it infers.

    The integral it builds is the one measure over `facets_of(domain)`,
    with every field on the domain read there through its parametrization. Call it
    to mark a subdomain, as `ds(3)`, which labels that domain.
    """

    def __call__(
        self, subdomain_id: Any = None, domain: AbstractParametrizedDomain | None = None
    ) -> ExteriorFacetMeasure:
        """Get this measure over a marked part of the facets, or over a named domain."""
        if isinstance(subdomain_id, AbstractParametrizedDomain):
            raise ValueError(
                "This is the legacy ds, whose argument is a subdomain id. Write "
                "ds(domain=...) for a domain."
            )
        return ExteriorFacetMeasure(subdomain_id, domain or self._domain)

    def facets_for(self, integrand: AbstractExpression) -> UnionDomain:
        """The domain of facets an integrand of this measure is over."""
        return facets_of(self.domain_of(integrand), ExteriorFacets(self._subdomain_id))

    def __rmul__(self, integrand: AbstractExpression) -> Integral:
        """Right multiply by an expression to form an integral over the facets."""
        if not isinstance(integrand, AbstractExpression):
            return NotImplemented
        facets = self.facets_for(integrand)
        domain = self.domain_of(integrand)
        integral = _restrict_fields(integrand, domain, facets) * measure_of(facets)
        assert isinstance(integral, Integral)
        return integral

    def __repr__(self) -> str:
        """Representation."""
        return f"ExteriorFacetMeasure({self._subdomain_id!r}, {self._domain!r})"


class InteriorFacetMeasure(CellMeasure):
    """`dS` as UFL writes it, which cannot be built yet."""

    def __call__(
        self, subdomain_id: Any = None, domain: AbstractParametrizedDomain | None = None
    ) -> InteriorFacetMeasure:
        """Get this measure over a marked part of the interface."""
        return InteriorFacetMeasure(subdomain_id, domain or self._domain)

    def __rmul__(self, integrand: AbstractExpression) -> Integral:
        """Refuse, saying what an interface integral needs.

        Raises:
            NotImplementedError: Always
        """
        raise NotImplementedError(
            "An interior facet integral is this measure over an interface, which is a "
            "facet seen from both sides. That needs a parametrization per side, which a union of "
            "pieces does not express: its pieces carry no relation to each other, where "
            "the two sides of one interface are two parametrizations of one set related by a "
            "transition map. INTEGRATION.md section 5 says what it costs."
        )

    def __repr__(self) -> str:
        """Representation."""
        return f"InteriorFacetMeasure({self._subdomain_id!r}, {self._domain!r})"


dx = CellMeasure()
ds = ExteriorFacetMeasure()
dS = InteriorFacetMeasure()


def facet_normal(domain: AbstractParametrizedDomain) -> AbstractGeometricQuantity:
    """`FacetNormal` as UFL writes it, which cannot be built yet.

    UFL's facet normal is the *outward* one. UFLx gives a facet a unit
    normal across it with no sign, since nothing yet says which side of
    the facet its cell is on, so there is no honest drop-in.

    What it costs is small and INTEGRATION.md section 4 says what:
    outwardness is about the interior and not about orientation, so a
    reference facet's outward normal plus transport by `J^-T` gives the
    outward one even under a reflection. Basix holds the reference
    normals already.

    Args:
        domain: The domain whose facets are integrated over

    Raises:
        NotImplementedError: Always
    """
    raise NotImplementedError(
        "A facet normal with a sign needs an orientation: which side of the facet the "
        "cell is on. UnitNormal(facets_of(domain)) is the normal without one, which is "
        "a direction and not an outward direction."
    )
