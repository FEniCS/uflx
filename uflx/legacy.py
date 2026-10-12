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

    facets = sub_entities_of(mesh, 1)
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

`dS` infers the mesh the same way and integrates over the interface of
its facets with its facets. A side is named by `plus(u)` and `minus(u)`,
and `jump`, `avg` and `facet_normal` are built on those with UFL's
signatures. A marked side is a placeholder until the measure knows the
interface, which is the inference again: which interface `u("+")` means
is read off the measure it ends up under, and nothing in `uflx` does
that for a field.

UFL's own spelling, `u("+")`, needs an expression node to be callable,
which the core's are not. `add_plus_minus_restriction(AbstractExpression)` makes
them so, putting a `__call__` on the class named that builds the same
marker, and `remove_plus_minus_restriction(AbstractExpression)` takes it off
again. That reaches into a core class rather than adding names beside
it, which is why it is a call the user makes, naming the class, and not
a side effect of importing this module: what a class does should not
depend on which modules happen to have been imported first, and which
class is changed should be written down where it is changed. A narrower
class, `Coefficient` say, is a smaller change. `plus(u)` is the spelling
a type checker accepts.
"""

from __future__ import annotations

from typing import Any, NamedTuple

from uflx.algorithms import replace
from uflx.domains import (
    AbstractCellularDomain,
    AbstractParametrizedDomain,
    UnionDomain,
    labelled_domain,
)
from uflx.expressions import AbstractExpression, RealScalar
from uflx.functions import AbstractFunction
from uflx.geometry import AbstractGeometricQuantity, OutwardNormal
from uflx.graphs import GraphNode, NodeOrder, as_graph
from uflx.integrals import Integral
from uflx.integrals import dx as measure_of
from uflx.interfaces import Interface
from uflx.operators import inner
from uflx.parametrizations import sub_entities_of
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
    """What a legacy interior facet measure means, as a domain's label."""

    subdomain_id: Any = None


class SideMarker(AbstractExpression):
    """A subtree marked as read on one side, before the interface is known.

    What `plus(u)` builds. It stands where a restriction will: the
    measure that ends up over it builds the interface and replaces it
    with the restriction to that interface's side, so until then it
    carries the values of what it marks and nothing else.
    """

    def __init__(self, expression: AbstractExpression, side: str):
        """Initialise.

        Args:
            expression: What is marked
            side: Which side, ``"+"`` or ``"-"``

        Raises:
            ValueError: If the side is neither, or the expression is
                marked already
        """
        if side not in ("+", "-"):
            raise ValueError(f"A side is '+' or '-', not {side!r}.")
        if any(isinstance(node, SideMarker) for node in as_graph(expression)):
            raise ValueError(f"{expression!r} is marked with a side already.")
        self._expression = expression
        self._side = side

    @property
    def expression(self) -> AbstractExpression:
        """What is marked."""
        return self._expression

    @property
    def side(self) -> str:
        """Which side."""
        return self._side

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The marked expression's shape."""
        return self._expression.value_shape

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._expression}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._expression, self._side

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return SideMarker(self._expression.component(*indices), self._side)

    def __repr__(self) -> str:
        """Representation."""
        return f"SideMarker({self._expression!r}, {self._side!r})"


def _restrict_to_side(self: AbstractExpression, side: str) -> SideMarker:
    """`u("+")` as UFL writes it, once `add_plus_minus_restriction` has made a call mean it."""
    if not isinstance(side, str):
        raise ValueError(
            f"An expression is called with a side, '+' or '-', not with {side!r}. "
            f"Evaluating one at a point is not what calling it means here."
        )
    return SideMarker(self, side)


def _own_call(cls: type) -> Any:
    """What a class itself says calling it means, if anything."""
    return vars(cls).get("__call__")


def add_plus_minus_restriction(cls: type[AbstractExpression]) -> None:
    """Make `u("+")` mean `plus(u)` for the given kind of expression node.

    Puts a `__call__` on the class named, so the change is written down
    where it is made rather than being a side effect of importing this
    module, and is as wide as the caller chooses: `AbstractExpression`
    for everything, `Coefficient` for coefficients alone. Doing it twice
    is fine.

    Args:
        cls: The class whose instances are to be callable with a side

    Raises:
        ValueError: If the class is not a kind of expression, or is
            callable already for some other reason, which this would
            silently override
    """
    if not isinstance(cls, type) or not issubclass(cls, AbstractExpression):
        raise ValueError(f"{cls!r} is not a kind of expression, so a side cannot be read of it.")
    if _own_call(cls) is _restrict_to_side:
        return
    if any("__call__" in vars(base) for base in cls.__mro__ if base is not object):
        raise ValueError(f"{cls.__name__} is callable already, so a side cannot claim it.")
    cls.__call__ = _restrict_to_side


def remove_plus_minus_restriction(cls: type[AbstractExpression]) -> None:
    """Take the `__call__` that `add_plus_minus_restriction` put on the given class off again.

    Nothing to do where it was never put on, and a `__call__` from
    anywhere else is left alone.

    Args:
        cls: The class `add_plus_minus_restriction` was given
    """
    if _own_call(cls) is _restrict_to_side:
        del cls.__call__


def plus(expression: AbstractExpression) -> SideMarker:
    """`u("+")` as UFL writes it: mark an expression as read on the ``+`` side."""
    return SideMarker(expression, "+")


def minus(expression: AbstractExpression) -> SideMarker:
    """`u("-")` as UFL writes it: mark an expression as read on the ``-`` side."""
    return SideMarker(expression, "-")


def jump(v: AbstractExpression, n: AbstractExpression | None = None) -> AbstractExpression:
    """`jump` as UFL writes it.

    ``v("+") - v("-")``, or with a normal ``v("+") n("+") + v("-") n("-")``,
    the product being an inner product where ``v`` is not a scalar.
    """
    if n is None:
        return plus(v) - minus(v)
    if v.value_shape == ():
        return plus(v) * plus(n) + minus(v) * minus(n)
    return inner(plus(v), plus(n)) + inner(minus(v), minus(n))


def avg(v: AbstractExpression) -> AbstractExpression:
    """`avg` as UFL writes it: ``(v("+") + v("-")) / 2``."""
    return RealScalar(0.5) * (plus(v) + minus(v))


def _inferred_domain(integrand: AbstractExpression) -> AbstractCellularDomain:
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
        if isinstance(node, AbstractFunction) and not node.is_on_an_entity
    ]
    distinct = list(dict.fromkeys(domains))
    if len(distinct) == 0:
        raise ValueError("No function in this integrand says which domain. Pass domain=.")
    if len(distinct) > 1:
        raise ValueError(
            f"This integrand's functions are on {len(distinct)} domains, so no one domain "
            f"can be inferred. Name it, as the native spelling does."
        )
    (domain,) = distinct
    if not isinstance(domain, AbstractCellularDomain):
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

    def __init__(self, subdomain_id: Any = None, domain: AbstractCellularDomain | None = None):
        """Initialise.

        Args:
            subdomain_id: What the consumer marked the cells with
            domain: The domain to integrate over, if it is not to be inferred
        """
        self._subdomain_id = subdomain_id
        self._domain = domain

    def __call__(
        self, subdomain_id: Any = None, domain: AbstractCellularDomain | None = None
    ) -> CellMeasure:
        """Get this measure over a marked subdomain, or over a named domain.

        Raises:
            ValueError: If the subdomain id is a domain, which is the
                native spelling of `dx`
        """
        if isinstance(subdomain_id, AbstractParametrizedDomain):
            raise ValueError(
                "The legacy dx takes a subdomain id. Write dx(domain=...), or use uflx.dx."
            )
        return CellMeasure(subdomain_id, domain or self._domain)

    @property
    def subdomain_id(self) -> Any:
        """What the consumer marked the cells with."""
        return self._subdomain_id

    def domain_of(self, integrand: AbstractExpression) -> AbstractCellularDomain:
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

    The integral it builds is the one measure over `sub_entities_of(domain, 1)`,
    with every field on the domain read there through its parametrization. Call it
    to mark a subdomain, as `ds(3)`, which labels that domain.
    """

    def __call__(
        self, subdomain_id: Any = None, domain: AbstractCellularDomain | None = None
    ) -> ExteriorFacetMeasure:
        """Get this measure over a marked part of the facets, or over a named domain."""
        if isinstance(subdomain_id, AbstractParametrizedDomain):
            raise ValueError("The legacy ds takes a subdomain id. Write ds(domain=...).")
        return ExteriorFacetMeasure(subdomain_id, domain or self._domain)

    def facets_for(self, integrand: AbstractExpression) -> UnionDomain:
        """The domain of facets an integrand of this measure is over."""
        return sub_entities_of(self.domain_of(integrand), 1, ExteriorFacets(self._subdomain_id))

    def __rmul__(self, integrand: AbstractExpression) -> Integral:
        """Right multiply by an expression to form an integral over the facets."""
        if not isinstance(integrand, AbstractExpression):
            return NotImplemented
        facets = self.facets_for(integrand)
        domain = self.domain_of(integrand)
        # The normal is read off the facets, so it is retargeted first and
        # then stops the fields around it being restricted as one subtree.
        integrand = _retarget_normals(integrand, domain, facets)
        integral = _restrict_fields(integrand, domain, facets) * measure_of(facets)
        assert isinstance(integral, Integral)
        return integral

    def __repr__(self) -> str:
        """Representation."""
        return f"ExteriorFacetMeasure({self._subdomain_id!r}, {self._domain!r})"


class InteriorFacetMeasure(CellMeasure):
    """`dS` as UFL writes it: an interior facet measure over the interface it infers.

    The integral it builds is the one measure over the interface the
    mesh's facets make with themselves, with every side marked by `plus`
    or `minus` read on the interface through that side. Call it to mark a
    subdomain, as `dS(3)`, which labels the facets of both sides.

    A field on the mesh that is not marked is refused, as UFL refuses an
    unrestricted field under `dS`, since a field reaches the interface
    only through a side.
    """

    def __call__(
        self, subdomain_id: Any = None, domain: AbstractCellularDomain | None = None
    ) -> InteriorFacetMeasure:
        """Get this measure over a marked part of the interface, or over a named domain."""
        if isinstance(subdomain_id, AbstractParametrizedDomain):
            raise ValueError("The legacy dS takes a subdomain id. Write dS(domain=...).")
        return InteriorFacetMeasure(subdomain_id, domain or self._domain)

    def interface_for(self, integrand: AbstractExpression) -> Interface:
        """The interface an integrand of this measure is over."""
        facets = sub_entities_of(self.domain_of(integrand), 1, InteriorFacets(self._subdomain_id))
        return Interface(facets, (("+", facets), ("-", facets)))

    def __rmul__(self, integrand: AbstractExpression) -> Integral:
        """Right multiply by an expression to form an integral over the interface."""
        if not isinstance(integrand, AbstractExpression):
            return NotImplemented
        interface = self.interface_for(integrand)
        domain = self.domain_of(integrand)
        integral = _resolve_sides(integrand, domain, interface) * measure_of(interface)
        assert isinstance(integral, Integral)
        return integral

    def __repr__(self) -> str:
        """Representation."""
        return f"InteriorFacetMeasure({self._subdomain_id!r}, {self._domain!r})"


def _retarget_normals(
    integrand: AbstractExpression,
    domain: AbstractCellularDomain,
    onto: AbstractParametrizedDomain,
) -> AbstractExpression:
    """Read the placeholder normal of a mesh off the domain it is a normal to."""
    retargeted: dict[GraphNode, GraphNode] = {
        node: node.on_domain(onto)
        for node in as_graph(integrand)
        if isinstance(node, OutwardNormal) and node.domain == domain
    }
    if len(retargeted) == 0:
        return integrand
    replaced = replace(integrand, retargeted)
    assert isinstance(replaced, AbstractExpression)
    return replaced


def _resolve_sides(
    integrand: AbstractExpression, domain: AbstractCellularDomain, interface: Interface
) -> AbstractExpression:
    """Replace each marked subtree by what it marks, read through the interface's side.

    The normal of the mesh becomes the normal out of that side, and
    every field in the subtree is restricted to the side, as a facet
    measure restricts them to the facets.
    """
    resolved: dict[GraphNode, GraphNode] = {}
    for node in as_graph(integrand):
        if not isinstance(node, SideMarker):
            continue
        side = interface.embedding(node.side)
        marked = _retarget_normals(node.expression, domain, side)
        resolved[node] = _restrict_fields(marked, domain, side)
    if len(resolved) == 0:
        return integrand
    replaced = replace(integrand, resolved)
    assert isinstance(replaced, AbstractExpression)
    return replaced


dx = CellMeasure()
ds = ExteriorFacetMeasure()
dS = InteriorFacetMeasure()


def facet_normal(domain: AbstractCellularDomain) -> AbstractGeometricQuantity:
    """`FacetNormal` as UFL writes it: the outward normal, read off the mesh.

    A normal is a quantity of the domain it is a normal to, which is a
    domain of facets and not the mesh, so this is a placeholder: the
    outward normal of the mesh's cells, which `ds` reads off the facets
    and `dS` reads off the side it is marked with, as `plus(n)`. The
    native spelling is `OutwardNormal(sub_entities_of(mesh, 1))`, or
    `OutwardNormal(interface.embedding("+"))`.

    Args:
        domain: The mesh whose cells the normal points out of

    Returns:
        The outward normal, to be retargeted by the measure it ends up
        under
    """
    return OutwardNormal(domain)
