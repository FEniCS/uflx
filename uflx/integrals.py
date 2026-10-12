# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
# Copyright (C) 2026 Jack S. Hale
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Measures and integrals."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable
from itertools import count
from typing import Any, cast

from uflx.algorithms import replace
from uflx.domains import (
    AbstractCellularDomain,
    AbstractDomain,
    AbstractParametrizedDomain,
    EntityDomain,
    parameter_regions,
    pieces_of,
    unlabelled,
)
from uflx.entities import AbstractEntity
from uflx.expressions import AbstractExpression
from uflx.functions import AbstractFunction
from uflx.geometry import AbstractGeometricQuantity, VolumeElement
from uflx.graphs import Graph, GraphNode, as_graph, generate_graph
from uflx.interfaces import Embedding, Interface
from uflx.restrictions import Restriction, evaluate_at
from uflx.variables import (
    AbstractVariable,
    Variable,
    create_variable,
)


class AbstractMeasure(ABC):
    """Abstract base class for an integral measure.

    A measure is the domain integrated over paired with the density used
    on it. Neither half is optional: a function is not a differential
    form, so integrating one needs a density, and a density is a density
    of something.

    The domain is a parametrized one, since the map is where the density
    comes from. It need not be made of cells: the Lebesgue measure of R^d
    is a measure of this kind, and so is one on a region or on its
    boundary.
    """

    def __rmul__(self, other: AbstractExpression) -> Integral:
        """Right multiply by an expression to form an integral."""
        if isinstance(other, AbstractExpression):
            return Integral(other, self)
        return NotImplemented

    @property
    @abstractmethod
    def domain(self) -> AbstractParametrizedDomain:
        """The domain this measure integrates over."""

    @property
    @abstractmethod
    def density(self) -> AbstractExpression:
        """The factor an integral picks up on being pulled back.

        In coordinates on the domain, the density written against the
        coordinate volume.
        """

    @abstractmethod
    def with_domain(self, domain: AbstractParametrizedDomain) -> AbstractMeasure:
        """Get this measure over another domain.

        Retargeting is what a pull back does to a measure, and what
        restricting an integral to one cell type does to it.

        Args:
            domain: The domain to integrate over instead

        Returns:
            The same kind of measure over that domain
        """

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node.

        A measure has none, the density included. A measure describes an
        integral rather than being part of its expression, and its
        density is not evaluated anywhere until a pull back multiplies it
        into the integrand, which is where it is told its point.
        """
        return set()

    @property
    @abstractmethod
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""

    def __eq__(self, other) -> bool:
        """Check for equality.

        Two measures are equal when they are the same kind over the same
        arguments. A measure is a description, so it has no identity of
        its own beyond them.
        """
        if not isinstance(other, AbstractMeasure):
            return NotImplemented
        return type(self) is type(other) and self.init_args == other.init_args

    def __hash__(self) -> int:
        """Hash."""
        return hash((type(self).__name__, self.init_args))


class AbstractIntegral(ABC):
    """Abstract base class for an integral."""

    @property
    @abstractmethod
    def integrand(self) -> AbstractExpression:
        """The integrand."""

    @property
    @abstractmethod
    def measure(self) -> AbstractMeasure:
        """The integral measure."""

    @property
    def domain(self) -> AbstractParametrizedDomain:
        """The domain this integral is over, which its measure names."""
        return self.measure.domain

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.integrand, self.measure}

    @property
    @abstractmethod
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""

    def __eq__(self, other):
        """Check for equality."""
        if isinstance(other, AbstractIntegral):
            return self.integrand == other.integrand and self.measure == other.measure
        return NotImplemented

    def __hash__(self):
        """Hash."""
        return hash((hash(self.integrand), hash(self.measure)))

    @property
    @abstractmethod
    def variable(self) -> AbstractVariable:
        """The dummy variable of this integral."""

    def split_by_piece(self) -> IntegralSum | None:
        """Split this integral into one over each piece of its domain.

        None unless there are pieces to split by, which an integral over
        a union of them overrides to say.
        """
        return None

    def split_by_cell_type(self) -> IntegralSum | None:
        """Split this integral into one over each cell type of its domain.

        None unless there is more than one cell type, which an integral
        over a domain of several overrides to say.
        """
        return None

    def __add__(self, other: Any) -> IntegralSum:
        """Add to another integral, or to a sum of them."""
        if isinstance(other, IntegralSum):
            return IntegralSum((self, *other.terms))
        if isinstance(other, AbstractIntegral):
            return IntegralSum((self, other))
        return NotImplemented

    @abstractmethod
    def __neg__(self) -> AbstractIntegral:
        """Negate, which negates the integrand."""

    def __sub__(self, other: Any) -> IntegralSum:
        """Subtract another integral, or a sum of them, by adding its negation."""
        if isinstance(other, AbstractIntegral | IntegralSum):
            return self + (-other)
        return NotImplemented


class Integral(AbstractIntegral):
    """An integral."""

    _n = count(0)

    def __init__(
        self,
        integrand: AbstractExpression,
        measure: AbstractMeasure,
        variable: AbstractVariable | None = None,
    ):
        """Initialise.

        Args:
            integrand: The expression being integrated
            measure: The measure to integrate it against, which says over
                what
            variable: The point the integrand is read at, if one is being
                carried over from the integral this one was built from

        Raises:
            ValueError: If a function in the integrand is not on the
                domain being integrated over, or a geometric quantity in
                it cannot be read there
        """
        self._measure = measure
        domain = measure.domain
        self._variable = create_variable(domain) if variable is None else variable

        # A restriction says where the field under it is read, so this
        # integral's variable is none of that subtree's business, and a
        # field on another domain is no error there either: reaching
        # another domain through a map is exactly what a restriction is.
        graph = as_graph(integrand)
        restricted = frozenset(
            descendant
            for node in graph
            if isinstance(node, Restriction)
            for descendant in graph.descendants(node) - {node}
        )

        for node in graph:
            if node in restricted:
                continue
            if isinstance(node, AbstractFunction) and node.variable is None:
                # Labels are taken off first: a label says which entities a
                # consumer means and nothing about the fields on them, so a
                # function on a domain is a function on any selection of it.
                if not _is_on_or_an_embedding_of(node.function_space.domain, domain):
                    raise ValueError(
                        f"Cannot integrate a function on {node.function_space.domain!r} "
                        f"over {domain!r}. A function reaches another domain only "
                        f"through a map, not by being integrated there."
                    )

        # A function or a Jacobian built during a pull back does not know
        # where it is evaluated. This integral's variable is that point.
        evaluated = evaluate_at(integrand, self._variable, restricted)
        assert isinstance(evaluated, AbstractExpression)

        # The same question of a geometric quantity, which until now was
        # asked of functions alone: a quantity of another domain in this
        # integrand built and pulled back and only failed on being
        # expanded. One that already carries its own point is a quantity a
        # pull back or a restriction put here, and says where it is read
        # itself.
        for node in as_graph(evaluated):
            if node in restricted or not isinstance(node, AbstractGeometricQuantity):
                continue
            if node.variable is None and not _is_on_or_an_embedding_of(node.domain, domain):
                raise ValueError(
                    f"Cannot read a quantity of {node.domain!r} over {domain!r}. A "
                    f"quantity of another domain reaches this one through a map, as a "
                    f"function does."
                )

        self._integrand = evaluated
        self._graph = generate_graph(self)

    @property
    def variable(self) -> AbstractVariable:
        """The dummy variable of this integral."""
        return self._variable

    @property
    def integrand(self) -> AbstractExpression:
        """The integrand."""
        return cast(AbstractExpression, self._integrand)

    @property
    def measure(self) -> AbstractMeasure:
        """The integral measure."""
        return self._measure

    @property
    def graph(self) -> Graph:
        """The graph that represents this object."""
        return self._graph

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._integrand, self._measure, self._variable

    @property
    def cellular_domain(self) -> AbstractCellularDomain:
        """The domain this integral is over, which must be made of cells.

        What the cell type fan out and the pull back both need, neither
        of them having anything to say about a domain parametrized some other
        way. Both want the cells and the map out of each, which is what a
        parametrized domain is.

        Raises:
            ValueError: If this integral's domain is not made of cells
        """
        domain = self.domain
        if not isinstance(domain, AbstractCellularDomain):
            raise ValueError(
                f"{domain!r} is not made of cells, so it has no cell types. Only a "
                f"domain of cells with a map out of each can be restricted to one of "
                f"them or pulled back onto one."
            )
        return domain

    def on_cell_type(self, cell: AbstractEntity) -> Integral:
        """Get this integral over the part of its domain made of one cell type.

        Args:
            cell: A cell type of this integral's domain

        Returns:
            The same integrand over that cell type alone
        """
        domain = self.cellular_domain.on_cell_type(cell)
        narrowed: dict[GraphNode, GraphNode] = {}
        for node in as_graph(self._integrand):
            if isinstance(node, AbstractFunction):
                on_cell = node.on_cell_type(cell)
                if on_cell is not node:
                    narrowed[node] = on_cell
        integrand = replace(self._integrand, narrowed) if narrowed else self._integrand
        assert isinstance(integrand, AbstractExpression)

        variable = self._variable
        if isinstance(variable, Variable):
            variable = variable.to_domain_coordinates(domain)
        return Integral(integrand, self._measure.with_domain(domain), variable)

    def on_piece(self, piece: AbstractParametrizedDomain) -> Integral:
        """Get this integral over one piece of its domain.

        The measure moves onto the piece, and so does everything in the
        integrand that was written against the whole union: a restriction
        is read through that piece's parametrization, and a geometric
        quantity is read off that piece's geometry.

        Args:
            piece: One piece of this integral's domain

        Returns:
            The same integrand over that piece alone
        """
        domain = self.domain
        replacements: dict[GraphNode, GraphNode] = {}
        for node in as_graph(self._integrand):
            if isinstance(node, Restriction):
                onto = _onto_piece(node.domain, domain, piece)
                if onto is not None:
                    replacements[node] = node.on_piece(onto)
            elif isinstance(node, AbstractGeometricQuantity):
                onto = _onto_piece(node.domain, domain, piece)
                if onto is not None:
                    replacements[node] = node.on_domain(onto)
        integrand = replace(self._integrand, replacements) if replacements else self._integrand
        assert isinstance(integrand, AbstractExpression)
        # The variable keeps its label, as it does on restricting to a cell
        # type, so that a function living on the domain itself and told the
        # union's variable is the same point after the pull back.
        variable = self._variable
        if isinstance(variable, Variable):
            variable = variable.to_domain_coordinates(piece)
        return Integral(integrand, self._measure.with_domain(piece), variable)

    def split_by_piece(self) -> IntegralSum | None:
        """Split this integral into one over each piece of its domain.

        Returns:
            A sum of integrals, one per piece, or None when the domain is
            not made of pieces

        Note:
            A union of one piece is still split. The union is a different
            domain from the piece, and it is the piece that names the
            parametrization a restriction reads a field through. A domain
            that is its own one piece, as an interface of two pieces is,
            has nothing to split.
        """
        domain = self.domain
        pieces = pieces_of(domain)
        if pieces == (domain,):
            return None
        return IntegralSum(tuple(self.on_piece(p) for p in pieces))

    def split_by_cell_type(self) -> IntegralSum | None:
        """Split this integral into one over each cell type of its domain.

        Returns:
            A sum of integrals, one per cell type, or None when the domain
            has one cell type and there is nothing to split
        """
        domain = self.domain
        if not isinstance(domain, AbstractCellularDomain):
            # Nothing to fan out over: cell types are what this splits by.
            return None
        cell_types = domain.cell_types
        if len(cell_types) == 1:
            return None
        return IntegralSum(tuple(self.on_cell_type(cell) for cell in cell_types))

    @property
    def region_entity(self) -> AbstractEntity:
        """The cell whose coordinates this integral is pulled back onto.

        Change of variables is onto the parameter region of the domain's
        one parametrization, which is INTEGRATION.md section 1 as stated.
        A domain of several parametrizations has no one region, and a
        region that is not a cell's coordinate domain has no cell to hand
        the integrand a variable of.

        Raises:
            ValueError: If this integral's domain has several
                parametrizations, or its one region is not a cell's
        """
        domain = self.domain
        sources = parameter_regions(domain)
        if len(sources) != 1:
            # Each parametrization has its own region, so there is no one set
            # of coordinates to pull back to. Split first.
            raise ValueError(
                f"Cannot pull an integral over {len(sources)} parameter regions back to "
                f"one of them. Split it by cell type, or by piece, first."
            )
        (source,) = sources
        if not isinstance(source, EntityDomain):
            raise ValueError(
                f"{domain!r} is not made of cells, so it is parametrized from "
                f"{source!r} and not from a cell's coordinates. Only an integral over a "
                f"domain of cells is pulled back onto one."
            )
        (cell,) = source.cell_types
        return cell

    def pull_back_to_entity(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the entity's coordinates."""
        integrand = node_map.get(self._integrand, self._integrand)
        cell = self.region_entity

        assert isinstance(integrand, AbstractExpression)

        # Change of variables moves the measure as well as the integrand.
        return Integral(
            self._measure.density * integrand,
            self._measure.with_domain(EntityDomain(cell)),
            self._variable.to_entity_coordinates(cell),
        )

    def __neg__(self) -> Integral:
        """Negate the integrand."""
        return Integral(-self._integrand, self._measure, self._variable)

    def __repr__(self) -> str:
        """Representation."""
        return f"Integral(variable={self._variable!r})"


def _bare(domain: AbstractDomain) -> AbstractDomain:
    """The domain with its labels taken off, and an interface as the domain it is over.

    Everything geometric about an interface is its domain's, the sides
    being structure on top, so a field living on that domain lives on
    the interface, which is where a multiplier on an interface lives.
    """
    domain = unlabelled(domain)
    if isinstance(domain, Interface):
        domain = unlabelled(domain.domain)
    return domain


def _is_on_or_an_embedding_of(where: AbstractDomain, domain: AbstractDomain) -> bool:
    """Whether a field or quantity of one domain is read in an integral over another.

    Labels are taken off first, a label saying which entities a consumer
    means and nothing about the geometry on them. A side of an interface
    is read in an integral over the interface, since the side is the
    interface as one domain presents it.
    """
    if _bare(where) == _bare(domain):
        return True
    return isinstance(where, Embedding) and _bare(where.interface) == _bare(domain)


def _onto_piece(
    where: AbstractDomain, domain: AbstractDomain, piece: AbstractParametrizedDomain
) -> AbstractParametrizedDomain | None:
    """Where something written against a whole domain moves when it is split.

    The piece itself for something written against the whole, and the
    piece's side of the same label for something written against a side
    of it. None for something written against anywhere else.
    """
    if where == domain:
        return piece
    if isinstance(where, Embedding) and where.interface == domain and isinstance(piece, Interface):
        return piece.embedding(where.label)
    return None


class IntegralSum:
    """A sum of integrals.

    An integral over a domain of several cell types is the sum of one
    integral per cell type, since each has its own coordinate domain, its
    own map out of it and so its own measure. This holds those terms.

    Not called a form: LANGUAGE.md defines a form as a multilinear
    functional and says a UFLx form need not be a sum of integrals.
    """

    def __init__(self, terms: tuple[AbstractIntegral, ...]):
        """Initialise.

        Args:
            terms: The integrals being added
        """
        if len(terms) == 0:
            raise ValueError("Cannot create an empty sum of integrals.")
        self._terms = terms

    @property
    def terms(self) -> tuple[AbstractIntegral, ...]:
        """The integrals being added."""
        return self._terms

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set(self._terms)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._terms,)

    def __add__(self, other: Any) -> IntegralSum:
        """Add another integral, or another sum of them, flattening the result."""
        if isinstance(other, IntegralSum):
            return IntegralSum((*self._terms, *other.terms))
        if isinstance(other, AbstractIntegral):
            return IntegralSum((*self._terms, other))
        return NotImplemented

    def __neg__(self) -> IntegralSum:
        """Negate every term."""
        return IntegralSum(tuple(-term for term in self._terms))

    def __sub__(self, other: Any) -> IntegralSum:
        """Subtract another integral, or another sum of them, by adding its negation."""
        if isinstance(other, AbstractIntegral | IntegralSum):
            return self + (-other)
        return NotImplemented

    def _split_terms(
        self, split: Callable[[AbstractIntegral], GraphNode | None]
    ) -> IntegralSum | None:
        """Replace each term by its own split, where it has one.

        A sum is where a fan out reaches the terms. Pulling a sum back
        walks it leaves first and so reaches each integral already, with
        no chance to split it on the way, and an integral over several
        pieces or several cell types cannot be pulled back whole.

        Args:
            split: What to ask of a term, giving None when it does not split

        Returns:
            The sum with each term split, or None when no term split
        """
        terms: list[AbstractIntegral] = []
        split_any = False
        for term in self._terms:
            result = split(term)
            if result is None:
                terms.append(term)
            elif isinstance(result, IntegralSum):
                terms.extend(result.terms)
                split_any = True
            else:
                assert isinstance(result, AbstractIntegral)
                terms.append(result)
                split_any = True
        if not split_any:
            return None
        return IntegralSum(tuple(terms))

    def split_by_piece(self) -> IntegralSum | None:
        """Split each term over a union of pieces into one integral per piece."""
        return self._split_terms(lambda term: term.split_by_piece())

    def split_by_cell_type(self) -> IntegralSum | None:
        """Split each term over several cell types into one integral per cell type."""
        return self._split_terms(lambda term: term.split_by_cell_type())

    def __eq__(self, other) -> bool:
        """Check for equality.

        The order of the terms counts, as it does for any other sum in
        UFLx before it is simplified.
        """
        return isinstance(other, IntegralSum) and self._terms == other._terms

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.IntegralSum", self._terms))

    def __repr__(self) -> str:
        """Representation."""
        return "IntegralSum(" + ", ".join(repr(i) for i in self._terms) + ")"


class Measure(AbstractMeasure):
    """Integration against the density a domain's own parametrization induces.

    The density is ``sqrt(det g)``, so that lengths and volumes measured
    in a parameter region's coordinates agree with the ambient ones. Where
    the map is the identity, as it is on R^d or on a region of it, that is
    one and the measure is the Lebesgue one. A measure that weighs its domain
    some other way is a class of its own rather than an argument
    defaulted here.
    """

    def __init__(self, domain: AbstractParametrizedDomain):
        """Initialise.

        Args:
            domain: The domain to integrate over
        """
        self._domain = domain

    @property
    def domain(self) -> AbstractParametrizedDomain:
        """The domain this measure integrates over."""
        return self._domain

    @property
    def density(self) -> AbstractExpression:
        """The volume element of this measure's domain."""
        return VolumeElement(self._domain)

    def with_domain(self, domain: AbstractParametrizedDomain) -> Measure:
        """Get the measure of another domain."""
        return Measure(domain)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._domain,)

    def __repr__(self) -> str:
        """Representation."""
        return f"{self.__class__.__name__}({self._domain!r})"


def dx(domain: AbstractParametrizedDomain) -> Measure:
    """Create the measure a domain's own parametrization induces.

    Args:
        domain: The domain to integrate over

    Returns:
        The measure of that domain
    """
    return Measure(domain)
