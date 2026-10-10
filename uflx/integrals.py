# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
# Copyright (C) 2026 Jack S. Hale
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Measures and integrals."""

from __future__ import annotations

from abc import ABC, abstractmethod
from itertools import count
from typing import Any, cast

from uflx.algorithms import replace
from uflx.domains import AbstractParametrizedDomain, EntityDomain
from uflx.entities import AbstractEntity
from uflx.expressions import AbstractExpression
from uflx.functions import (
    AbstractFunction,
    AbstractVariable,
    FiniteElementVariable,
    create_variable,
)
from uflx.geometry import AbstractGeometricQuantity, VolumeElement
from uflx.graphs import Graph, GraphNode, as_graph, generate_graph


class AbstractMeasure(ABC):
    """Abstract base class for an integral measure.

    A measure is the domain integrated over paired with the density used
    on it. Neither half is optional: a function is not a differential
    form, so integrating one needs a density, and a density is a density
    of something.

    The domain is one presented as the image of a map, since that map is
    where the density comes from.
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

    def __add__(self, other: Any) -> IntegralSum:
        """Add to another integral, or to a sum of them."""
        if isinstance(other, IntegralSum):
            return IntegralSum((self, *other.terms))
        if isinstance(other, AbstractIntegral):
            return IntegralSum((self, other))
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
                domain being integrated over
        """
        self._measure = measure
        domain = measure.domain
        self._variable = create_variable(domain) if variable is None else variable

        replacements: dict[GraphNode, GraphNode] = {}
        for node in as_graph(integrand):
            if isinstance(node, AbstractFunction) and node.variable is None:
                if node.function_space.domain != domain:
                    raise ValueError(
                        f"Cannot integrate a function on {node.function_space.domain!r} "
                        f"over {domain!r}. A function reaches another domain only "
                        f"through a map, not by being integrated there."
                    )
                replacements[node] = node.reconstruct_with_variable(self._variable)

        # A Jacobian built during a pull back does not know where it is
        # evaluated. This integral's variable is that point.
        for node in as_graph(integrand):
            if isinstance(node, AbstractGeometricQuantity) and node.point is None:
                evaluated = node.reconstruct_with_variable(self._variable)
                if evaluated is not node:
                    replacements[node] = evaluated

        if len(replacements) == 0:
            self._integrand = integrand
        else:
            self._integrand = replace(integrand, replacements)
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
    def domain(self) -> AbstractParametrizedDomain:
        """The domain this integral is over, which its measure names."""
        return self._measure.domain

    def restricted_to(self, cell: AbstractEntity) -> Integral:
        """Get this integral over the part of its domain made of one cell type.

        Args:
            cell: A cell type of this integral's domain

        Returns:
            The same integrand over that cell type alone
        """
        restricted_domain = self.domain.restricted_to(cell)
        restrictions: dict[GraphNode, GraphNode] = {}
        for node in as_graph(self._integrand):
            if isinstance(node, AbstractFunction):
                restricted = node.restricted_to(cell)
                if restricted is not node:
                    restrictions[node] = restricted
        integrand = replace(self._integrand, restrictions) if restrictions else self._integrand
        assert isinstance(integrand, AbstractExpression)

        variable = self._variable
        if isinstance(variable, FiniteElementVariable):
            variable = variable.to_ambient_coordinates(restricted_domain)
        return Integral(integrand, self._measure.with_domain(restricted_domain), variable)

    def split_by_cell_type(self) -> IntegralSum | None:
        """Split this integral into one over each cell type of its domain.

        Returns:
            A sum of integrals, one per cell type, or None when the domain
            has one cell type and there is nothing to split
        """
        cell_types = self.domain.cell_types
        if len(cell_types) == 1:
            return None
        return IntegralSum(tuple(self.restricted_to(cell) for cell in cell_types))

    def pull_back_to_entity(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the entity's coordinates."""
        integrand = node_map.get(self._integrand, self._integrand)
        domain = self.domain
        if len(domain.cell_types) != 1:
            # Each cell type has its own coordinate domain, so there is no one
            # set of coordinates to pull back to. Split first.
            raise ValueError(
                "Cannot pull an integral over several cell types back to one cell's "
                "coordinates. Split it by cell type first."
            )
        (cell,) = domain.cell_types

        assert isinstance(integrand, AbstractExpression)

        # Change of variables moves the measure as well as the integrand.
        return Integral(
            self._measure.density * integrand,
            self._measure.with_domain(EntityDomain(cell)),
            self._variable.to_entity_coordinates(cell),
        )

    def __repr__(self) -> str:
        """Representation."""
        return f"Integral(variable={self._variable!r})"


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

    The density is ``sqrt(det g)``, so that volumes measured in a cell's
    coordinates agree with the ambient ones. A measure that weighs its
    domain some other way is a class of its own rather than an argument
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
