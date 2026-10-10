# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
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
from uflx.domains import AbstractCellularDomain, AbstractParametrizedDomain
from uflx.entities import AbstractEntity
from uflx.expressions import AbstractExpression
from uflx.functions import (
    AbstractFunction,
    AbstractVariable,
    FiniteElementVariable,
    create_variable,
    extract_domain,
)
from uflx.geometry import AbstractGeometricQuantity, JacobianDeterminant
from uflx.graphs import Graph, GraphNode, as_graph, generate_graph


class AbstractMeasure(ABC):
    """Abstract base class for an integral measure."""

    def __rmul__(self, other: AbstractExpression) -> Integral:
        """Right multiply by an expression to form an integral."""
        if isinstance(other, AbstractExpression):
            return Integral(other, self)
        return NotImplemented

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    @abstractmethod
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""


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
        """Initialise."""
        self._measure = measure
        replacements: dict[GraphNode, GraphNode] = {}
        if variable is None:
            domain = None
            for node in as_graph(integrand):
                if isinstance(node, AbstractFunction) and node.variable is None:
                    if domain is None:
                        domain = node.function_space.domain
                        self._variable = create_variable(domain)
                    else:
                        assert domain == node.function_space.domain
                    replacements[node] = node.reconstruct_with_variable(self._variable)
            assert domain is not None
        else:
            self._variable = variable

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
    def domain(self) -> AbstractCellularDomain:
        """The domain this integral is over.

        Read off the functions in the integrand that are in ambient
        coordinates, which must all agree.
        """
        domain = extract_domain(self._integrand)
        assert isinstance(domain, AbstractCellularDomain)
        return domain

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
        return Integral(integrand, self._measure, variable)

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
        assert isinstance(domain, AbstractParametrizedDomain)
        if len(domain.cell_types) != 1:
            # Each cell type has its own coordinate domain, so there is no one
            # set of coordinates to pull back to. Split first.
            raise ValueError(
                "Cannot pull an integral over several cell types back to one cell's "
                "coordinates. Split it by cell type first."
            )
        (cell,) = domain.cell_types
        det = abs(JacobianDeterminant(domain))

        assert isinstance(integrand, AbstractExpression)

        return Integral(det * integrand, self._measure, self._variable.to_entity_coordinates(cell))

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
    """An integral measure."""

    def __init__(
        self, dim: int | None = None, codim: int | None = None, boundary_only: bool = False
    ):
        """Initialise."""
        self._dim = dim
        self._codim = codim
        self._boundary_only = boundary_only

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._dim, self._codim, self._boundary_only

    def __repr__(self) -> str:
        """Representation."""
        kwargs = {}
        if self._dim is not None:
            kwargs["dim"] = self._dim
        if self._codim is not None:
            kwargs["codim"] = self._codim
        if self._boundary_only:
            kwargs["boundary_only"] = self._boundary_only
        return (
            f"{self.__class__.__name__}("
            + ", ".join(f"{key}={value}" for key, value in kwargs.items())
            + ")"
        )


dx = Measure(codim=0)
