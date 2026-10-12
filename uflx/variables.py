# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
# Copyright (C) 2026 Jack S. Hale
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Variables.

A variable stands for a point of a domain without naming one. An
integral's dummy variable is the one every expression in its integrand is
read at, and pulling the integral back onto a cell's coordinates gives a
variable there instead.

This sits beside `uflx.domains` rather than in `uflx.functions`: a
variable is a point of a domain, so a domain's parametrizations take one, and the
geometry built on them does too. What a function does with a variable is
`uflx.functions`'s business. A domain and a point of it are peers and
refer to each other, which is why `uflx.domains` names this module only
in its type hints.
"""

from __future__ import annotations

from abc import abstractmethod
from itertools import count
from typing import Any, Protocol, runtime_checkable

from uflx.domains import (
    AbstractCoordinateDomain,
    AbstractDomain,
    AbstractParametrizedDomain,
    EntityDomain,
    parameter_regions,
)
from uflx.entities import AbstractEntity
from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode


@runtime_checkable
class ReadAtAVariable(GraphNode, Protocol):
    """An expression that is told where it is read after being built.

    A function, a geometric quantity and a restriction are all built
    before there is a point to read them at -- during a pull back there
    is none yet -- and are told later, by whatever knows the point: an
    integral hands over its dummy variable, and a restriction hands over
    the point its parametrization gives. `variable` is None until then, and
    reconstructing gives the node back unchanged when the variable is in
    coordinates that are none of its business.
    """

    @property
    def variable(self) -> AbstractVariable | None:
        """Where this is read, if it has been told yet."""

    def reconstruct_with_variable(self, variable: AbstractVariable) -> GraphNode:
        """Read this at the given variable."""


class AbstractVariable(AbstractExpression):
    """Base class for a variable that is the input to a function."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self.domain.geometric_dimension,)

    @property
    @abstractmethod
    def domain(self) -> AbstractDomain:
        """The domain that this variable is in."""

    @abstractmethod
    def __eq__(self, other) -> bool:
        """Check for equality."""

    @abstractmethod
    def __hash__(self) -> int:
        """Hash."""

    @property
    def is_on_an_entity(self) -> bool:
        """Whether this variable is a point of an entity's coordinate domain."""
        return isinstance(self.domain, EntityDomain)

    def to_entity_coordinates(self, cell: AbstractEntity) -> AbstractVariable:
        """Make a version of this variable in the given cell's coordinates.

        One already there is itself, a named point of a cell being as good
        a place to read a field as a variable standing for one.
        """
        if self.domain == EntityDomain(cell):
            return self
        raise ValueError("Cannot pull this variable back to a cell's coordinates.")

    def component(self, *indices: int) -> AbstractExpression:
        """Get one coordinate of the point this variable stands for.

        A coordinate domain's points are tuples of numbers, so each of
        those numbers is a scalar expression. A domain whose points are
        not tuples has no coordinates to name, which is why this asks
        the variable's domain rather than only its dimension.

        Args:
            indices: Which coordinate

        Returns:
            That coordinate of this variable

        Raises:
            ValueError: If the points of this variable's domain are not
                coordinate tuples, or there is no such coordinate
        """
        (index,) = indices
        domain = self.domain
        if not isinstance(domain, AbstractCoordinateDomain):
            raise ValueError(
                f"A point of {domain!r} is not a coordinate tuple, so it has no "
                f"coordinates to name. Split the integral by piece, and pull it back "
                f"onto a cell, first."
            )
        if not 0 <= index < domain.geometric_dimension:
            raise ValueError(
                f"A point of {domain!r} is named by {domain.geometric_dimension} "
                f"coordinates, so it has no coordinate {index}."
            )
        return PointComponent(self, index)


class PointComponent(AbstractExpression):
    """One coordinate of the point a variable stands for.

    A terminal, and the one every parametrization written in closed form bottoms
    out in: such a parametrization does arithmetic on the coordinates it is handed,
    where a finite element parametrization sums a basis against them and never
    reads them. A point that is named gives its own entry instead, so
    this stands for a coordinate nobody has written down yet.
    """

    def __init__(self, point: AbstractVariable, component: int | str):
        """Initialise.

        Args:
            point: The variable whose coordinate this is
            component: Which coordinate
        """
        self._point = point
        self._component = component

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression.

        Raises:
            ValueError: Always, a coordinate being a scalar
        """
        raise ValueError("Cannot get a component of a scalar expression")

    @property
    def component_index(self) -> int | str:
        """Get the component of the point."""
        return self._component

    @property
    def point(self) -> AbstractVariable:
        """Get the point."""
        return self._point

    def __repr__(self):
        """Representation."""
        return f"PointComponent({self._point}, {self._component})"

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._point}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._point, self._component

    def get_replacement(self, replacements: dict[GraphNode, GraphNode]) -> GraphNode | None:
        """Get the node to replace this node with, or None if no replacement can be made.

        Substituting a point for a variable is what turns a parametrization's
        arithmetic into numbers, so a coordinate of the variable becomes
        that coordinate of the point rather than being rebuilt around it.
        """
        replaced = replacements.get(self._point)
        if isinstance(replaced, AbstractVariable) and isinstance(self._component, int):
            return replaced.component(self._component)
        return None


class Variable(AbstractVariable):
    """A variable that is the input to a function.

    One in a cellular domain can be moved between that domain's own
    coordinates and one of its cells', which is what pulling an integral
    back does to its dummy variable. One in any other domain has nowhere
    to be moved to and refuses.
    """

    _n = count(0)

    def __init__(self, domain: AbstractDomain, label: str | None = None):
        """Initialise."""
        if label is None:
            self._label = f"variable-{next(self._n)}"
        else:
            self._label = label
        self._domain = domain

    @property
    def label(self) -> str:
        """The label of this variable."""
        return self._label

    @property
    def domain(self) -> AbstractDomain:
        """The domain that this variable is in."""
        return self._domain

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._domain, self._label)

    @property
    def _regions(self) -> set[AbstractCoordinateDomain]:
        """The parameter regions this variable's domain is parametrized from.

        A cell's coordinate domain is one, which is where a variable goes
        when an integral is pulled back onto that cell. The domain need
        not be made of cells to have one, an interface over a facet
        being parametrized from the reference facet as the facet is.

        Raises:
            ValueError: If the domain is not parametrized at all
        """
        if not isinstance(self._domain, AbstractParametrizedDomain):
            raise ValueError(
                f"{self._domain!r} is not made of cells, so a point of it is in no "
                f"cell's coordinates."
            )
        return parameter_regions(self._domain)

    def to_domain_coordinates(self, domain: AbstractParametrizedDomain) -> Variable:
        """Make a version of this variable in a domain's own coordinates.

        Args:
            domain: The domain to push the variable onto. A variable in a
                cell's coordinates does not know which domain it came
                from, and those coordinates belong to no domain in
                particular, so the target is named rather than recalled.
        """
        return Variable(domain, self._label)

    def to_entity_coordinates(self, cell: AbstractEntity) -> Variable:
        """Make a version of this variable in the given cell's coordinates.

        Args:
            cell: A cell type of this variable's domain. A sub-entity of
                one is not a cell of the domain, and a point in its
                coordinates reaches a cell element only through an
                inclusion, not through this.

        Raises:
            ValueError: If the cell is not a cell type of this variable's
                domain
        """
        if EntityDomain(cell) not in self._regions:
            raise ValueError(
                f"Cannot pull back to {cell}, which is not a cell type of this variable's domain."
            )
        return Variable(EntityDomain(cell), self._label)

    def __eq__(self, other) -> bool:
        """Check for equality.

        The domain counts: pulling a variable back keeps its label and
        changes only where it lives.
        """
        return (
            isinstance(other, Variable)
            and other.label == self.label
            and other.domain == self.domain
        )

    def __repr__(self) -> str:
        """Representation."""
        return f"Variable({self._label}, {self._domain!r})"

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.Variable", self._label, self._domain))


def create_variable(domain: AbstractDomain) -> Variable:
    """Create a new variable in a domain."""
    return Variable(domain)
