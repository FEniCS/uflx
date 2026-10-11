# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Functions.

A function is an item contained in a function space.
"""

from __future__ import annotations

from abc import abstractmethod
from itertools import count
from typing import Any, Self, cast

from uflx.domains import AbstractDomain, EntityDomain, ReadOnADomain
from uflx.entities import AbstractEntity
from uflx.expressions import AbstractExpression, Im, Re
from uflx.function_spaces import (
    AbstractFunctionSpace,
    AbstractMappedFunctionSpace,
    function_space,
)
from uflx.graphs import GraphNode
from uflx.maps import PushedForward
from uflx.tensors import zero
from uflx.variables import AbstractVariable


class AbstractFunction(AbstractExpression):
    """Base class for a function."""

    @property
    @abstractmethod
    def variable(self) -> AbstractVariable | None:
        """Get the variable that is this function's input."""

    @property
    def reads_on(self) -> AbstractDomain | None:
        """The domain this function is read on, until it is pulled back.

        A function on a cell's own coordinates has been pulled back and
        no longer says which domain it came from, which is what
        `is_on_an_entity` asks.
        """
        if self.is_on_an_entity:
            return None
        return self.function_space.domain

    @property
    def is_on_an_entity(self) -> bool:
        """Whether this function is on an entity's own coordinate domain.

        Which is to say whether it has been pulled back: a function on a
        mesh is not, and the same function on a cell is. Not the question
        of which coordinates its points are named in, which a domain's
        `coordinate_domain` answers, and which is a cell's for a facet of
        one without the function having been pulled back onto it.

        The variable answers when there is one, since a basis function
        keeps its mesh space and is told which coordinates by the point
        it is evaluated at. Otherwise the space's domain says.
        """
        if self.variable is not None:
            return self.variable.is_on_an_entity
        return isinstance(self.function_space.domain, EntityDomain)

    def on_cell_type(self, cell: AbstractEntity) -> AbstractFunction:
        """Get this function on the part of its space made of one cell type.

        Args:
            cell: A cell type of this function's space's domain

        Returns:
            The same function on that part of the space, or this function
            when its space does not change
        """
        space = self.function_space
        if not isinstance(space, AbstractMappedFunctionSpace):
            return self
        on_cell = space.on_cell_type(cell)
        if on_cell == space:
            return self
        return self.reconstruct_with_space(on_cell)

    @abstractmethod
    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""

    @abstractmethod
    def reconstruct_with_space(self, space: AbstractFunctionSpace) -> Self:
        """Reconstruct this function in another function space."""

    @abstractmethod
    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""

    @property
    @abstractmethod
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""

    @property
    @abstractmethod
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""

    @property
    @abstractmethod
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        if self.is_on_an_entity:
            assert isinstance(self.function_space, AbstractMappedFunctionSpace)
            return self.function_space.elements[0].entity_value_shape
        else:
            return self.function_space.value_shape

    @property
    def is_cellwise_constant(self) -> bool:
        """Whether this function's value is the same everywhere on every cell.

        True when every element of the function space is known to lie in
        the degree-0 Lagrange space (lagrange_superdegree == 0) and,
        for a function in its domain's own coordinates, every element's value
        map is known to preserve that constancy when pushed forward
        (see AbstractValueMap.preserves_constant_values).
        """
        assert isinstance(self.function_space, AbstractMappedFunctionSpace)
        elements = self.function_space.elements
        if any(e.lagrange_superdegree != 0 for e in elements):
            return False
        if self.is_on_an_entity:
            return True
        return all(e.value_map.preserves_constant_values for e in elements)

    @property
    def domain_size(self) -> int:
        """The size of the domain (ie the number of inputs to the function)."""
        assert self.function_space.domain.topological_dimension is not None
        return self.function_space.domain.topological_dimension

    @property
    def re(self) -> AbstractExpression:
        """Get real part."""
        if self.function_space.real_valued:
            return self
        else:
            return Re(self)

    @property
    def im(self) -> AbstractExpression:
        """Get imaginary part."""
        if self.function_space.real_valued:
            return zero(self.value_shape)
        else:
            return Im(self)


def _pulled_back_to_entity(function: AbstractFunction) -> GraphNode:
    """Read a function on the element of the cell its space is on.

    The same for every function built from a space, which is why it is
    written once here: the space moves onto the cell's own coordinate
    domain, the variable moves with it, and the element's value map
    carries the values. Which function it is, the function says, through
    `reconstruct_with_space`.

    Args:
        function: The function to pull back

    Returns:
        That function in the entity's coordinates, pushed forward

    Raises:
        ValueError: If the function is on an entity already
    """
    if function.is_on_an_entity:
        raise ValueError("This function is on an entity already.")
    space = function.function_space
    assert isinstance(space, AbstractMappedFunctionSpace)
    element = space.elements[0]
    pulled_back = function.reconstruct_with_space(
        function_space(EntityDomain(element.cell), element)
    )
    if function.variable is not None:
        pulled_back = pulled_back.reconstruct_with_variable(
            function.variable.to_entity_coordinates(element.cell)
        )
    return PushedForward(element.value_map, pulled_back)


class Argument(AbstractFunction):
    """A function that is a dimension of the tensor to be assembled."""

    def __init__(
        self,
        space: AbstractFunctionSpace,
        component: int,
        variable: AbstractVariable | None = None,
    ):
        """Initialise.

        Args:
            space: The function space that this function lives in
            component: The component of the finite element tensor
                       to be assembled that this function represents
            variable: The variable that is this argument's input
        """
        self._space = space
        self._variable = variable
        self._component = component

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""
        return self.__class__(self._space, self._component, variable)

    def reconstruct_with_space(self, space: AbstractFunctionSpace) -> Self:
        """Reconstruct this function in another function space."""
        return self.__class__(space, self._component, self._variable)

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space

    @property
    def variable(self) -> AbstractVariable | None:
        """Get the variable that is this function's input."""
        return self._variable

    @property
    def component_index(self) -> int:
        """The component of the finite element tensor that this function represents."""
        return self._component

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self._component, self._variable

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError()

    def get_replacement(self, replacements: dict[GraphNode, GraphNode]) -> GraphNode | None:
        """Get the node to replace this node with, or None if no replacement can be made."""
        for old, new in replacements.items():
            if (
                isinstance(old, Argument)
                and old.function_space == self.function_space
                and old.component_index == self.component_index
            ):
                if isinstance(new, Argument) and self.variable is not None and new.variable is None:
                    return new.reconstruct_with_variable(self.variable)
                return new

    def pull_back_to_entity(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the entity's coordinates."""
        return _pulled_back_to_entity(self)


class Coefficient(AbstractFunction):
    """A coefficient.

    A Coefficient represents a known function, such as a previous solution, a
    material property, or any other field supplied at assembly time.
    """

    _n = count(0)

    def __init__(
        self,
        space: AbstractFunctionSpace,
        coefficient_label: str | None = None,
        variable: AbstractVariable | None = None,
    ):
        """Initialise.

        Args:
            space: The function space that this function lives in
            coefficient_label: The label for this coefficient
            variable: The variable that is this argument's input
        """
        self._space = space
        self._variable = variable
        if coefficient_label is None:
            self._label = f"coefficient-{next(self._n)}"
        else:
            self._label = coefficient_label

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""
        return self.__class__(self._space, self._label, variable)

    def reconstruct_with_space(self, space: AbstractFunctionSpace) -> Self:
        """Reconstruct this function in another function space."""
        return self.__class__(space, self._label, self._variable)

    @property
    def function_space(self) -> AbstractFunctionSpace:
        """The function space that this function lives in."""
        return self._space

    @property
    def label(self) -> str:
        """The unique label of this coefficient."""
        return self._label

    @property
    def variable(self) -> AbstractVariable | None:
        """Get the variable that is this function's input."""
        return self._variable

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self._label, self._variable

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        if not 0 <= index < self.domain_size:
            raise ValueError(
                f"Derivative index {index} out of range for domain size {self.domain_size}"
            )
        if self.is_cellwise_constant:
            # A cellwise constant's derivative is a plain zero
            # Tensor/RealScalar -- not itself an AbstractFunction -- so
            # this is the one place that distinction has to be cast away
            # rather than widening AbstractFunction.diff's own contract
            # (which would break chained .diff().diff() calls elsewhere,
            # eg test_basis_functions.py, whose intermediate values are
            # statically typed as bare AbstractFunction).
            return cast(AbstractFunction, zero(self.value_shape))
        raise NotImplementedError()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError()

    def get_replacement(self, replacements: dict[GraphNode, GraphNode]) -> GraphNode | None:
        """Get the node to replace this node with, or None if no replacement can be made."""
        for old, new in replacements.items():
            if (
                isinstance(old, Coefficient)
                and old.function_space == self.function_space
                and old.label == self.label
            ):
                if (
                    isinstance(new, Coefficient)
                    and self.variable is not None
                    and new.variable is None
                ):
                    return new.reconstruct_with_variable(self.variable)
                return new

    def pull_back_to_entity(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the entity's coordinates."""
        return _pulled_back_to_entity(self)


class TestFunction(Argument):
    """A test function."""

    __test__ = False

    def __init__(
        self,
        space: AbstractFunctionSpace,
        variable: AbstractVariable | None = None,
    ):
        """Initialise."""
        super().__init__(space, 0, variable)

    def __repr__(self) -> str:
        """Representation."""
        return f"TestFunction({self.variable!r})"

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""
        return self.__class__(self._space, variable)

    def reconstruct_with_space(self, space: AbstractFunctionSpace) -> Self:
        """Reconstruct this function in another function space."""
        return self.__class__(space, self._variable)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self.variable

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        raise NotImplementedError()


class TrialFunction(Argument):
    """A trial function."""

    def __init__(
        self,
        space: AbstractFunctionSpace,
        variable: AbstractVariable | None = None,
    ):
        """Initialise."""
        super().__init__(space, 1, variable)

    def __repr__(self) -> str:
        """Representation."""
        return f"TrialFunction({self.variable!r})"

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""
        return self.__class__(self._space, variable)

    def reconstruct_with_space(self, space: AbstractFunctionSpace) -> Self:
        """Reconstruct this function in another function space."""
        return self.__class__(space, self._variable)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self.variable

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        raise NotImplementedError()


def _domains_read_on(node: GraphNode) -> list[AbstractDomain]:
    """The domains the parts of an expression say they are read on.

    A node that says is not looked inside: a restriction says the domain
    it restricts to, and the field under it is read through a map, which
    is the restriction's business and not this one's.
    """
    if isinstance(node, ReadOnADomain):
        reads_on = node.reads_on
        if reads_on is not None:
            return [reads_on]
    found: list[AbstractDomain] = []
    for successor in node.successors:
        for domain in _domains_read_on(successor):
            if domain not in found:
                found.append(domain)
    return found


def extract_domain(expression: GraphNode) -> AbstractDomain:
    """Get the domain an expression is read on.

    Read off the parts of it that say: a function from its space, a
    geometric quantity from the domain it differentiates, a restriction
    from the domain it restricts to. A part that has been positioned
    already says nothing -- a function pulled back onto a cell, a
    quantity told where it is read -- and arithmetic and tensors have no
    domain of their own, so they are asked of their successors.

    Args:
        expression: The expression to read

    Returns:
        The domain its parts agree on

    Raises:
        ValueError: If no part of it says which domain, or they disagree
    """
    found = _domains_read_on(expression)
    if len(found) == 0:
        raise ValueError(
            f"No part of {expression!r} says which domain it is read on. A domain is "
            f"read off the functions, geometric quantities and restrictions in an "
            f"expression, and this has none of them."
        )
    if len(found) > 1:
        raise ValueError(
            f"{expression!r} is read on {len(found)} domains, {found!r}. An expression "
            f"spanning two of them reaches the second through a map, which is what a "
            f"restriction is."
        )
    (domain,) = found
    return domain
