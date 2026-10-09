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

from uflx.domains import AbstractDomain, AbstractFiniteElementDomain, EntityDomain
from uflx.expressions import AbstractExpression, Im, Re
from uflx.function_spaces import AbstractFunctionSpace, AbstractMappedFunctionSpace
from uflx.graphs import GraphNode
from uflx.maps import PushedForward
from uflx.tensors import zero


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
    def in_entity_coordinates(self) -> bool:
        """Check if this variable's components are in an entity's coordinates."""
        return isinstance(self.domain, EntityDomain)

    def to_entity_coordinates(self) -> FiniteElementVariable:
        """Make a version of this variable in the entity's coordinates."""
        raise ValueError("Cannot pull this variable back to the entity's coordinates.")


class Variable(AbstractVariable):
    """A variable that is the input to a function."""

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

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return isinstance(other, Variable) and other.label == self.label

    def __repr__(self) -> str:
        """Representation."""
        return f"Variable({self._label})"

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.Variable", self._label))

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._domain, self._label)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError()


class FiniteElementVariable(AbstractVariable):
    """A variable that is the input to a function."""

    _n = count(0)

    def __init__(
        self,
        domain: AbstractFiniteElementDomain,
        label: str | None = None,
        in_entity_coordinates: bool = False,
    ):
        """Initialise."""
        if label is None:
            self._label = f"variable-{next(self._n)}"
        else:
            self._label = label
        self._domain = domain
        self._in_entity_coordinates = in_entity_coordinates

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._domain, self._label, self._in_entity_coordinates)

    @property
    def label(self) -> str:
        """The label of this variable."""
        return self._label

    @property
    def domain(self) -> AbstractDomain:
        """The domain that this variable is in."""
        return self._domain

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, FiniteElementVariable)
            and other.label == self.label
            and other.in_entity_coordinates == self.in_entity_coordinates
        )

    def __repr__(self) -> str:
        """Representation."""
        if self._in_entity_coordinates:
            return f"FiniteElementVariable({self._label}, in_entity_coordinates=True)"
        else:
            return f"FiniteElementVariable({self._label})"

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.FiniteElementVariable", self._label, self._in_entity_coordinates))

    @property
    def in_entity_coordinates(self) -> bool:
        """Check if this variable's components are in an entity's coordinates.

        Stored rather than derived from the domain, because
        to_entity_coordinates keeps the same domain for now.
        """
        return self._in_entity_coordinates

    def to_entity_coordinates(self) -> FiniteElementVariable:
        """Make a version of this variable in the entity's coordinates."""
        return FiniteElementVariable(self._domain, self._label, True)

    def to_ambient_coordinates(self) -> FiniteElementVariable:
        """Make a version of this variable in ambient coordinates."""
        return FiniteElementVariable(self._domain, self._label, False)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError()


class AbstractFunction(AbstractExpression):
    """Base class for a function."""

    @property
    @abstractmethod
    def variable(self) -> AbstractVariable | None:
        """Get the variable that is this function's input."""

    @property
    def in_entity_coordinates(self) -> bool:
        """Check if this function's components are in an entity's coordinates."""
        if self.variable is not None:
            return self.variable.in_entity_coordinates
        return False

    @abstractmethod
    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""

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
        if self.in_entity_coordinates:
            assert isinstance(self.function_space, AbstractMappedFunctionSpace)
            return self.function_space.elements[0].entity_value_shape
        else:
            return self.function_space.value_shape

    @property
    def is_cellwise_constant(self) -> bool:
        """Whether this function's value is the same everywhere on every cell.

        True when every element of the function space is known to lie in
        the degree-0 Lagrange space (lagrange_superdegree == 0) and,
        for a function in ambient coordinates, every element's value
        map is known to preserve that constancy when pushed forward
        (see AbstractValueMap.preserves_constant_values).
        """
        assert isinstance(self.function_space, AbstractMappedFunctionSpace)
        elements = self.function_space.elements
        if any(e.lagrange_superdegree != 0 for e in elements):
            return False
        if self.in_entity_coordinates:
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


class Argument(AbstractFunction):
    """A function that is a dimension of the tensor to be assembled."""

    def __init__(
        self,
        space: AbstractFunctionSpace,
        component: int,
        in_entity_coordinates: bool = False,
        variable: AbstractVariable | None = None,
    ):
        """Initialise.

        Args:
            space: The function space that this function lives in
            component: The component of the finite element tensor
                       to be assembled that this function represents
            in_entity_coordinates: Are this argument's components in the entity's coordinates?
            variable: The variable that is this argument's input
        """
        if variable is not None:
            assert in_entity_coordinates == variable.in_entity_coordinates
        self._space = space
        self._in_entity_coordinates = in_entity_coordinates
        self._variable = variable
        self._component = component

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""
        return self.__class__(self._space, self._component, self._in_entity_coordinates, variable)

    @property
    def in_entity_coordinates(self) -> bool:
        """Check if this function's components are in an entity's coordinates."""
        return self._in_entity_coordinates

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
        return self._space, self._component, self._in_entity_coordinates, self._variable

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
        in_entity_coordinates: bool = False,
        variable: AbstractVariable | None = None,
    ):
        """Initialise.

        Args:
            space: The function space that this function lives in
            coefficient_label: The label for this coefficient
            in_entity_coordinates: Are this argument's components in the entity's coordinates?
            variable: The variable that is this argument's input
        """
        if variable is not None:
            if isinstance(variable, FiniteElementVariable):
                assert in_entity_coordinates == variable.in_entity_coordinates
            else:
                assert not in_entity_coordinates
        self._space = space
        self._in_entity_coordinates = in_entity_coordinates
        self._variable = variable
        if coefficient_label is None:
            self._label = f"coefficient-{next(self._n)}"
        else:
            self._label = coefficient_label

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""
        return self.__class__(self._space, self._label, self._in_entity_coordinates, variable)

    @property
    def in_entity_coordinates(self) -> bool:
        """Check if this function's components are in an entity's coordinates."""
        return self._in_entity_coordinates

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
        return self._space, self._label, self._in_entity_coordinates, self._variable

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
        if self.in_entity_coordinates:
            raise ValueError("Cannot pull back a function already in the entity's coordinates")
        assert isinstance(self._space, AbstractMappedFunctionSpace)
        if self._variable is None:
            return PushedForward(
                self._space.elements[0].value_map,
                Coefficient(self._space, self._label, True, None),
            )
        else:
            assert isinstance(self._variable, FiniteElementVariable)
            return PushedForward(
                self._space.elements[0].value_map,
                Coefficient(self._space, self._label, True, self._variable.to_entity_coordinates()),
            )


class TestFunction(Argument):
    """A test function."""

    __test__ = False

    def __init__(
        self,
        space: AbstractFunctionSpace,
        in_entity_coordinates: bool = False,
        variable: AbstractVariable | None = None,
    ):
        """Initialise."""
        super().__init__(space, 0, in_entity_coordinates, variable)

    def __repr__(self) -> str:
        """Representation."""
        return f"TestFunction({self.variable!r})"

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""
        return self.__class__(self._space, self._in_entity_coordinates, variable)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self.in_entity_coordinates, self.variable

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        raise NotImplementedError()

    def pull_back_to_entity(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the entity's coordinates."""
        if self.in_entity_coordinates:
            raise ValueError("Cannot pull back a function already in the entity's coordinates")
        assert isinstance(self._space, AbstractMappedFunctionSpace)
        if self.variable is None:
            return PushedForward(
                self._space.elements[0].value_map,
                TestFunction(self._space, True, None),
            )
        else:
            assert isinstance(self.variable, FiniteElementVariable)
            return PushedForward(
                self._space.elements[0].value_map,
                TestFunction(self._space, True, self.variable.to_entity_coordinates()),
            )


class TrialFunction(Argument):
    """A trial function."""

    def __init__(
        self,
        space: AbstractFunctionSpace,
        in_entity_coordinates: bool = False,
        variable: AbstractVariable | None = None,
    ):
        """Initialise."""
        super().__init__(space, 1, in_entity_coordinates, variable)

    def __repr__(self) -> str:
        """Representation."""
        return f"TrialFunction({self.variable!r})"

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Reconstruct this function taking the input variable as input."""
        return self.__class__(self._space, self._in_entity_coordinates, variable)

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._space, self.in_entity_coordinates, self.variable

    def diff(self, index: int) -> AbstractFunction:
        """Take a derivative of this function."""
        raise NotImplementedError()

    def pull_back_to_entity(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the entity's coordinates."""
        if self.in_entity_coordinates:
            raise ValueError("Cannot pull back a function already in the entity's coordinates")
        assert isinstance(self._space, AbstractMappedFunctionSpace)
        if self.variable is None:
            return PushedForward(
                self._space.elements[0].value_map,
                TrialFunction(self._space, True, None),
            )
        else:
            assert isinstance(self.variable, FiniteElementVariable)
            return PushedForward(
                self._space.elements[0].value_map,
                TrialFunction(self._space, True, self.variable.to_entity_coordinates()),
            )


def create_variable(domain: AbstractDomain) -> AbstractVariable:
    """Create a new variable in a domain."""
    if isinstance(domain, AbstractFiniteElementDomain):
        return FiniteElementVariable(domain)
    else:
        return Variable(domain)
