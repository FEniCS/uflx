"""Push forward and pull back maps.

A push forward and a pull back act on fields. A field on an entity is pushed forward to
one on the ambient coordinates, and a field on the ambient coordinates is pulled back to
one on the entity, so which direction a value map runs in is the whole of what it says.

Which map a field needs is fixed by what kind of field it is, which in differential
geometry is the degree of the form it represents, and each case is one of the Jacobian
quantities `uflx.geometry` names:

- a 0-form, a function, pushes forward by composition alone, so the map is the identity.
  This is the H1 conforming case, Lagrange and its relatives;
- a 1-form pushes forward by ``J^-T``, which is the covariant Piola map. This is the
  H(curl) conforming case, Nedelec and its relatives;
- an (n-1)-form pushes forward by ``J / det J``, which is the contravariant Piola map.
  This is the H(div) conforming case, Raviart-Thomas and its relatives. The determinant
  here is the signed one, :class:`uflx.geometry.JacobianDeterminant`, not the density
  :class:`uflx.geometry.VolumeElement`, which is why this map and not the others needs to
  know whether a chart inverts.

UFLx names the direction and leaves the expression to the element, so the maps here are
abstract and a consumer supplies the Piola it needs.

Applying a chart to a point is a different operation, being an image rather than a push
forward; `uflx.geometry`'s :class:`uflx.geometry.ImagePoint` and
:class:`uflx.geometry.PreimagePoint` are named for that.
"""

from abc import ABC, abstractmethod
from typing import Any, Protocol, runtime_checkable

from uflx.algorithms import replace
from uflx.expressions import AbstractExpression
from uflx.graphs import GraphNode, as_graph


class AbstractValueMap(ABC):
    """Abstract base class for value maps."""

    @abstractmethod
    def push_forward(self, function: AbstractExpression) -> AbstractExpression:
        """Map values from an entity's coordinates to ambient coordinates."""

    @abstractmethod
    def pull_back(self, function: AbstractExpression) -> AbstractExpression:
        """Map values from ambient coordinates to an entity's coordinates."""

    @abstractmethod
    def ambient_value_shape(
        self, entity_value_shape: tuple[int, ...], geometric_dimension: int
    ) -> tuple[int, ...]:
        """The shape a value takes once pushed forward.

        Args:
            entity_value_shape: The shape a value has on the entity
            geometric_dimension: The number of ambient coordinates

        Returns:
            The shape the pushed forward value has
        """

    @property
    def preserves_constant_values(self) -> bool:
        """Whether a function constant in an entity's coordinates stays constant once mapped.

        Conservative by default (False): a subclass overrides this only
        when it is known to preserve cellwise-constant values.

        """
        return False

    @property
    def is_identity(self) -> bool:
        """Whether this map leaves values untouched.

        Conservative by default (False): a subclass overrides this only
        when its push forward and pull back are both the identity.
        """
        return False


class IdentityValueMap(AbstractValueMap):
    """Identity map."""

    def push_forward(self, function: AbstractExpression) -> AbstractExpression:
        """Map values from an entity's coordinates to ambient coordinates."""
        return function

    def pull_back(self, function: AbstractExpression) -> AbstractExpression:
        """Map values from ambient coordinates to an entity's coordinates."""
        return function

    def ambient_value_shape(
        self, entity_value_shape: tuple[int, ...], geometric_dimension: int
    ) -> tuple[int, ...]:
        """The identity map leaves the shape alone."""
        return entity_value_shape

    @property
    def preserves_constant_values(self) -> bool:
        """The identity map trivially preserves constant values."""
        return True

    @property
    def is_identity(self) -> bool:
        """The identity map leaves values untouched."""
        return True


class BlockedValueMap(AbstractValueMap):
    """Map for blocked element."""

    def __init__(
        self,
        component_map: AbstractValueMap,
        shape: tuple[int, ...],
    ):
        """Initialise."""
        self._component_map = component_map
        self._shape = shape

    def push_forward(self, function: AbstractExpression) -> AbstractExpression:
        """Map values from an entity's coordinates to ambient coordinates."""
        return function  # TODO

    def pull_back(self, function: AbstractExpression) -> AbstractExpression:
        """Map values from ambient coordinates to an entity's coordinates."""
        return function  # TODO

    def ambient_value_shape(
        self, entity_value_shape: tuple[int, ...], geometric_dimension: int
    ) -> tuple[int, ...]:
        """A blocked value takes the shape of its block."""
        return self._shape

    @property
    def preserves_constant_values(self) -> bool:
        """A blocked map preserves constants iff its component map does."""
        return self._component_map.preserves_constant_values

    @property
    def is_identity(self) -> bool:
        """A blocked map leaves values untouched iff its component map does."""
        return self._component_map.is_identity


class SymmetricValueMap(AbstractValueMap):
    """Symmetric map."""

    def __init__(
        self,
        component_map: AbstractValueMap,
        shape: tuple[int, ...],
        symmetry_map: dict[tuple[int, ...], int],
    ):
        """Initialise."""
        self._component_map = component_map
        self._shape = shape
        self._symmetry_map = symmetry_map

    def push_forward(self, function: AbstractExpression) -> AbstractExpression:
        """Map values from an entity's coordinates to ambient coordinates."""
        return function  # TODO

    def pull_back(self, function: AbstractExpression) -> AbstractExpression:
        """Map values from ambient coordinates to an entity's coordinates."""
        return function  # TODO

    def ambient_value_shape(
        self, entity_value_shape: tuple[int, ...], geometric_dimension: int
    ) -> tuple[int, ...]:
        """A symmetric value takes the shape it is stored with."""
        return self._shape

    @property
    def preserves_constant_values(self) -> bool:
        """A symmetric map preserves constants iff its component map does."""
        return self._component_map.preserves_constant_values


class MixedValueMap(AbstractValueMap):
    """Map for a mixed element."""

    def __init__(
        self,
        sub_maps: list[AbstractValueMap],
        shapes: list[tuple[int, ...]],
    ):
        """Initialise."""
        self._sub_maps = sub_maps
        self._shapes = shapes

    def push_forward(self, function: AbstractExpression) -> AbstractExpression:
        """Map values from an entity's coordinates to ambient coordinates."""
        return function  # TODO

    def pull_back(self, function: AbstractExpression) -> AbstractExpression:
        """Map values from ambient coordinates to an entity's coordinates."""
        return function  # TODO

    def ambient_value_shape(
        self, entity_value_shape: tuple[int, ...], geometric_dimension: int
    ) -> tuple[int, ...]:
        """A mixed value's shape is its sub-maps' shapes laid end to end."""
        shape: tuple[int, ...] = ()
        for s in self._shapes:
            shape += s
        return shape

    @property
    def preserves_constant_values(self) -> bool:
        """A mixed map preserves constants iff every sub-map does."""
        return all(m.preserves_constant_values for m in self._sub_maps)


@runtime_checkable
class IsPushedForward(Protocol):
    """An object that has been pushed forward."""

    def apply_push_forward(self) -> GraphNode:
        """Apply the push forward."""


@runtime_checkable
class IsPulledBack(Protocol):
    """An object that has been pulled back."""

    def apply_pull_back(self) -> GraphNode:
        """Apply the pull back."""


class PushedForward(AbstractExpression):
    """A function in an entity's coordinates that has been mapped to ambient coordinates."""

    def __init__(self, map: AbstractValueMap, function: AbstractExpression):
        """Initialise."""
        self.map = map
        self.function = function

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.function.value_shape

    def __repr__(self):
        """Representation."""
        return f"PushedForward({self.map})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.function}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.map, self.function

    def apply_push_forward(self) -> AbstractExpression:
        """Apply the push forward."""
        return self.map.push_forward(self.function)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError("Cannot take a 'component' of a PushedForward object")


class PulledBack(AbstractExpression):
    """A function in ambient coordinates that has been mapped to an entity's coordinates."""

    def __init__(self, map: AbstractValueMap, function: AbstractExpression):
        """Initalise."""
        self.map = map
        self.function = function

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self.function.value_shape

    def __repr__(self):
        """Representation."""
        return f"PulledBack({self.map})"

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self.function}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.map, self.function

    def apply_pull_back(self) -> AbstractExpression:
        """Apply the pull back."""
        return self.map.pull_back(self.function)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise NotImplementedError("Cannot take a 'component' of a PulledBack object")


def apply_push_forwards(
    expression: GraphNode,
) -> GraphNode:
    """Apply push forward maps to functions."""
    return replace(
        expression,
        {
            node: node.apply_push_forward()
            for node in as_graph(expression)
            if isinstance(node, IsPushedForward) and isinstance(node, GraphNode)
        },
    )


def apply_pull_backs(
    expression: GraphNode,
) -> GraphNode:
    """Apply pull back maps to functions."""
    return replace(
        expression,
        {
            node: node.apply_pull_back()
            for node in as_graph(expression)
            if isinstance(node, IsPulledBack) and isinstance(node, GraphNode)
        },
    )
