"""Geometry."""

from typing import Any, Protocol, Self, runtime_checkable

from uflx.algorithms import replace
from uflx.domains import RD, AbstractCoordinateDomain, AbstractParametrization, EntityDomain
from uflx.expressions import AbstractExpression
from uflx.functions import AbstractVariable
from uflx.graphs import GraphNode, as_graph
from uflx.points import AbstractPoint, Point
from uflx.tensors import Identity, Matrix


@runtime_checkable
class ExpandableGeometry(Protocol):
    """Geometry that can be expanded into components."""

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""


class SingleSpatialCoordinate(AbstractExpression):
    """A variable representing a component of R^d."""

    def __init__(self, dimension: int, component: int):
        """Initialise."""
        self._dimension = dimension
        self._component = component

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._dimension, self._component

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class SpatialCoordinate(AbstractExpression):
    """A variable on R^d."""

    def __init__(self, dimension: int):
        """Initialise."""
        self._dimension = dimension

    def __getitem__(self, component: int) -> SingleSpatialCoordinate:
        """Get item."""
        if component < 0 or component >= self._dimension:
            raise IndexError("coordinate index out of range")
        return SingleSpatialCoordinate(self._dimension, component)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self._dimension,)

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return (self._dimension,)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        (i,) = indices
        return SingleSpatialCoordinate(self._dimension, i)


def as_matrix(jacobian: AbstractExpression) -> Matrix:
    """Densify a Jacobian so a determinant, inverse or transpose can be taken.

    An identity map's Jacobian is an Identity, which stays symbolic in a
    matrix product but has to be written out to be inverted.
    """
    if isinstance(jacobian, Identity):
        size = jacobian.size
        return Matrix([[jacobian.component(i, j) for j in range(size)] for i in range(size)])
    assert isinstance(jacobian, Matrix)
    return jacobian


class PushedForwardPoint(AbstractPoint):
    """A point in an entity's coordinates, mapped through a parametrization."""

    def __init__(self, point: AbstractPoint, parametrization: AbstractParametrization):
        """Initialise.

        Args:
            point: The point in the entity's coordinates
            parametrization: The map to push the point forward through
        """
        assert isinstance(point.domain, EntityDomain)
        self._point = point
        self._parametrization = parametrization

    @property
    def entity_point(self) -> AbstractPoint:
        """The point in the entity's coordinates."""
        return self._point

    @property
    def parametrization(self) -> AbstractParametrization:
        """The map this point is pushed forward through."""
        return self._parametrization

    @property
    def domain(self) -> RD:
        """The domain.

        A pushed forward point lies in the coordinates the map lands in.
        """
        return RD(self._parametrization.target_dimension)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self.dim,)

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        return self._parametrization.target_dimension

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._point}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._point, self._parametrization

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        value = self._parametrization.value(self._point)

        # The map's values are ambient coordinates, so what comes out is a
        # point of R^gdim rather than one of the domain it parametrizes.
        return Point([value.component(j) for j in range(self.dim)], self.domain)

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, PushedForwardPoint)
            and self._point == other._point
            and self._parametrization == other._parametrization
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.PushedForwardPoint", self._point, self._parametrization))

    @property
    def index(self) -> int | str:
        """The point's index in the set of points."""
        return self._point.index


class PulledBackPoint(AbstractPoint):
    """A point in ambient coordinates, mapped to an entity's coordinates."""

    def __init__(self, point: AbstractPoint, parametrization: AbstractParametrization):
        """Initialise.

        Args:
            point: The point in ambient coordinates
            parametrization: The map the point is pulled back through
        """
        assert not isinstance(point.domain, EntityDomain)
        self._point = point
        self._parametrization = parametrization

    @property
    def ambient_point(self) -> AbstractPoint:
        """The point in ambient coordinates."""
        return self._point

    @property
    def domain(self) -> AbstractCoordinateDomain:
        """The domain.

        A pulled back point lies in the coordinates the map starts from,
        not in the ones it lands in.
        """
        return self._parametrization.source

    @property
    def parametrization(self) -> AbstractParametrization:
        """The map this point is pulled back through."""
        return self._parametrization

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self.dim,)

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        return self._parametrization.source.geometric_dimension

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._point}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._point, self._parametrization

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, PulledBackPoint)
            and self._point == other._point
            and self._parametrization == other._parametrization
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.PulledBackPoint", self._point, self._parametrization))

    @property
    def index(self) -> int | str:
        """The point's index in the set of points."""
        return self._point.index


class AbstractJacobian(AbstractExpression):
    """Base class for the derivative of a parametrization, and quantities built on it.

    A Jacobian is evaluated somewhere, but when one is created during a
    pull back there is no point to evaluate it at yet. It is therefore
    built without one and told later, by the integral whose dummy
    variable stands for the point.
    """

    def __init__(
        self, parametrization: AbstractParametrization, point: AbstractVariable | None = None
    ):
        """Initialise.

        Args:
            parametrization: The map being differentiated
            point: Where to differentiate it, if that is known yet
        """
        self.parametrization = parametrization
        self.point = point

    @property
    def _jacobian(self) -> "Jacobian":
        """The Jacobian this quantity is built from."""
        return Jacobian(self.parametrization, self.point)

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.parametrization, self.point

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Evaluate this quantity at the given variable.

        The variable stands for a point of the map's source, so one in
        any other coordinates is not this quantity's to take.
        """
        if variable.domain != self.parametrization.source:
            return self
        return self.__class__(self.parametrization, variable)


class Jacobian(AbstractJacobian):
    """The Jacobian."""

    @property
    def _jacobian(self) -> "Jacobian":
        """The Jacobian is its own."""
        return self

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (
            self.parametrization.target_dimension,
            self.parametrization.source.geometric_dimension,
        )

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        assert self.point is not None
        return self.parametrization.jacobian(self.point)

    def __repr__(self) -> str:
        """Representation."""
        return f"Jacobian({self.parametrization!r}, {self.point!r})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        if (
            isinstance(other, JacobianInverse)
            and self.parametrization == other.parametrization
            and self.point == other.point
        ):
            return Identity(self.value_shape[0])


class JacobianDeterminant(AbstractJacobian):
    """The determinant of the Jacobian."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = as_matrix(self._jacobian.expand_geometry())
        return abs(j.compute_determinant())

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class JacobianInverse(AbstractJacobian):
    """The inverse of the Jacobian."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._jacobian.value_shape[::-1]

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = as_matrix(self._jacobian.expand_geometry())
        return j.compute_inverse()

    def __repr__(self) -> str:
        """Representation."""
        return f"JacobianInverse({self.parametrization!r}, {self.point!r})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        if (
            isinstance(other, Jacobian)
            and self.parametrization == other.parametrization
            and self.point == other.point
        ):
            return Identity(self.value_shape[0])


class JacobianTranspose(AbstractJacobian):
    """The transpose of the Jacobian."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._jacobian.value_shape[::-1]

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = as_matrix(self._jacobian.expand_geometry())
        return j.transpose()

    def __repr__(self) -> str:
        """Representation."""
        return f"JacobianTranspose({self.parametrization!r}, {self.point!r})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        if (
            isinstance(other, JacobianInverseTranspose)
            and self.parametrization == other.parametrization
            and self.point == other.point
        ):
            return Identity(self.value_shape[0])


class JacobianInverseTranspose(AbstractJacobian):
    """The inverse transpose of the Jacobian."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._jacobian.value_shape

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = as_matrix(self._jacobian.expand_geometry())
        return j.compute_inverse().transpose()

    def __repr__(self) -> str:
        """Representation."""
        return f"JacobianInverseTranspose({self.parametrization!r}, {self.point!r})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        if (
            isinstance(other, JacobianTranspose)
            and self.parametrization == other.parametrization
            and self.point == other.point
        ):
            return Identity(self.value_shape[0])


def expand_geometry(
    expression: GraphNode,
) -> GraphNode:
    """Replace jacobians with evaluations of the derivatives of finite elements."""
    to_replace: dict[GraphNode, GraphNode] = {}

    for node in as_graph(expression):
        if isinstance(node, GraphNode) and isinstance(node, ExpandableGeometry):
            to_replace[node] = node.expand_geometry()

    return replace(expression, to_replace)
