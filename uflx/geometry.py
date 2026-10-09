"""Geometry."""

from typing import Any, Protocol, runtime_checkable

from uflx.algorithms import replace
from uflx.basis_functions import EvaluatedBasisFunction
from uflx.domains import RD, AbstractParametrizedDomain, EntityDomain
from uflx.expressions import AbstractExpression, expression_sum
from uflx.function_spaces import function_space
from uflx.graphs import GraphNode, as_graph
from uflx.points import AbstractPoint, Point
from uflx.tensors import FlattenedTensorMap, Identity, Matrix


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


class PushedForwardPoint(AbstractPoint):
    """A point in an entity's coordinates, mapped to ambient coordinates."""

    def __init__(self, point: AbstractPoint, domain: AbstractParametrizedDomain):
        """Initialise."""
        assert isinstance(point.domain, EntityDomain)
        self._point = point
        self._domain = domain

    @property
    def entity_point(self) -> AbstractPoint:
        """The point in the entity's coordinates."""
        return self._point

    @property
    def domain(self) -> AbstractParametrizedDomain:
        """The domain."""
        return self._domain

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self._domain.geometric_dimension,)

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        return self._domain.geometric_dimension

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._point}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._point, self._domain

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        if len(self.domain.cell_types) != 1:
            raise NotImplementedError("Only domains with exactly on element supported for now.")
        element = self.domain.parametrization(self.domain.cell_types[0])
        (dim,) = element.entity_value_shape

        components = [
            expression_sum(
                FlattenedTensorMap((i // dim, i % dim), (dim,))
                * EvaluatedBasisFunction(
                    function_space(self.domain, element), i, self.entity_point, component=j
                )
                for i in range(element.dim)
            )
            for j in range(dim)
        ]

        # The expansion is an interpolation sum of coordinate dofs, so what
        # comes out is explicit ambient coordinates.
        return Point(components, RD(dim))

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, PushedForwardPoint)
            and self._point == other._point
            and self._domain == other._domain
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.PushedForwardPoint", hash(self._point), hash(self._domain)))

    @property
    def index(self) -> int | str:
        """The point's index in the set of points."""
        return self._point.index


class PulledBackPoint(AbstractPoint):
    """A point in ambient coordinates, mapped to an entity's coordinates."""

    def __init__(self, point: AbstractPoint, domain: AbstractParametrizedDomain):
        """Initialise.

        Args:
            point: The point in ambient coordinates
            domain: The parametrized domain the point is pulled back through
        """
        assert not isinstance(point.domain, EntityDomain)
        self._point = point
        self._parametrized_domain = domain

    @property
    def ambient_point(self) -> AbstractPoint:
        """The point in ambient coordinates."""
        return self._point

    @property
    def domain(self) -> EntityDomain:
        """The domain.

        A pulled back point lies in the entity's coordinates, not on the
        parametrized domain it came from.
        """
        (cell,) = self._parametrized_domain.cell_types
        return EntityDomain(cell)

    @property
    def parametrized_domain(self) -> AbstractParametrizedDomain:
        """The parametrized domain that this point was pulled back through."""
        return self._parametrized_domain

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self.dim,)

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        tdim = self._parametrized_domain.topological_dimension
        if tdim is None:
            raise NotImplementedError(
                "Points in a cell's coordinates are not supported on domains with "
                "cells of several topological dimensions."
            )
        return tdim

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._point}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._point, self._parametrized_domain

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, PulledBackPoint)
            and self._point == other._point
            and self._parametrized_domain == other._parametrized_domain
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.PulledBackPoint", self._point, self._parametrized_domain))

    @property
    def index(self) -> int | str:
        """The point's index in the set of points."""
        return self._point.index


class Jacobian(AbstractExpression):
    """The Jacobian."""

    def __init__(self, domain: AbstractParametrizedDomain, point: AbstractPoint | None = None):
        """Initalise."""
        self.domain = domain
        self.point = point

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        assert self.domain.topological_dimension is not None
        return (self.domain.geometric_dimension, self.domain.topological_dimension)

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.domain, self.point

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        gdim, tdim = self.value_shape
        if len(self.domain.cell_types) > 1:
            raise NotImplementedError()
        (cell,) = self.domain.cell_types
        element = self.domain.parametrization(cell)

        assert self.point is not None

        return Matrix(
            [
                [
                    expression_sum(
                        FlattenedTensorMap((i // gdim, i % gdim), (gdim,))
                        * EvaluatedBasisFunction(
                            function_space(self.domain, element),
                            i,
                            self.point,
                            derivative=tuple(1 if d == col else 0 for d in range(tdim)),
                            component=row,
                        )
                        for i in range(element.dim)
                    )
                    for col in range(tdim)
                ]
                for row in range(gdim)
            ]
        )

    def __repr__(self) -> str:
        """Representation."""
        return f"Jacobian({self.domain!r}, {self.point!r})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        if (
            isinstance(other, JacobianInverse)
            and self.domain == other.domain
            and self.point == other.point
        ):
            return Identity(self.value_shape[0])


class JacobianDeterminant(AbstractExpression):
    """The determinant of the Jacobian."""

    def __init__(self, domain: AbstractParametrizedDomain, point: AbstractPoint | None = None):
        """Initialise."""
        self._jacobian = Jacobian(domain, point)
        self.domain = domain
        self.point = point

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
        return self.domain, self.point

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = self._jacobian.expand_geometry()
        assert isinstance(j, Matrix)
        return abs(j.compute_determinant())

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class JacobianInverse(AbstractExpression):
    """The inverse of the Jacobian."""

    def __init__(self, domain: AbstractParametrizedDomain, point: AbstractPoint | None = None):
        """Initalise."""
        self._jacobian = Jacobian(domain, point)
        self.domain = domain
        self.point = point

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._jacobian.value_shape[::-1]

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.domain, self.point

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = self._jacobian.expand_geometry()
        assert isinstance(j, Matrix)
        return j.compute_inverse()

    def __repr__(self) -> str:
        """Representation."""
        return f"JacobianInverse({self.domain!r}, {self.point!r})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        if (
            isinstance(other, Jacobian)
            and self.domain == other.domain
            and self.point == other.point
        ):
            return Identity(self.value_shape[0])


class JacobianTranspose(AbstractExpression):
    """The transpose of the Jacobian."""

    def __init__(self, domain: AbstractParametrizedDomain, point: AbstractPoint | None = None):
        """Initalise."""
        self._jacobian = Jacobian(domain, point)
        self.domain = domain
        self.point = point

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._jacobian.value_shape[::-1]

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.domain, self.point

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = self._jacobian.expand_geometry()
        assert isinstance(j, Matrix)
        return j.transpose()

    def __repr__(self) -> str:
        """Representation."""
        return f"JacobianTranspose({self.domain!r}, {self.point!r})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        if (
            isinstance(other, JacobianInverseTranspose)
            and self.domain == other.domain
            and self.point == other.point
        ):
            return Identity(self.value_shape[0])


class JacobianInverseTranspose(AbstractExpression):
    """The inverse transpose of the Jacobian."""

    def __init__(self, domain: AbstractParametrizedDomain, point: AbstractPoint | None = None):
        """Initalise."""
        self._jacobian = Jacobian(domain, point)
        self.domain = domain
        self.point = point

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._jacobian.value_shape

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self.domain, self.point

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = self._jacobian.expand_geometry()
        assert isinstance(j, Matrix)
        return j.compute_inverse().transpose()

    def __repr__(self) -> str:
        """Representation."""
        return f"JacobianInverseTranspose({self.domain!r}, {self.point!r})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        if (
            isinstance(other, JacobianTranspose)
            and self.domain == other.domain
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
