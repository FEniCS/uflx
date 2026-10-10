"""Geometry."""

from __future__ import annotations

from typing import Any, Protocol, Self, runtime_checkable

from uflx.algorithms import replace
from uflx.domains import (
    RD,
    AbstractCoordinateDomain,
    AbstractParametrization,
    AbstractParametrizedDomain,
    EntityDomain,
)
from uflx.expressions import AbstractExpression, Sqrt, expression_sum
from uflx.functions import AbstractVariable
from uflx.graphs import GraphNode, as_graph
from uflx.points import AbstractPoint, Point
from uflx.tensors import Identity, Matrix, Vector


@runtime_checkable
class ExpandableGeometry(Protocol):
    """Geometry that can be expanded into components."""

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""


def _as_dense_matrix(jacobian: AbstractExpression) -> Matrix:
    """Write a Jacobian out entry by entry, so it can be inverted or reduced.

    A map's Jacobian need not already be a Matrix: an identity map gives
    an Identity, and a composed map gives the MatrixProduct of the chain
    rule. Both have to be written out before a determinant, an inverse or
    a transpose can be taken of them.
    """
    if isinstance(jacobian, Matrix):
        return jacobian
    rows, cols = jacobian.value_shape
    return Matrix([[jacobian.component(i, j) for j in range(cols)] for i in range(rows)])


class PushedForwardPoint(AbstractPoint):
    """A point in an entity's coordinates, mapped through a parametrization.

    The point level action of a map: the same thing
    :class:`SpatialCoordinate` gives as an expression, but being a point
    it can be the variable a basis function is evaluated at. An element
    defined on the physical cell rather than on a reference one needs
    that, having no map to compose from.

    It takes the domain rather than the map because the point it starts
    from already lies in a cell's coordinate domain, and so names the
    cell whose map carries it. Contrast :class:`PulledBackPoint`.
    """

    def __init__(self, point: AbstractPoint, domain: AbstractParametrizedDomain):
        """Initialise.

        Args:
            point: The point in the entity's coordinates
            domain: The domain to push the point forward onto

        Raises:
            ValueError: If the point is not in a cell's coordinates
        """
        if not isinstance(point.domain, EntityDomain):
            raise ValueError(
                f"A point is pushed forward from a cell's coordinates, not from {point.domain!r}."
            )
        self._point = point
        self._parametrized_domain = domain

    @property
    def entity_point(self) -> AbstractPoint:
        """The point in the entity's coordinates."""
        return self._point

    @property
    def parametrized_domain(self) -> AbstractParametrizedDomain:
        """The domain this point is pushed forward onto."""
        return self._parametrized_domain

    @property
    def parametrization(self) -> AbstractParametrization:
        """The map this point is pushed forward through.

        The point already lies in a cell's coordinate domain, so it names
        the cell whose map carries it.
        """
        source = self._point.domain
        assert isinstance(source, EntityDomain)
        (cell,) = source.cell_types
        return self._parametrized_domain.parametrization(cell)

    @property
    def domain(self) -> RD:
        """The domain.

        A pushed forward point lies in the coordinates the map lands in.
        """
        return RD(self._parametrized_domain.geometric_dimension)

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self.dim,)

    @property
    def dim(self) -> int:
        """The dimension of the point."""
        return self._parametrized_domain.geometric_dimension

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._point}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._point, self._parametrized_domain

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        # Expanding has to reach the bottom in one go, since the walk that
        # calls this does not revisit what it returns.
        value = SpatialCoordinate(self._parametrized_domain, self._point).expand_geometry()

        # The map's values are ambient coordinates, so what comes out is a
        # point of R^gdim rather than one of the domain it parametrizes.
        return Point([value.component(j) for j in range(self.dim)], self.domain)

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, PushedForwardPoint)
            and self._point == other._point
            and self._parametrized_domain == other._parametrized_domain
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.PushedForwardPoint", self._point, self._parametrized_domain))

    @property
    def index(self) -> int | str:
        """The point's index in the set of points."""
        return self._point.index


class PulledBackPoint(AbstractPoint):
    """A point in ambient coordinates, mapped to an entity's coordinates.

    A terminal, deliberately: it has no expansion and will not get one.
    A parametrization offers a value and a derivative and no inverse,
    because inverting a finite element map is a Newton solve rather than
    anything symbolic. So this names a point a consumer computes, the way
    a coordinate dof names a number a mesh holds.

    It takes the map rather than the domain because a point of the ambient
    coordinates does not say which cell it should land in; that is part of
    the question being asked. Contrast :class:`PushedForwardPoint`.

    What wants it is evaluation at a physical location: a reference basis
    read at ``phi^-1(x)``, which is what a point evaluation needs.
    """

    def __init__(self, point: AbstractPoint, parametrization: AbstractParametrization):
        """Initialise.

        Args:
            point: The point in ambient coordinates
            parametrization: The map the point is pulled back through, which says
                which cell's coordinates it lands in

        Raises:
            ValueError: If the point is already in a cell's coordinates
        """
        if isinstance(point.domain, EntityDomain):
            raise ValueError(
                f"A point in {point.domain!r} is already in a cell's coordinates, so "
                f"there is nothing to pull back."
            )
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


class AbstractGeometricQuantity(AbstractExpression):
    """Base class for a quantity read off a domain's parametrization at a point.

    A geometric quantity is evaluated somewhere, but when one is created
    during a pull back there is no point to evaluate it at yet. It is
    therefore built without one and told later, by the integral whose
    dummy variable stands for the point.
    """

    def __init__(self, domain: AbstractParametrizedDomain, point: AbstractVariable | None = None):
        """Initialise.

        Args:
            domain: The domain whose geometry is being differentiated
            point: Where to differentiate it, if that is known yet
        """
        self._domain = domain
        self._point = point

    @property
    def domain(self) -> AbstractParametrizedDomain:
        """The domain whose geometry this quantity differentiates."""
        return self._domain

    @property
    def point(self) -> AbstractVariable | None:
        """Where it is differentiated, if it has been told yet."""
        return self._point

    @property
    def parametrization(self) -> AbstractParametrization:
        """The map this quantity differentiates, for the cell its point lies in.

        A point lies in a cell's coordinate domain, so it names the cell.
        Until a point arrives this quantity is generic over the domain's
        cell types, which is what lets it be built during a pull back.
        """
        if self.point is None:
            raise ValueError(
                "This quantity has not been told where it is evaluated, so the cell "
                "whose map it differentiates is not known."
            )
        source = self.point.domain
        if not isinstance(source, EntityDomain):
            raise ValueError(
                f"This quantity is evaluated at a point of a cell's coordinate domain, "
                f"not of {source!r}."
            )
        (cell,) = source.cell_types
        return self.domain.parametrization(cell)

    @property
    def _jacobian(self) -> Jacobian:
        """The Jacobian this quantity is built from."""
        return Jacobian(self.domain, self.point)

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._domain, self._point

    def __repr__(self) -> str:
        """Representation."""
        return f"{self.__class__.__name__}({self._domain!r}, {self._point!r})"

    def _is_partnered_by(self, other: GraphNode, partner: type[AbstractGeometricQuantity]) -> bool:
        """Whether `other` is a quantity of the given kind on the same geometry.

        Args:
            other: The quantity this one is multiplied by
            partner: The class it must be for the product to simplify
        """
        return (
            isinstance(other, partner)
            and self._domain == other.domain
            and self._point == other.point
        )

    def _onto_the_tangent_space(self) -> GraphNode:
        """What a map composed with its own pseudo-inverse comes to.

        Only the identity when the map is square. Otherwise it projects the
        ambient coordinates onto the tangent space, keeping tdim of the
        gdim directions.
        """
        gdim, tdim = self._jacobian.value_shape
        if gdim == tdim:
            return Identity(gdim)
        return TangentialProjector(self._domain, self._point)

    def _at_point(self, point: AbstractVariable) -> Self:
        """This quantity at a point, rebuilt.

        A subclass carrying more than a domain and a point overrides this.
        """
        return self.__class__(self._domain, point)

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Evaluate this quantity at the given variable.

        The variable stands for a point of one of this domain's cells, so
        one in any other coordinates is not this quantity's to take.
        """
        source = variable.domain
        if not isinstance(source, EntityDomain):
            return self
        (cell,) = source.cell_types
        if cell not in self.domain.cell_types:
            return self
        return self._at_point(variable)


class SpatialCoordinate(AbstractGeometricQuantity):
    """The ambient coordinates a domain's parametrization lands on.

    ``x = phi(X)``, the map's own value: where a point of a cell's
    coordinate domain sits in the ambient coordinates. Writing a source
    term as ``sin(x[0])`` is what this is for.
    """

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self._domain.geometric_dimension,)

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        parametrization = self.parametrization
        assert self._point is not None
        return parametrization.value(self._point)

    def __getitem__(self, component: int) -> SingleSpatialCoordinate:
        """Get one coordinate."""
        return self.component(component)

    def component(self, *indices: int) -> SingleSpatialCoordinate:
        """Get a component of the expression."""
        (i,) = indices
        (dimension,) = self.value_shape
        if i < 0 or i >= dimension:
            raise IndexError("coordinate index out of range")
        return SingleSpatialCoordinate(self._domain, i, self._point)


class SingleSpatialCoordinate(AbstractGeometricQuantity):
    """One of the ambient coordinates a domain's parametrization lands on.

    Kept as a node of its own rather than expanded on sight, so that
    ``sin(x[0])`` stays readable until geometry is expanded.
    """

    def __init__(
        self,
        domain: AbstractParametrizedDomain,
        component: int,
        point: AbstractVariable | None = None,
    ):
        """Initialise.

        Args:
            domain: The domain whose coordinates these are
            component: Which coordinate, in ``range(geometric_dimension)``
            point: Where to evaluate it, if that is known yet
        """
        super().__init__(domain, point)
        self._component = component

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._domain, self._component, self._point

    def _at_point(self, point: AbstractVariable) -> Self:
        """This coordinate at a point, rebuilt."""
        return self.__class__(self._domain, self._component, point)

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        return (
            SpatialCoordinate(self._domain, self._point)
            .expand_geometry()
            .component(self._component)
        )

    def __repr__(self) -> str:
        """Representation."""
        return f"SingleSpatialCoordinate({self._domain!r}, {self._component}, {self._point!r})"

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class Jacobian(AbstractGeometricQuantity):
    """The Jacobian."""

    @property
    def _jacobian(self) -> Jacobian:
        """The Jacobian is its own."""
        return self

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        tdim = self.domain.topological_dimension
        if tdim is None:
            raise NotImplementedError(
                "A Jacobian is not supported on a domain whose cells have several "
                "topological dimensions."
            )
        return (self.domain.geometric_dimension, tdim)

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        parametrization = self.parametrization
        assert self.point is not None
        return parametrization.jacobian(self.point)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        if self._is_partnered_by(other, JacobianInverse):
            return self._onto_the_tangent_space()
        return None


class MetricTensor(AbstractGeometricQuantity):
    """The metric a parametrization induces, also called the first fundamental form.

    The pull back of the Euclidean metric on the ambient coordinates by
    the map, ``g = J^T J``, so that lengths and angles measured in a
    cell's coordinates agree with the ambient ones. It is symmetric, and
    positive definite wherever the map is an immersion.

    Its determinant is the squared volume scaling, so ``sqrt(det g)`` is
    the factor an integral picks up on being pulled back, which is what
    :class:`JacobianDeterminant` gives.
    """

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        tdim = self._jacobian.value_shape[1]
        return (tdim, tdim)

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = _as_dense_matrix(self._jacobian.expand_geometry())
        return j.transpose().matmat(j)

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)


class TangentialProjector(AbstractGeometricQuantity):
    """The orthogonal projection of the ambient coordinates onto the tangent space.

    ``P = J g^-1 J^T``, which is what a map composed with its own
    pseudo-inverse comes to. It is symmetric and idempotent, and its rank
    is the topological dimension: it keeps the directions a cell can move
    in and discards the rest. For a hypersurface it is ``I - n (x) n``.

    Where the map is square it is the identity, since then the tangent
    space is the whole of the ambient coordinates.
    """

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        gdim, _ = self._jacobian.value_shape
        return (gdim, gdim)

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        gdim, tdim = self._jacobian.value_shape
        if gdim == tdim:
            return Identity(gdim)
        j = _as_dense_matrix(self._jacobian.expand_geometry())
        return j.matmat(j.compute_inverse())

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        if self._is_partnered_by(other, TangentialProjector):
            # A projection projected again is the same projection.
            return self
        return None


class UnitNormal(AbstractGeometricQuantity):
    """The unit normal to a domain of codimension one in its ambient coordinates.

    Where ``gdim == tdim + 1`` the directions orthogonal to every column of
    the Jacobian form a line, so the tangent space fixes the normal up to
    sign with no further information. The component form is the
    generalised cross product of the Jacobian's columns, ``n_i = (-1)^i``
    times the determinant of the Jacobian with row ``i`` removed,
    normalised.

    Which of the two directions is meant is a convention this does not
    fix. A consumer that needs the outward one of a cell supplies that
    orientation itself, since nothing here knows which side the cell is
    on.

    Hand it a facet's domain and it is the facet normal; hand it a surface
    mesh's own domain and it is that surface's normal.
    """

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        gdim, tdim = self._jacobian.value_shape
        if gdim != tdim + 1:
            raise ValueError(
                f"A normal is a single direction only where a domain has codimension "
                f"one, but this one has topological dimension {tdim} in {gdim} "
                f"coordinates."
            )
        if tdim == 0:
            raise NotImplementedError(
                "The normal to a domain of topological dimension zero is a sign, which "
                "is a convention rather than something to compute."
            )
        return (gdim,)

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        (gdim,) = self.value_shape
        j = _as_dense_matrix(self._jacobian.expand_geometry())
        tdim = gdim - 1

        components = []
        for i in range(gdim):
            minor = Matrix(
                [[j.component(r, c) for c in range(tdim)] for r in range(gdim) if r != i]
            )
            determinant = minor.compute_determinant()
            components.append(determinant if i % 2 == 0 else -determinant)

        norm = Sqrt(expression_sum(c * c for c in components))
        return Vector([c / norm for c in components])

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)


class JacobianDeterminant(AbstractGeometricQuantity):
    """The determinant of the Jacobian."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = _as_dense_matrix(self._jacobian.expand_geometry())
        return abs(j.compute_determinant())

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        raise ValueError("Cannot get a component of a scalar expression")


class JacobianInverse(AbstractGeometricQuantity):
    """The inverse of the Jacobian."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._jacobian.value_shape[::-1]

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = _as_dense_matrix(self._jacobian.expand_geometry())
        return j.compute_inverse()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        if self._is_partnered_by(other, Jacobian):
            # J+ J is the identity on the cell's own coordinates.
            return Identity(self.value_shape[0])
        return None


class JacobianTranspose(AbstractGeometricQuantity):
    """The transpose of the Jacobian."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._jacobian.value_shape[::-1]

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = _as_dense_matrix(self._jacobian.expand_geometry())
        return j.transpose()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        if self._is_partnered_by(other, JacobianInverseTranspose):
            # J+ J is the identity on the cell's own coordinates.
            return Identity(self.value_shape[0])
        return None


class JacobianInverseTranspose(AbstractGeometricQuantity):
    """The inverse transpose of the Jacobian."""

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return self._jacobian.value_shape

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        j = _as_dense_matrix(self._jacobian.expand_geometry())
        return j.compute_inverse().transpose()

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression."""
        return self.expand_geometry().component(*indices)

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """
        if self._is_partnered_by(other, JacobianTranspose):
            return self._onto_the_tangent_space()
        return None


def expand_geometry(
    expression: GraphNode,
) -> GraphNode:
    """Replace geometric quantities by expressions in a domain's parametrizations.

    A Jacobian becomes the derivative of the map it differentiates and a
    pushed forward point becomes the coordinates it lands on. What those
    expressions contain is the map's business: a finite element one gives
    a sum over its basis, a closed form one gives its own expression.
    """
    to_replace: dict[GraphNode, GraphNode] = {}

    for node in as_graph(expression):
        # The GraphNode check is for the type checker, which loses that fact
        # on the protocol check.
        if isinstance(node, GraphNode) and isinstance(node, ExpandableGeometry):
            to_replace[node] = node.expand_geometry()

    return replace(expression, to_replace)
