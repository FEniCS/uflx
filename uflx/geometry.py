"""Geometry."""

from __future__ import annotations

from abc import abstractmethod
from typing import Any, Protocol, Self, runtime_checkable

from uflx.algorithms import replace
from uflx.domains import (
    AbstractCoordinateDomain,
    AbstractParametrization,
    AbstractParametrizedDomain,
    EntityDomain,
    PresentsAnInterface,
    StatesAnOutwardNormal,
    unlabelled,
)
from uflx.expressions import AbstractExpression, RealScalar, Sqrt, expression_sum
from uflx.graphs import GraphNode, as_graph
from uflx.parametrizations import included_sub_entity, sub_entity_in_cell
from uflx.points import AbstractPoint, Point
from uflx.tensors import Identity, Matrix, Vector
from uflx.variables import AbstractVariable


@runtime_checkable
class ExpandableGeometry(Protocol):
    """Geometry that can be expanded into components."""

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""


def _parametrizations_out_of(
    domain: AbstractParametrizedDomain, source: AbstractCoordinateDomain
) -> tuple[AbstractParametrization, ...]:
    """The domain's parametrizations that start from the given region.

    A domain holds its parametrizations rather than keying them by their region,
    since several may share one. Needing exactly one is what evaluating a
    quantity at a point asks for, so it is asked here.
    """
    return tuple(
        parametrization
        for parametrization in domain.parametrizations
        if parametrization.source == source
    )


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


class ImagePoint(AbstractPoint):
    """Where a parametrization sends a point of the region it starts from.

    Applying a map to a point, which is an image and not a pushforward:
    a pushforward and a pullback act on fields, and `uflx.maps` is where
    those live. The same thing
    :class:`SpatialCoordinate` gives as an expression, but being a point
    it can be the variable a basis function is evaluated at. An element
    defined on the physical cell rather than on a reference one needs
    that, having no map to compose from.

    It takes the domain rather than the map because the point it starts
    from already lies in a parameter region of the domain, a cell's
    coordinate domain or a region of R^d, and so names the map that
    carries it. Contrast :class:`PreimagePoint`.

    Where it lands is the domain's coordinate domain, which is the
    ambient one for a mesh and the cell's own for a sub-entity of a
    reference cell. That is what lets a cell's basis be read at a point
    of one of its facets.
    """

    def __init__(self, point: AbstractVariable, domain: AbstractParametrizedDomain):
        """Initialise.

        Args:
            point: The point in a parameter region of the domain
            domain: The domain to push the point forward onto

        Raises:
            ValueError: If the domain has no parametrization out of the
                point's region
        """
        source = point.domain
        if not isinstance(source, AbstractCoordinateDomain) or not _parametrizations_out_of(
            domain, source
        ):
            raise ValueError(
                f"A point is pushed forward from a cell's coordinates, or from another "
                f"region {domain!r} is parametrized from, not from {source!r}."
            )
        self._point = point
        self._parametrized_domain = domain

    @property
    def entity_point(self) -> AbstractVariable:
        """The point in the entity's coordinates."""
        return self._point

    @property
    def parametrized_domain(self) -> AbstractParametrizedDomain:
        """The domain this point is pushed forward onto."""
        return self._parametrized_domain

    @property
    def parametrization(self) -> AbstractParametrization:
        """The map this point is pushed forward through.

        The point lies in a parameter region of the domain, so it names
        the map out of that region.

        Raises:
            ValueError: If several of the domain's maps start from that
                region, so the point names no one of them
        """
        source = self._point.domain
        assert isinstance(source, AbstractCoordinateDomain)
        parametrizations = _parametrizations_out_of(self._parametrized_domain, source)
        if len(parametrizations) != 1:
            raise ValueError(
                f"{self._parametrized_domain!r} has {len(parametrizations)} parametrizations "
                f"out of {source!r}, so a point there names no one map. Split the "
                f"integral by piece first."
            )
        (parametrization,) = parametrizations
        return parametrization

    @property
    def domain(self) -> AbstractCoordinateDomain:
        """The domain.

        An image lies in the coordinates the map lands in, which the
        domain being mapped onto names.
        """
        return self._parametrized_domain.coordinate_domain

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
            isinstance(other, ImagePoint)
            and self._point == other._point
            and self._parametrized_domain == other._parametrized_domain
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.ImagePoint", self._point, self._parametrized_domain))

    @property
    def index(self) -> int | str:
        """The point's index in the set of points.

        Raises:
            NotImplementedError: If the point mapped is a variable
                standing for a point rather than one of a set
        """
        if not isinstance(self._point, AbstractPoint):
            raise NotImplementedError(f"{self._point!r} is not a point of a set of them.")
        return self._point.index


class PreimagePoint(AbstractPoint):
    """The point of a parametrization's region that lands on a given ambient point.

    The inverse image of a point, and not a pullback: a pullback acts on
    fields, and `uflx.maps` is where those live.

    A terminal, deliberately: it has no expansion and will not get one.
    A parametrization offers a value and a derivative and no inverse, the
    map back being an algorithm the surrounding library provides. So this
    names a point that library computes, the way a coordinate dof names a
    number a mesh holds.

    It takes the map rather than the domain because a point of the ambient
    coordinates does not say which cell it should land in; that is part of
    the question being asked. Contrast :class:`ImagePoint`.

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
            isinstance(other, PreimagePoint)
            and self._point == other._point
            and self._parametrization == other._parametrization
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.PreimagePoint", self._point, self._parametrization))

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

    def __init__(
        self, domain: AbstractParametrizedDomain, variable: AbstractVariable | None = None
    ):
        """Initialise.

        Args:
            domain: The domain whose geometry is being differentiated
            variable: The point to differentiate it at, if that is known yet
        """
        self._domain = domain
        self._variable = variable

    @property
    def domain(self) -> AbstractParametrizedDomain:
        """The domain whose geometry this quantity differentiates."""
        return self._domain

    @property
    def variable(self) -> AbstractVariable | None:
        """Where it is differentiated, if it has been told yet."""
        return self._variable

    @property
    def reads_on(self) -> AbstractParametrizedDomain | None:
        """The domain this quantity is read on, until it is told its point.

        One already told was put there by whatever knew the point -- a
        pull back, or a restriction -- so it says nothing about the
        domain the expression around it is on.
        """
        if self._variable is not None:
            return None
        return self._domain

    @property
    def parametrization(self) -> AbstractParametrization:
        """The map this quantity differentiates, being the one its point lies in.

        A point lies in the parameter region of one of the domain's
        parametrizations, so it names that map. Until a point arrives this
        quantity is generic over the domain's parametrizations, which is
        what lets it be built during a pull back.
        """
        if self.variable is None:
            raise ValueError(
                "This quantity has not been told where it is evaluated, so the map "
                "it differentiates is not known."
            )
        source = self.variable.domain
        if not isinstance(source, AbstractCoordinateDomain):
            raise ValueError(
                f"This quantity is evaluated at a point of a coordinate domain, not of {source!r}."
            )
        parametrizations = _parametrizations_out_of(self.domain, source)
        if len(parametrizations) == 0:
            raise ValueError(f"{self.domain!r} has no parametrization out of {source!r}.")
        if len(parametrizations) > 1:
            raise ValueError(
                f"{self.domain!r} has {len(parametrizations)} parametrizations out "
                f"of {source!r}, so a point there names no one map. Split the integral "
                f"by piece first."
            )
        (parametrization,) = parametrizations
        return parametrization

    @property
    def _jacobian(self) -> Jacobian:
        """The Jacobian this quantity is built from."""
        return Jacobian(self.domain, self.variable)

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._domain, self._variable

    def __repr__(self) -> str:
        """Representation."""
        return f"{self.__class__.__name__}({self._domain!r}, {self._variable!r})"

    @abstractmethod
    def expand_geometry(self) -> AbstractExpression:
        """Expand this quantity into an expression in the domain's parametrization."""

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression.

        Every quantity is read by expanding it, so this is asked of the
        expansion rather than answered by each quantity in turn. A scalar
        one has no component to get.

        Raises:
            ValueError: If this quantity is a scalar
        """
        if self.value_shape == ():
            raise ValueError("Cannot get a component of a scalar expression")
        return self.expand_geometry().component(*indices)

    def _is_partnered_by(self, other: GraphNode, partner: type[AbstractGeometricQuantity]) -> bool:
        """Whether `other` is a quantity of the given kind on the same geometry.

        Args:
            other: The quantity this one is multiplied by
            partner: The class it must be for the product to simplify
        """
        return (
            isinstance(other, partner)
            and self._domain == other.domain
            and self._variable == other.variable
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
        return TangentialProjector(self._domain, self._variable)

    def _rebuilt(
        self, domain: AbstractParametrizedDomain, variable: AbstractVariable | None
    ) -> Self:
        """This quantity on a domain, read at a variable.

        The one place a quantity is rebuilt, so a subclass carrying more
        than a domain and a point overrides this alone and neither of the
        two callers has to know what else it carries.
        """
        return self.__class__(domain, variable)

    def on_domain(self, domain: AbstractParametrizedDomain) -> Self:
        """This quantity on another domain, not yet told where it is read.

        What splitting an integral over a union of pieces does to the
        quantities in it: the normal of a domain of facets is the normal
        of the piece the integral is split onto, read at that piece's own
        point rather than at the union's.

        Args:
            domain: The domain to read this quantity off instead

        Returns:
            The same kind of quantity there
        """
        return self._rebuilt(domain, None)

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Self:
        """Evaluate this quantity at the given variable.

        The variable stands for a point of one of this domain's parameter
        regions, so one in any other coordinates is not this quantity's to
        take. A quantity of a side of an interface is read at the point
        of the side's own chart that the interface's point corresponds
        to.
        """
        source = variable.domain
        if not isinstance(source, AbstractCoordinateDomain):
            return self
        if len(_parametrizations_out_of(self._domain, source)) == 0:
            return self
        if isinstance(self._domain, PresentsAnInterface):
            variable = self._domain.chart_point(variable)
        return self._rebuilt(self._domain, variable)


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
        assert self._variable is not None
        return parametrization.value(self._variable)

    def __getitem__(self, component: int) -> SingleSpatialCoordinate:
        """Get one coordinate."""
        return self.component(component)

    def component(self, *indices: int) -> SingleSpatialCoordinate:
        """Get a component of the expression."""
        (i,) = indices
        (dimension,) = self.value_shape
        if i < 0 or i >= dimension:
            raise IndexError("coordinate index out of range")
        return SingleSpatialCoordinate(self._domain, i, self._variable)


class SingleSpatialCoordinate(AbstractGeometricQuantity):
    """One of the ambient coordinates a domain's parametrization lands on.

    Kept as a node of its own rather than expanded on sight, so that
    ``sin(x[0])`` stays readable until geometry is expanded.
    """

    def __init__(
        self,
        domain: AbstractParametrizedDomain,
        component: int,
        variable: AbstractVariable | None = None,
    ):
        """Initialise.

        Args:
            domain: The domain whose coordinates these are
            component: Which coordinate, in ``range(geometric_dimension)``
            variable: The point to evaluate it at, if that is known yet
        """
        super().__init__(domain, variable)
        self._component = component

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._domain, self._component, self._variable

    def _rebuilt(
        self, domain: AbstractParametrizedDomain, variable: AbstractVariable | None
    ) -> Self:
        """This coordinate, rebuilt, carrying which coordinate it is."""
        return self.__class__(domain, self._component, variable)

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        return (
            SpatialCoordinate(self._domain, self._variable)
            .expand_geometry()
            .component(self._component)
        )

    def __repr__(self) -> str:
        """Representation."""
        return f"SingleSpatialCoordinate({self._domain!r}, {self._component}, {self._variable!r})"


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
        assert self.variable is not None
        return parametrization.jacobian(self.variable)

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
    :class:`VolumeElement` gives.
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


class VolumeElement(AbstractGeometricQuantity):
    """The density a parametrization induces, ``sqrt(det g)``.

    The factor an integral picks up on being pulled back to a cell's
    coordinates: the volume the map gives a unit volume of the cell. On a
    manifold it is ``sqrt(det(J^T J))``, and where the map is square it
    reduces to ``abs(det J)``, which is cheaper.

    Never negative, being a density rather than a volume form. Integrating
    a function needs only this. The sign a density discards is the
    orientation the map gives the cell, which :class:`JacobianDeterminant`
    keeps where the map is square.
    """

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return ()

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        gdim, tdim = self._jacobian.value_shape
        determinant = _as_dense_matrix(self._jacobian.expand_geometry()).compute_determinant()
        if gdim == tdim:
            return abs(determinant)
        # A non-square map's pseudo-determinant is a square root already.
        return determinant


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


def _normalised(components: list[AbstractExpression]) -> Vector:
    """A vector scaled to unit length."""
    norm = Sqrt(expression_sum(c * c for c in components))
    return Vector([c / norm for c in components])


class OutwardNormal(AbstractGeometricQuantity):
    """The unit normal pointing out of a domain along a facet of it.

    The outward conormal of a manifold with boundary, which needs no
    orientation: outwardness is about which side of the facet the
    interior lies on, and the map carries interior to interior. So the
    reference facet's outward normal, which the cell supplies, transported
    by the cell map's inverse transpose, is outward even under a
    reflection, where the cross product of the mapped facet's own
    Jacobian is not. INTEGRATION.md section 4 is the argument. On a
    manifold the pseudo-inverse gives the conormal within the surface.

    Hand it a facet of a cell mapped on by the cell's own map, or a side
    of an interface presented by one, and it is the normal out of that
    cell. A domain that states its own outward normal, as a face of a box
    does, is taken at its word. :class:`UnitNormal` is the direction a
    domain of codimension one has on its own, with no cell to be out of,
    and the two agree up to sign.
    """

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        return (self._domain.geometric_dimension,)

    @property
    def _presentation(self) -> AbstractParametrizedDomain:
        """The facet this points out of, which a side presents."""
        domain = self._domain
        if isinstance(domain, PresentsAnInterface):
            domain = domain.presentation
        presentation = unlabelled(domain)
        assert isinstance(presentation, AbstractParametrizedDomain)
        return presentation

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry.

        Raises:
            ValueError: If the domain is not a facet of a cell, or has not
                been told where it is read
        """
        presentation = self._presentation
        if isinstance(presentation, StatesAnOutwardNormal):
            return Vector([RealScalar(c) for c in presentation.outward_normal])
        inner, cell, cell_map = sub_entity_in_cell(presentation)
        sub_entity = included_sub_entity(presentation)
        if sub_entity is None or sub_entity[0] != cell.topological_dimension - 1:
            raise ValueError(
                f"{presentation!r} is not a facet of {cell}, so there is no out of the "
                f"cell across it. The normal of a domain presented some other way is "
                f"UnitNormal of it, with the sign a side names."
            )
        _, index = sub_entity
        if self._variable is None:
            raise ValueError("This normal has not been told where it is read.")
        reference = [RealScalar(c) for c in cell.facet_outward_normal(index)]
        # Expanding has to reach the bottom in one go, so the included point
        # is written out, for a map whose Jacobian reads its coordinates, and
        # the cell map's Jacobian is written out rather than left as a quantity.
        included = ImagePoint(self._variable, inner).expand_geometry()
        assert isinstance(included, AbstractVariable)
        jacobian = _as_dense_matrix(cell_map.jacobian(included))
        inverse_transpose = jacobian.compute_inverse().transpose()
        rows, cols = inverse_transpose.value_shape
        transported = [
            expression_sum(inverse_transpose.component(i, j) * reference[j] for j in range(cols))
            for i in range(rows)
        ]
        return _normalised(transported)


class JacobianDeterminant(AbstractGeometricQuantity):
    """The signed determinant of a square Jacobian.

    Defined only where the map is square, a non-square matrix having no
    determinant. The factor an integral picks up on being pulled back is
    :class:`VolumeElement`, which is defined either way and is never
    negative.

    The sign is the orientation the map gives the cell. Nothing consumes
    it yet, and it is named because it is what a density throws away.
    """

    def _check_square(self) -> None:
        """Check that the map has a determinant at all.

        Raises:
            ValueError: If the map is not square
        """
        gdim, tdim = self._jacobian.value_shape
        if gdim != tdim:
            raise ValueError(
                f"A map from {tdim} coordinates into {gdim} has no determinant. The "
                f"factor an integral picks up on being pulled back is VolumeElement."
            )

    @property
    def value_shape(self) -> tuple[int, ...]:
        """The value shape of the expression."""
        self._check_square()
        return ()

    def expand_geometry(self) -> AbstractExpression:
        """Expand geometry."""
        self._check_square()
        j = _as_dense_matrix(self._jacobian.expand_geometry())
        return j.compute_determinant()


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
