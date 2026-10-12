# Copyright (C) 2026 Jack S. Hale
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Interfaces: a domain two others meet along, presented from each side.

Two domains ``Omega_+`` and ``Omega_-`` meeting along ``Gamma`` of
codimension one, with an embedding of ``Gamma`` into the closure of each.
A field on either side reaches ``Gamma`` by pullback along its side's
embedding, and the jump ``iota_+^* u_+ - iota_-^* u_-`` typechecks
because both terms are fields on the one ``Gamma``. There is one
integration domain, one measure, and two maps named in the integrand.
INTEGRATION.md section 5 is the argument.

In the vocabulary here, a presentation of ``Gamma`` from one side is a
piece: a domain in a chart of that side's domain, carried by that
domain's own map. A facet of a cell mapped on by the cell's map is one,
and so is a face of a box. An interface is a domain and its
presentations, each named by a side, that are declared to cover one
set. That is the atlas a union of pieces is not: a union's pieces make
no claim to be related, where the sides of an interface are two charts
of one set with a transition map between them.

The transition map is never written as an object. The point of a side's
parameter region that a point ``X`` of the interface corresponds to is
``psi_s^-1(psi_Gamma(X))``, which is the composite of an image through
the interface's parametrization and a preimage through the side's, and
both of those exist already. The preimage is a terminal the surrounding
library evaluates: combinatorially, as a facet permutation, when the
interface is a facet of both sides' cells, and by point location when it
is not. No derivative of it is ever needed, since a restriction is the
pullback of a function and the measure is the interface's own.

The two sides of one mesh are parametrized by one map, structurally, so a
preimage through the ``-`` side's map would not say which side it was
through. A side's parametrization is therefore its presentation's map
named by the side, :class:`SideParametrization`, which is to a map what
a labelled domain is to a domain: identity and nothing geometric.

None of that needs the interface to be a facet of either side. The map
``iota_s`` written in charts is ``phi_s^-1 . T_s . psi_Gamma``, from the
interface's region into a cell's coordinates, and a facet inclusion is
the case where that is affine and an entity can write it down. In
general it is :class:`PreimageParametrization`: its value is a preimage,
the terminal again, and its Jacobian is the chain rule in closed form.
``T_s`` is a map the user states, the identity of sets unless said, which
is what makes a periodic boundary an interface whose other side is a
translation away. A presentation built from it is parametrized from the
interface's own region, so there is no chart transition to take. The
class lives here and not in `uflx.parametrizations` because it is built
from an image and a preimage, which `uflx.geometry` holds.

Which pieces of the interface lie in one cell of each side is a union
whose membership is external, as a mesh's cells are, and within a piece
the map is smooth. That is where a supermesh lives, and the language
says nothing about it.
"""

from __future__ import annotations

from itertools import product
from typing import Any

from uflx.domains import (
    RD,
    AbstractCoordinateDomain,
    AbstractParametrization,
    AbstractParametrizedDomain,
    CarriesADomain,
    EntityDomain,
    MadeOfPieces,
)
from uflx.expressions import AbstractExpression, MatrixProduct
from uflx.geometry import ImagePoint, PreimagePoint, _as_dense_matrix
from uflx.parametrizations import ComposedParametrization
from uflx.points import AbstractPoint, Point
from uflx.variables import AbstractVariable


class PreimageParametrization(AbstractParametrization):
    """A domain carried into a cell's coordinates, ``phi^-1 . T . psi``.

    The map of an interface's side written in charts, in general: a
    point of the domain's region is carried onto the domain by its own
    parametrization ``psi``, moved by a stated map ``T`` if there is
    one, and pulled back through a cell's map ``phi`` into that cell's
    coordinates. The value is a preimage, which the surrounding library
    answers by point location, or combinatorially where the point lands
    on a facet; the Jacobian is the chain rule, ``J_phi^-1 J_T J_psi``,
    in closed form, with the pseudo-inverse where the cell is a surface.

    Presenting the domain from a side is then
    ``composed_domain(composed_domain(region, this), phi)``: the domain
    in the cell's coordinates, mapped on by the cell's own map, which is
    the same shape a facet of a mesh has and is read the same way.
    Implements :class:`uflx.domains.CarriesADomain`.
    """

    def __init__(
        self,
        domain: AbstractParametrizedDomain,
        through: AbstractParametrization,
        by: AbstractParametrization | None = None,
    ):
        """Initialise.

        Args:
            domain: The domain carried, parametrized from one cell's
                coordinates
            through: The map of the cell it is carried into, out of that
                cell's coordinates
            by: The map applied to the domain's points first, from the
                coordinates they are named in to those the cell map lands
                in, or None for the identity of sets

        Raises:
            ValueError: If the domain is parametrized from several regions
                or from a region that is not a cell's, the cell map is
                not out of a cell, or the stated map does not join them
        """
        sources = {parametrization.source for parametrization in domain.parametrizations}
        if len(sources) != 1 or (isinstance(domain, MadeOfPieces) and domain.pieces != (domain,)):
            raise ValueError(
                f"{domain!r} is parametrized piece by piece. Carry one piece of it at a time."
            )
        (source,) = sources
        if not isinstance(source, EntityDomain):
            raise ValueError(
                f"{domain!r} is parametrized from {source!r}, which is not a cell's "
                f"coordinates, and a point is carried onto a domain only from a cell's."
            )
        if not isinstance(through.source, EntityDomain):
            raise ValueError(
                f"{through!r} is not the map of a cell, being out of {through.source!r}."
            )
        if by is None:
            if domain.geometric_dimension != through.target_dimension:
                raise ValueError(
                    f"{domain!r} is named in {domain.geometric_dimension} coordinates and "
                    f"{through!r} lands in {through.target_dimension}, so neither is a "
                    f"point of the other's without a map between them."
                )
        else:
            starts = (domain.coordinate_domain, RD(domain.geometric_dimension))
            lands = (through.target, RD(through.target_dimension))
            if by.source not in starts or by.target not in lands:
                raise ValueError(
                    f"{by!r} does not take the coordinates of {domain!r} to those "
                    f"{through!r} lands in."
                )
        self._domain = domain
        self._through = through
        self._by = by
        self._source = source

    @property
    def carries(self) -> AbstractParametrizedDomain:
        """The domain carried."""
        return self._domain

    @property
    def through(self) -> AbstractParametrization:
        """The map of the cell the domain is carried into."""
        return self._through

    @property
    def by(self) -> AbstractParametrization | None:
        """The map applied to the domain's points first, if one is stated."""
        return self._by

    @property
    def source(self) -> AbstractCoordinateDomain:
        """The region the domain is parametrized from."""
        return self._source

    @property
    def target(self) -> AbstractCoordinateDomain:
        """The cell's coordinates."""
        return self._through.source

    def _image(self, point: AbstractVariable) -> AbstractPoint:
        """Where the point sits once carried onto the domain and moved by the stated map."""
        image: AbstractPoint = ImagePoint(point, self._domain)
        if self._by is not None:
            moved = self._by.value(image)
            image = Point(
                [moved.component(j) for j in range(self._by.target_dimension)], self._by.target
            )
        return image

    def value(self, point: AbstractVariable) -> AbstractExpression:
        """The point of the cell that lands where the domain's point does."""
        return PreimagePoint(self._image(point), self._through)

    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """The chain rule, with the cell map inverted at the preimage.

        The inverse is written out rather than left as a quantity, so
        that expanding reaches the bottom in one go.
        """
        (psi,) = [p for p in self._domain.parametrizations if p.source == point.domain]
        preimage = self.value(point)
        assert isinstance(preimage, AbstractVariable)
        inverse = _as_dense_matrix(self._through.jacobian(preimage)).compute_inverse()
        along = psi.jacobian(point)
        if self._by is not None:
            # A product of two at a time, which is what a product's component reads.
            along = MatrixProduct([self._by.jacobian(ImagePoint(point, self._domain)), along])
        return MatrixProduct([inverse, along])

    def __repr__(self) -> str:
        """Representation."""
        return f"PreimageParametrization({self._domain!r}, {self._through!r}, {self._by!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, PreimageParametrization)
            and self._domain == other._domain
            and self._through == other._through
            and self._by == other._by
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.PreimageParametrization", self._domain, self._through, self._by))


class SideParametrization(AbstractParametrization):
    """A presentation's map, named by the side it presents.

    Delegates everything to the map it names. What it adds is identity:
    the two sides of one mesh are parametrized by one map, so a preimage
    through the ``-`` side's map has to say which side it is through,
    and this is how it says.
    """

    def __init__(self, parametrization: AbstractParametrization, side: Side):
        """Initialise.

        Args:
            parametrization: The map being named
            side: The side naming it
        """
        self._parametrization = parametrization
        self._side = side

    @property
    def parametrization(self) -> AbstractParametrization:
        """The map being named."""
        return self._parametrization

    @property
    def side(self) -> Side:
        """The side naming it."""
        return self._side

    @property
    def source(self) -> AbstractCoordinateDomain:
        """The map's source."""
        return self._parametrization.source

    @property
    def target(self) -> AbstractCoordinateDomain:
        """The map's target."""
        return self._parametrization.target

    def value(self, point: AbstractVariable) -> AbstractExpression:
        """The map's value."""
        return self._parametrization.value(point)

    def jacobian(self, point: AbstractVariable) -> AbstractExpression:
        """The map's Jacobian."""
        return self._parametrization.jacobian(point)

    @property
    def is_affine(self) -> bool:
        """Whether the map is affine."""
        return self._parametrization.is_affine

    @property
    def is_identity(self) -> bool:
        """Whether the map is the identity."""
        return self._parametrization.is_identity

    def __repr__(self) -> str:
        """Representation."""
        return f"SideParametrization({self._parametrization!r}, {self._side!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, SideParametrization)
            and self._parametrization == other._parametrization
            and self._side == other._side
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.SideParametrization", self._parametrization, self._side))


class Side(AbstractParametrizedDomain):
    """An interface as one of its sides presents it.

    The same set as the interface, parametrized by the side's own
    presentation. A field on the side's domain is restricted to the
    interface through this, and the normal pointing out of that domain is
    a quantity of this. Implements :class:`uflx.domains.PresentsAnInterface`.
    """

    def __init__(self, interface: Interface, label: Any):
        """Initialise.

        Args:
            interface: The interface this is a side of
            label: Which side, as the interface names it
        """
        if label not in interface.labels:
            raise ValueError(f"{interface!r} has no side {label!r}.")
        self._interface = interface
        self._label = label

    @property
    def interface(self) -> Interface:
        """The interface this is a side of."""
        return self._interface

    @property
    def label(self) -> Any:
        """Which side."""
        return self._label

    @property
    def presentation(self) -> AbstractParametrizedDomain:
        """The interface as this side's domain presents it."""
        return self._interface.presentation(self._label)

    @property
    def geometric_dimension(self) -> int:
        """The coordinates the presentation is named in."""
        return self.presentation.geometric_dimension

    @property
    def topological_dimension(self) -> int | None:
        """The presentation's topological dimension."""
        return self.presentation.topological_dimension

    @property
    def parametrizations(self) -> tuple[AbstractParametrization, ...]:
        """The presentation's maps, each named by this side."""
        return tuple(
            SideParametrization(parametrization, self)
            for parametrization in self.presentation.parametrizations
        )

    def chart_point(self, variable: AbstractVariable) -> AbstractVariable:
        """Where in this side's parameter region a point of the interface is.

        The variable itself where this is the interface's primary side,
        the one whose presentation the domain is, since then the regions
        are one, and where the presentation carries the domain, since
        then the transition is inside its map. Otherwise
        the chart transition, ``psi_s^-1(psi_Gamma(X))``: the point
        carried onto the interface by its parametrization and pulled back
        through this side's, which is a terminal the surrounding library
        evaluates.

        Args:
            variable: A point of a parameter region of the interface. One
                of anywhere else is given back untouched, being none of
                this side's business.

        Returns:
            The point of this side's parameter region that lands on it

        Raises:
            ValueError: If this side has several parametrizations, so no
                one of them to pull back through
        """
        domain = self._interface.domain
        if variable.domain not in {p.source for p in domain.parametrizations}:
            return variable
        if self._label == self._interface.primary or _carried_by(self.presentation) == domain:
            return variable
        parametrizations = self.parametrizations
        if len(parametrizations) != 1:
            raise ValueError(
                f"{self!r} has {len(parametrizations)} parametrizations, so a point of the "
                f"interface names no one of them to pull back through. Split the "
                f"integral by piece first."
            )
        (parametrization,) = parametrizations
        return PreimagePoint(ImagePoint(variable, domain), parametrization)

    def __repr__(self) -> str:
        """Representation."""
        return f"Side({self._interface!r}, {self._label!r})"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, Side)
            and self._interface == other._interface
            and self._label == other._label
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.Side", self._interface, self._label))


class Interface(AbstractParametrizedDomain):
    """A domain two others meet along, with a presentation of it from each side.

    The domain is what is integrated over: its parametrization gives the
    measure, the metric and the unsigned normal direction. Where the
    interface is a facet of both sides' cells, it is one side's
    presentation, and the other induces the same density, the two
    differing by an isometry of the reference facet.

    The sides are ordered, and the first is first: ``jump`` subtracts the
    second from it. The primary side is the first whose presentation is
    the domain itself, and is read at the integral's own point; every
    other side is read through the chart transition, a second side with
    the same presentation included, since that is the interface of a
    mesh with itself and its two sides are two cells.

    Implements :class:`uflx.domains.MadeOfPieces`: an interface of unions
    of pieces, which the facets of a mesh are, is integrated over as the
    sum over the pairs of them.
    """

    def __init__(
        self,
        domain: AbstractParametrizedDomain,
        sides: tuple[tuple[Any, AbstractParametrizedDomain], ...],
    ):
        """Initialise.

        Args:
            domain: The domain integrated over
            sides: A label and a presentation of the domain per side, in
                order

        A presentation need not be named in the domain's coordinates: a
        side's map may take the domain somewhere else first, as a
        periodic boundary's does, so what the sides share is the domain
        and not an ambient space.

        Raises:
            ValueError: If there are fewer than two sides, two share a
                label, or a presentation is of a different topological
                dimension
        """
        if len(sides) < 2:
            raise ValueError("An interface is where two domains meet, so it has two sides.")
        labels = [label for label, _ in sides]
        if len(set(labels)) != len(labels):
            raise ValueError(f"The sides of an interface are told apart by label, not {labels}.")
        for label, presentation in sides:
            if presentation.topological_dimension != domain.topological_dimension:
                raise ValueError(
                    f"Side {label!r} has topological dimension "
                    f"{presentation.topological_dimension} where {domain!r} has "
                    f"{domain.topological_dimension}."
                )
        self._domain = domain
        self._sides = sides

    @property
    def domain(self) -> AbstractParametrizedDomain:
        """The domain integrated over."""
        return self._domain

    @property
    def labels(self) -> tuple[Any, ...]:
        """The sides' labels, in order."""
        return tuple(label for label, _ in self._sides)

    @property
    def primary(self) -> Any:
        """The side the domain is the presentation of, or None if it is its own.

        The first side whose presentation is the domain. One whose domain
        is none of its sides' is the general interface, read on from every
        side through the chart transition.
        """
        for label, presentation in self._sides:
            if presentation == self._domain:
                return label
        return None

    def presentation(self, label: Any) -> AbstractParametrizedDomain:
        """The domain as one side presents it.

        Args:
            label: Which side

        Returns:
            That side's presentation
        """
        for side, presentation in self._sides:
            if side == label:
                return presentation
        raise ValueError(f"{self!r} has no side {label!r}.")

    def side(self, label: Any) -> Side:
        """One side of this interface.

        Args:
            label: Which side

        Returns:
            That side, as a domain
        """
        return Side(self, label)

    @property
    def sides(self) -> tuple[Side, ...]:
        """The sides, in order."""
        return tuple(Side(self, label) for label in self.labels)

    @property
    def geometric_dimension(self) -> int:
        """The coordinates the domain is named in."""
        return self._domain.geometric_dimension

    @property
    def topological_dimension(self) -> int | None:
        """The domain's topological dimension."""
        return self._domain.topological_dimension

    @property
    def parametrizations(self) -> tuple[AbstractParametrization, ...]:
        """The domain's own maps, which are what the measure is read off."""
        return self._domain.parametrizations

    @property
    def pieces(self) -> tuple[AbstractParametrizedDomain, ...]:
        """One interface per combination of the sides' pieces.

        The primary side follows the domain's piece, and so does a piece
        that carries it, pairing with no other. Combinations whose pieces
        start from different regions are left out, a triangle's facet
        being no interface with a tetrahedron's: the sides of one
        interface are charted from one region. An interface with no
        pieces to combine is its own one piece.
        """
        if not isinstance(self._domain, MadeOfPieces) and not any(
            isinstance(presentation, MadeOfPieces) for _, presentation in self._sides
        ):
            return (self,)
        pieces: list[AbstractParametrizedDomain] = []
        primary = self.primary
        for gamma in _pieces_of(self._domain):
            choices = [
                (gamma,)
                if label == primary
                else [p for p in _pieces_of(presentation) if _carried_by(p) in (None, gamma)]
                for label, presentation in self._sides
            ]
            for combination in product(*choices):
                regions = {
                    parametrization.source
                    for piece in (gamma, *combination)
                    for parametrization in piece.parametrizations
                }
                if len(regions) != 1:
                    continue
                pieces.append(Interface(gamma, tuple(zip(self.labels, combination, strict=True))))
        return tuple(pieces)

    def __repr__(self) -> str:
        """Representation."""
        sides = ", ".join(f"({label!r}, {presentation!r})" for label, presentation in self._sides)
        return f"Interface({self._domain!r}, ({sides}))"

    def __eq__(self, other) -> bool:
        """Check for equality."""
        return (
            isinstance(other, Interface)
            and self._domain == other._domain
            and self._sides == other._sides
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(("uflx.Interface", self._domain, self._sides))


def _carried_by(presentation: AbstractParametrizedDomain) -> AbstractParametrizedDomain | None:
    """The domain a presentation's maps carry, where they all carry one.

    The map of a presentation built from a :class:`PreimageParametrization`
    is that followed by the cell's own, so the carrier is the inner map
    of a composition.
    """
    carried: set[AbstractParametrizedDomain | None] = set()
    for parametrization in presentation.parametrizations:
        inner = parametrization
        while isinstance(inner, ComposedParametrization):
            inner = inner.inner
        carried.add(inner.carries if isinstance(inner, CarriesADomain) else None)
    if len(carried) == 1 and None not in carried:
        return carried.pop()
    return None


def _pieces_of(domain: AbstractParametrizedDomain) -> tuple[AbstractParametrizedDomain, ...]:
    """The pieces of a domain, or the domain itself when it is one piece."""
    if isinstance(domain, MadeOfPieces):
        return domain.pieces
    return (domain,)


def interface(
    domain: AbstractParametrizedDomain,
    sides: tuple[tuple[Any, AbstractParametrizedDomain], ...],
) -> Interface:
    """Create the interface a domain is, presented from each of its sides.

    Args:
        domain: The domain integrated over
        sides: A label and a presentation of the domain per side, in order

    Returns:
        That interface
    """
    return Interface(domain, sides)
