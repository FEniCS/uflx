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
"""

from __future__ import annotations

from itertools import product
from typing import Any

from uflx.domains import (
    AbstractCellularDomain,
    AbstractCoordinateDomain,
    AbstractParametrization,
    AbstractParametrizedDomain,
    MadeOfPieces,
)
from uflx.expressions import AbstractExpression
from uflx.geometry import ImagePoint, PreimagePoint
from uflx.parametrizations import sub_entities_of
from uflx.variables import AbstractVariable


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
        are one. Otherwise
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
        if self._label == self._interface.primary:
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

        Raises:
            ValueError: If there are fewer than two sides, two share a
                label, or a presentation is of a different dimension
        """
        if len(sides) < 2:
            raise ValueError("An interface is where two domains meet, so it has two sides.")
        labels = [label for label, _ in sides]
        if len(set(labels)) != len(labels):
            raise ValueError(f"The sides of an interface are told apart by label, not {labels}.")
        for label, presentation in sides:
            if presentation.geometric_dimension != domain.geometric_dimension:
                raise ValueError(
                    f"Side {label!r} is named in {presentation.geometric_dimension} "
                    f"coordinates where {domain!r} is named in {domain.geometric_dimension}."
                )
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

        The primary side follows the domain's piece. Combinations whose pieces start from different
        regions are left out, a triangle's facet being no interface with
        a tetrahedron's: the sides of one interface are charted from one
        region. An interface with no pieces to combine is its own one
        piece.
        """
        if not isinstance(self._domain, MadeOfPieces) and not any(
            isinstance(presentation, MadeOfPieces) for _, presentation in self._sides
        ):
            return (self,)
        pieces: list[AbstractParametrizedDomain] = []
        primary = self.primary
        for gamma in _pieces_of(self._domain):
            choices = [
                (gamma,) if label == primary else _pieces_of(presentation)
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


def interface_of(
    plus: AbstractCellularDomain, minus: AbstractCellularDomain, label: Any = None
) -> Interface:
    """Create the interface along which the facets of two domains of cells meet.

    UFL's ``dS``, as a domain: the one measure over this is the interior
    facet integral, with the fields read on it through the sides ``"+"``
    and ``"-"``. Passing one domain twice is the interface of a mesh with
    itself. Which facets two cells actually share is external, as the
    cells are; this says only how each pair of them is parametrized.

    Args:
        plus: The domain whose facets are the ``"+"`` side
        minus: The domain whose facets are the ``"-"`` side
        label: Anything hashable the consumer selects these facets by, or
            None for no selection

    Returns:
        The interface, over the ``"+"`` side's facets
    """
    facets = sub_entities_of(plus, 1, label)
    return Interface(facets, (("+", facets), ("-", sub_entities_of(minus, 1, label))))
