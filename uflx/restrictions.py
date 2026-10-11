# Copyright (C) 2026 Jack S. Hale
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT
"""Restriction of a field to a domain it is not defined on.

A field on a domain reaches a piece of that domain's boundary the way a
field reaches anywhere else: through a map. If ``iota`` parametrizations a facet of
a cell, the restriction of ``u`` to that facet is the pullback
``iota^* u``, the field read at the point the parametrization gives. Written out,
``(iota^* u)(X) = u(iota(X))``.

That is why an exterior facet integral needs no new measure and no new
kind of function space. The integral is over the facet, so its dummy
variable is a point of the reference facet; the field is still the
cell's, read at the included point. A gradient inside a restriction is
therefore still the cell's gradient, carried by the cell's inverse
Jacobian, rather than the facet's tangential one -- which is the whole
reason the restriction is written in the integrand and not folded into
the measure.

Only the value is left alone. A restriction here is the pullback of a
function and not of a form, so it does not carry a value map; what the
restriction of a 1-form field should be is `uflx.maps`'s question.

One map, and so never nested. Pulling back twice is pulling back along
the composition, ``iota_2^*(iota_1^* u) = (iota_1 . iota_2)^* u``, so
iterating buys no expressiveness and costs identity: the edges of a
tetrahedron's faces are twelve parametrizations for six ridges, and none of them
says which ridge it is, where ``sub_entities_of(mesh, 2)`` gives each one
once in the cell's own numbering. So a field is restricted to a
sub-entity of its own cells, and a sub-entity of one of those is reached
by asking the cell for the codimension rather than by restricting again.

What that leaves out is not a nested restriction. A field living on a
facet -- a multiplier on one -- restricted to that facet's own boundary
is a single map out of the domain it lives on, like any other, since the
facet is where the field is.
"""

from __future__ import annotations

from typing import Any

from uflx.algorithms import replace
from uflx.domains import (
    AbstractCellularDomain,
    AbstractParametrization,
    AbstractParametrizedDomain,
    ComposesADomain,
    EntityDomain,
    UnionDomain,
    unlabelled,
)
from uflx.entities import AbstractEntity
from uflx.expressions import AbstractExpression
from uflx.functions import extract_domain
from uflx.geometry import ImagePoint
from uflx.graphs import GraphNode, as_graph
from uflx.variables import AbstractVariable, ReadAtAVariable


def evaluate_at(
    expression: GraphNode, variable: AbstractVariable, skip: frozenset[GraphNode] = frozenset()
) -> GraphNode:
    """Tell the parts of an expression that do not know it where they are read.

    A function or a geometric quantity can be built before there is a
    point to read it at, and is told later. This is that telling, done
    once for an integral's dummy variable and once for the point a
    restriction's parametrization gives.

    Args:
        expression: The expression to evaluate
        variable: The point its parts are read at
        skip: Nodes to leave alone, being read somewhere else

    Returns:
        The expression with its parts evaluated at the variable
    """
    replacements: dict[GraphNode, GraphNode] = {}
    for node in as_graph(expression):
        if node in skip:
            continue
        if isinstance(node, ReadAtAVariable) and node.variable is None:
            # A node whose coordinates these are not gives itself back.
            read = node.reconstruct_with_variable(variable)
            if read is not node:
                replacements[node] = read
    if len(replacements) == 0:
        return expression
    return replace(expression, replacements)


class Restriction(AbstractExpression):
    """A field read at the point a domain's parametrization gives it, ``iota^* u``.

    The domain restricted to is a piece of a sub-entity domain, or the
    union of such pieces before an integral over it has been split. Until
    it is one piece there is no one inclusion to read the field through,
    so a union answers nothing about its parametrization and the pull back refuses.
    """

    def __init__(
        self,
        field: AbstractExpression,
        domain: AbstractParametrizedDomain,
        variable: AbstractVariable | None = None,
    ):
        """Initialise.

        Args:
            field: The field being restricted, on the domain the pieces
                are pieces of
            domain: The domain it is restricted to
            variable: Where the restriction is read, if that is known yet

        Raises:
            ValueError: If the domain is not made of sub-entities of the
                field's own domain
        """
        if isinstance(field, Restriction):
            # Now that a restriction says which domain it is read on, this
            # would typecheck: the composition is a map like any other. It
            # is refused because one map is the rule and not because the
            # check cannot see it, which is the module docstring's argument.
            raise ValueError(
                f"{field!r} is already read through a map, and a field is read through "
                f"one. Restricting along the composition is what this would be: take "
                f"the sub-entities of the cell, with sub_entities_of(domain, codim)."
            )
        field_domain = extract_domain(field)
        if not isinstance(field_domain, AbstractCellularDomain):
            raise ValueError(f"{field_domain!r} has no cells to include a sub-entity into.")
        for piece in _pieces_of(domain):
            _, cell, cell_map = _inclusion_of(piece)
            if cell not in field_domain.cell_types or cell_map != field_domain.parametrization(
                cell
            ):
                raise ValueError(
                    f"{piece!r} is not a sub-entity of a cell of {field_domain!r}, mapped "
                    f"on by that cell's own map. A field is read through one map, so a "
                    f"sub-entity of a sub-entity is reached by taking it of the cell: "
                    f"sub_entities_of(domain, codim), at the codimension the cell counts "
                    f"in."
                )
        self._field = field
        self._domain = domain
        self._variable = variable

    @property
    def field(self) -> AbstractExpression:
        """The field being restricted."""
        return self._field

    @property
    def domain(self) -> AbstractParametrizedDomain:
        """The domain it is restricted to."""
        return self._domain

    @property
    def variable(self) -> AbstractVariable | None:
        """Where the restriction is read, if it has been told yet."""
        return self._variable

    @property
    def reads_on(self) -> AbstractParametrizedDomain:
        """The domain this restriction is read on, which is what it restricts to.

        Always, where a function or a quantity stops saying once it has
        been positioned: a restriction is the map, so the field under it
        being on another domain is the whole point rather than a
        disagreement.
        """
        return self._domain

    @property
    def value_shape(self) -> tuple[int, ...]:
        """A pullback of a function's values leaves them as they are."""
        return self._field.value_shape

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""
        return {self._field}

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""
        return self._field, self._domain, self._variable

    def component(self, *indices: int) -> AbstractExpression:
        """Get a component of the expression.

        Reading a component and reading through a parametrization commute, the
        values being left alone.
        """
        return Restriction(self._field.component(*indices), self._domain, self._variable)

    def reconstruct_with_variable(self, variable: AbstractVariable) -> Restriction:
        """Read this restriction at the given variable."""
        return Restriction(self._field, self._domain, variable)

    def on_piece(self, piece: AbstractParametrizedDomain) -> Restriction:
        """Get this restriction onto one piece of the domain it restricts to.

        What the split of an integral over a union does to it: the
        inclusion a field is read through is a piece's and not a union's.

        Args:
            piece: One piece of this restriction's domain

        Returns:
            The same field restricted to that piece
        """
        return Restriction(self._field, piece)

    @property
    def image_point(self) -> ImagePoint:
        """The point of the cell this restriction reads its field at.

        The restriction's own point lies in the reference sub-entity's
        coordinates, and the inclusion carries it into the cell's, where
        the field's basis is defined.

        Raises:
            ValueError: If this restriction has not been told where it is
                read, or is over a union rather than one piece
        """
        if self._variable is None:
            raise ValueError("This restriction has not been told where it is read.")
        if isinstance(self._domain, UnionDomain):
            raise ValueError(
                f"{self._domain!r} has a parametrization per piece. Split the "
                f"integral by piece first."
            )
        inclusion, _, _ = _inclusion_of(self._domain)
        (sub_entity,) = inclusion.cell_types
        return ImagePoint(self._variable.to_entity_coordinates(sub_entity), inclusion)

    def pull_back_to_entity(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Read the pulled back field at the included point.

        The field is pulled back onto its cell by the walk that reaches
        this node, and the restriction is what says where on that cell it
        is read. Having said it, the restriction is done: what is left is
        the field itself, evaluated.
        """
        field = node_map.get(self._field, self._field)
        return evaluate_at(field, self.image_point)

    def __repr__(self) -> str:
        """Representation."""
        return f"Restriction({self._field!r}, {self._domain!r}, {self._variable!r})"


def _pieces_of(domain: AbstractParametrizedDomain) -> tuple[AbstractParametrizedDomain, ...]:
    """The pieces of a union, or the domain itself when it is one piece."""
    if isinstance(domain, UnionDomain):
        return domain.pieces
    return (domain,)


def _inclusion_of(
    domain: AbstractParametrizedDomain,
) -> tuple[AbstractCellularDomain, AbstractEntity, AbstractParametrization]:
    """Read a piece as a sub-entity sitting in a cell, mapped on by the cell's own map.

    A sub-entity sitting in a cell is a domain whose points are named in
    that cell's coordinates, which is what its coordinate domain says.

    Args:
        domain: One piece of the domain a field is restricted to

    Returns:
        The domain it is parametrized from, the cell its points are named in,
        and the map applied to it

    Raises:
        ValueError: If the piece is not a sub-entity of a cell mapped on
    """
    composite = unlabelled(domain)
    if not isinstance(composite, ComposesADomain):
        raise ValueError(f"{domain!r} is not a sub-entity of a cell, mapped on.")
    inner = composite.inner_domain
    cell_map = composite.map
    sits_in = inner.coordinate_domain
    if not isinstance(sits_in, EntityDomain):
        raise ValueError(f"{inner!r} is named in {sits_in!r}, not in a cell's coordinates.")
    (cell,) = sits_in.cell_types
    return inner, cell, cell_map


def restrict(field: AbstractExpression, domain: AbstractParametrizedDomain) -> Restriction:
    """Restrict a field to a domain made of sub-entities of its own.

    Args:
        field: The field being restricted
        domain: The domain it is restricted to

    Returns:
        The field read at the point that domain's parametrization gives
    """
    return Restriction(field, domain)
