# Copyright (C) 2025 Matthew Scroggs and Garth N. Wells
#
# This file is part of UFLx (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

"""Entities.

A entity is a single item in a mesh (usually a polytope). Entities of (topological)
dimension 0, 1, 2 and 3 are called vertices, edges, faces and volumes (respectively).

In a mesh, the entities of the highest dimension are called cells. The codimension
of an entity in this mesh is equal to the dimension of the entity subtracted from
the dimension of the cell. Entities of codimension 0, 1, 2 and 3 are called cells,
facets, ridges and peaks (respectively).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from uflx.domains import AbstractParametrization


class AbstractEntity(ABC):
    """Abstract base class for entities."""

    @abstractmethod
    def __eq__(self, other) -> bool:
        """Check if this entity is equal to another entity."""

    @property
    @abstractmethod
    def topological_dimension(self) -> int:
        """Topological dimension of the entity."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Name of the entity type."""

    @abstractmethod
    def sub_entities(self, dim: int) -> list[AbstractEntity]:
        """Get a list of sub-entities of a given dimension.

        Args:
            dim: Dimension of the sub-entities to get.

        Returns:
            A list of sub-entities of the given dimension.
        """

    @abstractmethod
    def sub_entity_vertices(self, dim: int) -> list[list[int]]:
        """Get lists of the vertices of sub-entities of a given dimension.

        Args:
            dim: Dimension of the sub-entities to get.

        Returns:
            A list of lists of vertices of sub-entities of the given dimension.
        """

    def sub_entity_inclusion(self, dim: int, index: int) -> AbstractParametrization:
        """Get the map that carries a sub-entity's coordinates into this entity's.

        Reference geometry, which this class otherwise does not hold:
        sub_entity_vertices says which vertices a sub-entity has and not
        where they are, so the inclusion of facet 0 of a triangle as
        ``s -> (1 - s, s)`` cannot be read off the combinatorics. It is
        supplied by whoever defines the entity, as the coordinates of an
        entity domain's points are.

        Args:
            dim: Dimension of the sub-entity.
            index: Which sub-entity of that dimension.

        Returns:
            The map from the sub-entity's coordinate domain into this one's.

        Raises:
            NotImplementedError: If this entity does not say where its
                sub-entities are.
        """
        raise NotImplementedError(
            f"{self!r} does not say where its sub-entities are. Where they sit is "
            f"reference geometry, which the definer of the entity supplies."
        )

    def sub_entity_count(self, dim: int) -> int:
        """Get the number of sub-entities of a given dimension.

        Args:
            dim: Dimension of the sub-entities to get.

        Returns:
            The number of sub-entities of the given dimension.
        """
        return len(self.sub_entities(dim))

    @property
    def is_simplex(self) -> bool:
        """Whether this entity is a simplex, ie has exactly dim + 1 vertices.

        Computed generically from sub_entities(0) and topological_dimension, so
        subclasses do not need to declare this themselves.
        """
        return len(self.sub_entities(0)) == self.topological_dimension + 1

    @abstractmethod
    def __hash__(self):
        """Hash."""
