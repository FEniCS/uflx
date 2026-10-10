"""Algorithms related to push forward and pull back maps."""

from typing import Protocol, runtime_checkable

from uflx.algorithms.reconstruct import reconstruct_node
from uflx.graphs import GraphNode, as_graph


@runtime_checkable
class PullBackToEntity(Protocol):
    """Pull a node back to the entity's coordinates."""

    def pull_back_to_entity(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the entity's coordinates."""


@runtime_checkable
class SplitByCellType(Protocol):
    """Split a node into one term per cell type of its domain."""

    def split_by_cell_type(self) -> GraphNode | None:
        """Split the node, or give None when it is over one cell type already."""


def pull_back_to_entity(
    expression: GraphNode,
) -> GraphNode:
    """Pull terms in integrals back to the entity's coordinates.

    An integral over a domain of several cell types is split first, into
    one integral per cell type. Each cell type has its own coordinate
    domain, so pulling back without splitting would have to pick one of
    them. The split has to come first because this walk is leaves first:
    the functions in an integrand are pulled back, onto the element of
    the cell they are on, before the integral itself is reached.
    """
    if isinstance(expression, SplitByCellType):
        split = expression.split_by_cell_type()
        if split is not None:
            return pull_back_to_entity(split)

    node_map: dict[GraphNode, GraphNode] = {}
    for node in as_graph(expression).ordered_nodes():
        if isinstance(node, PullBackToEntity):
            node_map[node] = node.pull_back_to_entity(node_map)
        elif any(a in node_map for a in node.successors):
            node_map[node] = reconstruct_node(node, node_map)

    return node_map.get(expression, expression)
