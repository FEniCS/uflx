"""Algorithms related to push forward and pull back maps."""

from typing import Protocol, runtime_checkable

from uflx.algorithms.reconstruct import reconstruct_node
from uflx.graphs import GraphNode, as_graph


@runtime_checkable
class PullBackToEntity(Protocol):
    """Pull a node back to the entity's coordinates."""

    def pull_back_to_entity(self, node_map: dict[GraphNode, GraphNode]) -> GraphNode:
        """Pull the node back to the entity's coordinates."""


def pull_back_to_entity(
    expression: GraphNode,
) -> GraphNode:
    """Pull terms in integrals back to the entity's coordinates."""
    node_map: dict[GraphNode, GraphNode] = {}
    for node in as_graph(expression).ordered_nodes():
        if isinstance(node, PullBackToEntity):
            node_map[node] = node.pull_back_to_entity(node_map)
        elif any(a in node_map for a in node.successors):
            node_map[node] = reconstruct_node(node, node_map)

    return node_map.get(expression, expression)
