"""Simplifying expressions."""

from collections.abc import Sequence
from itertools import pairwise
from typing import Any, Protocol, runtime_checkable

from uflx.algorithms.reconstruct import reconstruct_node
from uflx.graphs import GraphNode, as_graph


@runtime_checkable
class Simplifiable(Protocol):
    """An expression that can be simplified."""

    def simplify(self) -> GraphNode:
        """Simplify this expression.

        This function should return None if no simplification can be made.
        """


@runtime_checkable
class SimplifiableInProduct(Protocol):
    """An expression that can be combined with others within a multiplication."""

    def simplified_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified product.

        This function should return None if no simplification can be made.
        """


def simplify_product_items(items: Sequence[GraphNode]) -> list[GraphNode]:
    """Simplify a list of items in a product."""
    items = list(items)
    while True:
        for i, item in enumerate(items):
            if isinstance(item, SimplifiableInProduct):
                for j, item2 in enumerate(items):
                    if i != j and (s := item.simplified_product(item2)) is not None:
                        items = [it for k, it in enumerate(items) if k not in [i, j]] + [s]
                        break
                else:
                    continue
                break
        else:
            break
    return items


@runtime_checkable
class SimplifiableInSum(Protocol):
    """An expression that can be combined with others within a multiplication."""

    def simplified_sum(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified sum.

        This function should return None if no simplification can be made.
        """


def simplify_sum_items(items: Sequence[GraphNode]) -> list[GraphNode]:
    """Simplify a list of items in a sum."""
    items = list(items)
    while True:
        for i, item in enumerate(items):
            if isinstance(item, SimplifiableInSum):
                for j, item2 in enumerate(items):
                    if i != j and (s := item.simplified_sum(item2)) is not None:
                        items = [it for k, it in enumerate(items) if k not in [i, j]] + [s]
                        break
                else:
                    continue
                break
        else:
            break
    return items


@runtime_checkable
class Commutative(Protocol):
    """An expression whose operands can be reordered without changing its value."""

    @property
    def is_commutative(self) -> bool:
        """Whether the operands of this expression can be reordered."""

    @property
    def operands(self) -> tuple[GraphNode, ...]:
        """The operands."""

    def with_operands(self, operands: Sequence[GraphNode]) -> GraphNode:
        """Create the same expression with different operands."""


def sort_key(value: Any) -> tuple:
    """A structural key that orders expressions in the same way in every run.

    The key does not use hash, as Python randomises the hashes of strings.
    """
    if isinstance(value, tuple | list):
        return ("tuple", tuple(sort_key(i) for i in value))
    if isinstance(value, GraphNode):
        return (type(value).__name__, sort_key(value.init_args))
    return (type(value).__name__, repr(value))


def sort_operands(node: GraphNode) -> GraphNode:
    """Sort the operands of a commutative expression into a canonical order.

    Expressions that differ only in the order of these operands are then equal.
    """
    if not (isinstance(node, Commutative) and node.is_commutative):
        return node
    operands = sorted(node.operands, key=sort_key)
    if operands == list(node.operands):
        return node
    return node.with_operands(operands)


def simplify(expression: GraphNode) -> GraphNode:
    """Apply simplifications to an expression."""
    graph = as_graph(expression)
    assert graph.is_dag()

    node_map: dict[GraphNode, GraphNode] = {}
    for node in graph.ordered_nodes():
        new_node = node
        if isinstance(node, Simplifiable):
            new_node = reconstruct_node(node, node_map)
            assert isinstance(new_node, Simplifiable)
            new_node = new_node.simplify()
        elif any(a in node_map for a in node.successors):
            new_node = reconstruct_node(node, node_map)
        new_node = sort_operands(new_node)
        if new_node is not node:
            node_map[node] = new_node

    return node_map.get(graph.root, graph.root)


@runtime_checkable
class SimplifiableInMatrixProduct(Protocol):
    """An expression that can be combined with others within a matrix product."""

    def simplified_matrix_product(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """


@runtime_checkable
class RightSimplifiableInMatrixProduct(Protocol):
    """An expression that can be combined with others within a matrix product."""

    def simplified_matrix_product_right(self, other: GraphNode) -> GraphNode | None:
        """Return a single expression representing the simplified matrix product.

        This function should return None if no simplification can be made.
        """


def simplify_matrix_product_items(items: Sequence[GraphNode]) -> list[GraphNode]:
    """Simplify a list of items in a matrix product."""
    items = list(items)
    while True:
        for i, (item, item2) in enumerate(pairwise(items)):
            if (
                isinstance(item, SimplifiableInMatrixProduct)
                and (s := item.simplified_matrix_product(item2)) is not None
            ):
                items = [*items[:i], s, *items[i + 2 :]]
                break
            if (
                isinstance(item2, RightSimplifiableInMatrixProduct)
                and (s := item2.simplified_matrix_product_right(item)) is not None
            ):
                items = [*items[:i], s, *items[i + 2 :]]
                break
        else:
            break
    return items
