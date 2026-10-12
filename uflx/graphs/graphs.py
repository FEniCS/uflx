"""Graphs."""

from __future__ import annotations

import functools
from collections.abc import Callable, Hashable, Iterable
from enum import Enum
from typing import Any, Protocol, cast, runtime_checkable

import networkx as nx


class Immutable:
    """A value that is built once, never changed and hashed once.

    Everything in a graph is hashed and compared by what it was built
    from. A node that changed after being built would sit under the
    wrong hash and compare unequal to itself. This mixin refuses to set
    or delete an attribute once the outermost ``__init__`` has returned.
    It also keeps the hash the class computes on the instance after the
    first call. A changed node is a new node, which
    :func:`uflx.algorithms.reconstruct_node` builds.

    A derived value may still be cached on an instance through
    ``object.__setattr__``, since it changes no argument. The kept hash
    is cached this way.

    A value's identity is its class and its ``init_args``. A class that
    defines ``init_args`` gets its equality, hash and repr from here. A
    class whose identity is something else writes its own and says why.
    """

    _hash: int

    def __eq__(self, other: Any) -> bool:
        """Check for equality.

        Two values are equal when they have the same class and equal
        arguments. The kept hashes are compared first, which tells most
        unequal values apart without walking their arguments.
        """
        if self is other:
            return True
        if type(self) is not type(other) or hash(self) != hash(other):
            return False
        return _arguments_of(self) == _arguments_of(other)

    def __repr__(self) -> str:
        """Representation.

        This is the class name and its arguments. A class without
        ``init_args`` prints as a plain object.
        """
        try:
            arguments = _arguments_of(self)
        except TypeError:
            return object.__repr__(self)
        return f"{type(self).__name__}({', '.join(repr(a) for a in arguments)})"

    def __hash__(self) -> int:
        """Hash the class and arguments, and keep the result."""
        try:
            return self.__dict__["_hash"]
        except KeyError:
            value = hash((type(self), *_arguments_of(self)))
            object.__setattr__(self, "_hash", value)
            return value

    def __init_subclass__(cls, **kwargs: Any) -> None:
        """Wrap the class's own ``__init__`` to freeze and ``__hash__`` to keep.

        Only methods defined in the class body are wrapped. An inherited
        method was already wrapped in its parent. Only the outermost
        ``__init__`` freezes, so an ``__init__`` that calls its parent's
        can still set attributes afterwards.
        """
        super().__init_subclass__(**kwargs)
        init = cls.__dict__.get("__init__")
        if init is not None:
            setattr(cls, "__init__", _freezing(init))
        hash_ = cls.__dict__.get("__hash__")
        if hash_ is not None:
            setattr(cls, "__hash__", _kept(hash_))

    def __setattr__(self, name: str, value: Any) -> None:
        """Set an attribute during construction only."""
        if self.__dict__.get("_frozen", False):
            raise AttributeError(
                f"{self.__class__.__name__} is immutable. Build a new one instead."
            )
        object.__setattr__(self, name, value)

    def __delattr__(self, name: str) -> None:
        """Refuse to delete an attribute."""
        raise AttributeError(f"{self.__class__.__name__} is immutable.")


def _arguments_of(value: Any) -> tuple[Any, ...]:
    """Get the arguments a value was built from.

    Raises:
        TypeError: If the value has no ``init_args``.
    """
    try:
        return value.init_args
    except AttributeError:
        raise TypeError(
            f"{type(value).__name__} has no init_args. Define init_args, or write "
            f"__eq__ and __hash__."
        ) from None


def _freezing(init: Callable[..., None]) -> Callable[..., None]:
    """Wrap ``__init__`` to freeze the instance when the outermost call returns."""

    @functools.wraps(init)
    def freezing_init(self: Any, *args: Any, **kwargs: Any) -> None:
        outermost = "_initialising" not in self.__dict__
        if outermost:
            object.__setattr__(self, "_initialising", True)
        init(self, *args, **kwargs)
        if outermost:
            del self.__dict__["_initialising"]
            object.__setattr__(self, "_frozen", True)

    return freezing_init


def _kept(hash_: Callable[[Any], int]) -> Callable[[Any], int]:
    """Wrap ``__hash__`` to compute it once and keep it on the instance."""

    @functools.wraps(hash_)
    def kept_hash(self: Any) -> int:
        try:
            return self.__dict__["_hash"]
        except KeyError:
            value = hash_(self)
            object.__setattr__(self, "_hash", value)
            return value

    return kept_hash


class NodeOrder(Enum):
    """A node order."""

    roots_first = 0
    leaves_first = 1


def _print_tree(graph: Graph, node: GraphNode, prefix: str, is_last: bool):
    connector = "└── " if is_last else "├── "
    print(prefix + connector + repr(node))
    successors = list(graph.successors(node))
    extension = "    " if is_last else "│   "
    for i, child in enumerate(successors):
        _print_tree(graph, child, prefix + extension, i == len(successors) - 1)


def print_node(graph: Graph, node: GraphNode, indentation: int = 0):
    """Print a graph using the node as the root node."""
    indent = "    " * indentation
    print(indent + repr(node))
    successors = list(graph.successors(node))
    for i, child in enumerate(successors):
        _print_tree(graph, child, indent, i == len(successors) - 1)


@runtime_checkable
class GraphNode(Hashable, Protocol):
    """A node in a graph."""

    @property
    def successors(self) -> set[GraphNode]:
        """The successors of this node."""

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object."""


class Graph(nx.DiGraph):
    """An acyclic directed graph."""

    _root: GraphNode | None

    def __init__(self, *args, **kwargs):
        """Initialise."""
        self._root = None
        super().__init__(*args, **kwargs)
        assert nx.is_directed_acyclic_graph(self)

    @property
    def root(self) -> GraphNode:
        """Get the root node of the graph."""
        assert self._root is not None
        return self._root

    @property
    def has_root(self) -> bool:
        """Check if this graph has a root node."""
        return self._root is not None

    def set_root(self, node: GraphNode):
        """Set the root node of the graph."""
        self._root = node

    def add_root_node(self, node: GraphNode, *args, **kwargs):
        """Add a new node to the graph and set it to be the root."""
        self.add_node(node, *args, **kwargs)
        self._root = node

    def print(self):
        """Print a graph."""
        print_node(self, self.root)

    def subgraph_of_node(self, node: GraphNode) -> Graph:
        """Get the subgraph with the input node as the root node."""
        return as_graph(node)

    def ordered_nodes(self, order: NodeOrder = NodeOrder.leaves_first) -> Iterable[GraphNode]:
        """Iterate through the ordered graph nodes."""
        nodes = cast(list[GraphNode], list(nx.topological_sort(self)))
        match order:
            case NodeOrder.roots_first:
                return nodes
            case NodeOrder.leaves_first:
                return list(reversed(nodes))
            case _:
                raise ValueError("Invalid node order")

    def descendants(self, node: GraphNode) -> set[GraphNode]:
        """Get all descendants of a node."""
        return cast(set[GraphNode], {node, *nx.descendants(self, node)})

    def is_dag(self) -> bool:
        """Check if this graph is a directed acyclic graph."""
        return nx.is_directed_acyclic_graph(self)


@runtime_checkable
class RepresentedByGraph(Protocol):
    """An object whose construction is represented by a graph."""

    @property
    def graph(self) -> Graph:
        """The graph that represents this object."""


def generate_graph(node: GraphNode) -> Graph:
    """Generate the graph that represents the construction of a node."""
    graph = Graph()

    added_nodes = {node}
    graph.add_root_node(node)
    while len(added_nodes) > 0:
        for n in added_nodes:
            for successor in n.successors:
                graph.add_node(successor)
                graph.add_edge(n, successor)
        added_nodes = set().union(*[n.successors for n in added_nodes])

    return graph


def as_graph(node: GraphNode) -> Graph:
    """Convert to item to a graph."""
    if isinstance(node, RepresentedByGraph):
        return node.graph
    return generate_graph(node)
