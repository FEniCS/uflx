"""Graphs."""

from __future__ import annotations

import functools
import inspect
import math
import weakref
from abc import ABCMeta
from collections import deque
from collections.abc import Callable, Hashable, Iterable, Mapping
from collections.abc import Set as AbstractSet
from enum import Enum
from hashlib import blake2b
from typing import Any, Protocol, TypeVar, cast, runtime_checkable

import networkx as nx

_T = TypeVar("_T")


class _Interning(ABCMeta):
    """Return the value already built from the same class and arguments.

    The key is read after the outermost ``__init__`` returns. A
    constructor that normalises its arguments interns on the normalised
    ones. The table holds its values weakly and keeps nothing alive.
    """

    _values: weakref.WeakValueDictionary[Any, Any] = weakref.WeakValueDictionary()

    def __call__(cls: type[_T], *args: Any, **kwargs: Any) -> _T:
        """Build a value or return the equal one already built."""
        value = super().__call__(*args, **kwargs)
        # Fingerprint first so that an argument that cannot be digested is refused with a reason.
        value.fingerprint
        return _Interning._values.setdefault((cls, _interning_key(value.init_args)), value)


def _interning_key(value: Any) -> Any:
    """Key an argument so that only identical values share a key.

    A value is already interned and keys as itself. A number keys with
    its type and a float also with its sign. This keeps ``1`` apart from
    ``1.0`` and ``0.0`` apart from ``-0.0``. Otherwise a value would
    print and generate code as whichever was built first.
    """
    if isinstance(value, tuple | list):
        return (type(value), *(_interning_key(item) for item in value))
    if isinstance(value, float):
        return (type(value), value, math.copysign(1.0, value))
    if isinstance(value, int | complex):
        return (type(value), value)
    return value


class Immutable(metaclass=_Interning):
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

    A value's identity is its class and its ``init_args``. Building a
    value from the same class and arguments as a live one returns that
    one. Equal values are the same object and equality is ``is``. Hash
    and repr also come from here. ``init_args`` is read back from
    the attributes the constructor stores its arguments in. A class
    should avoid defining ``__eq__``, ``__hash__`` or ``__repr__`` itself.
    Any other identity is one that rebuilding a node can silently change.
    """

    _hash: int
    _init_parameters: tuple[str, ...] | str = ()

    @property
    def init_args(self) -> tuple[Any, ...]:
        """The arguments used to initialise this object.

        Each constructor argument ``name`` is read from the attribute
        ``_name``. A class that stores an argument any other way defines
        this itself.
        """
        try:
            return self.__dict__["_init_args"]
        except KeyError:
            pass
        parameters = type(self)._init_parameters
        if isinstance(parameters, str):
            raise TypeError(parameters)
        arguments = tuple(_stored_argument(self, name) for name in parameters)
        # Attributes may still change during construction, so cache only once frozen.
        if self.__dict__.get("_frozen", False):
            object.__setattr__(self, "_init_args", arguments)
        return arguments

    def __repr__(self) -> str:
        """Representation.

        This is :func:`describe` at its default depth of one.
        """
        return describe(self)

    @property
    def successors(self) -> set[GraphNode]:
        """The values this one was built from.

        These are the values among the arguments. A graph has an edge to
        each, so a walk or a rewrite reaches everything a value depends on.
        """
        return set(_values_in(self.init_args))

    @property
    def fingerprint(self) -> bytes:
        """A digest of the class and arguments that is the same in every process.

        Python salts the hash of a string in each process, and a class
        hashes by its address. The fingerprint uses neither, so MPI ranks
        agree on it whatever their hash seeds.
        """
        unfinished: list[Immutable] = [self]
        while unfinished:
            value = unfinished[-1]
            if "_fingerprint" in value.__dict__:
                unfinished.pop()
                continue
            # Fingerprint the arguments first so that a deep graph does not recurse.
            arguments = [a for a in _values_in(value.init_args) if "_fingerprint" not in a.__dict__]
            if arguments:
                unfinished.extend(arguments)
                continue
            object.__setattr__(
                value,
                "_fingerprint",
                digest(_qualified_name(type(value)), fingerprint(value.init_args)),
            )
            unfinished.pop()
        return self.__dict__["_fingerprint"]

    def __hash__(self) -> int:
        """Hash by the fingerprint, and keep the result.

        Equality is identity. The hash still comes from the arguments so
        that sets and dicts of values iterate in the same order in every
        process.
        """
        try:
            return self.__dict__["_hash"]
        except KeyError:
            value = int.from_bytes(self.fingerprint[:8], "little", signed=True)
            object.__setattr__(self, "_hash", value)
            return value

    def __reduce__(self) -> tuple[type, tuple[Any, ...]]:
        """Copy or unpickle by building again to get the interned value."""
        return type(self), self.init_args

    def __init_subclass__(cls, **kwargs: Any) -> None:
        """Wrap the class's own ``__init__`` to freeze and ``__hash__`` to keep.

        Only methods defined in the class body are wrapped. An inherited
        method was already wrapped in its parent. Only the outermost
        ``__init__`` freezes, so an ``__init__`` that calls its parent's
        can still set attributes afterwards. The names of the
        constructor's arguments are recorded here for ``init_args``.
        """
        super().__init_subclass__(**kwargs)
        init = cls.__dict__.get("__init__")
        if init is not None:
            setattr(cls, "__init__", _freezing(init))
            cls._init_parameters = _parameters_of(cls, init)
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


def _parameters_of(cls: type, init: Callable[..., None]) -> tuple[str, ...] | str:
    """Get the names of a constructor's arguments.

    Returns:
        The names in order, or a message saying why ``init_args`` cannot
        be read back. Only positional arguments can be read back.
    """
    names = []
    for parameter in list(inspect.signature(init).parameters.values())[1:]:
        if parameter.kind is not parameter.POSITIONAL_OR_KEYWORD:
            return (
                f"{cls.__name__}.__init__ takes {parameter.name} as {parameter.kind.description}, "
                f"so its init_args cannot be read back. Define init_args."
            )
        names.append(parameter.name)
    return tuple(names)


def digest(*parts: str | bytes) -> bytes:
    """Digest parts, each prefixed by its length so that parts cannot run together."""
    h = blake2b(digest_size=16)
    for part in parts:
        data = part.encode() if isinstance(part, str) else part
        h.update(len(data).to_bytes(8, "little"))
        h.update(data)
    return h.digest()


def fingerprint(value: Any) -> bytes:
    """Digest a value or an argument the same way in every process.

    A value is digested by its class and arguments. A tuple or list is
    digested item by item in order. Any other object is digested by its
    class and its own repr, or its own hash if it has no repr of its own.
    Either must be the same in every process.

    Raises:
        TypeError: If an object is a set or a mapping, or has neither its
            own repr nor its own hash.
    """
    if isinstance(value, Immutable):
        return value.fingerprint
    if isinstance(value, tuple | list):
        return digest(type(value).__name__, *map(fingerprint, value))
    cls = type(value)
    if isinstance(value, AbstractSet | Mapping):
        raise TypeError(
            f"{_qualified_name(cls)} prints in an order that differs between processes "
            f"or between equal values. Store a tuple instead."
        )
    if cls.__repr__ is not object.__repr__:
        return digest(_qualified_name(cls), repr(value))
    if cls.__hash__ is not None and cls.__hash__ is not object.__hash__:
        return digest(_qualified_name(cls), str(hash(value)))
    raise TypeError(
        f"{_qualified_name(cls)} prints and hashes by its address, so it cannot be "
        f"digested the same way in every process. Give it a repr or a hash."
    )


def _values_in(argument: Any) -> Iterable[Immutable]:
    if isinstance(argument, Immutable):
        yield argument
    elif isinstance(argument, tuple | list):
        for item in argument:
            yield from _values_in(item)


def _qualified_name(cls: type) -> str:
    return f"{cls.__module__}.{cls.__qualname__}"


def describe(value: Immutable, depth: int = 1) -> str:
    """Print a value as its class and its arguments, to a given depth.

    An argument nested deeper than ``depth`` prints as ``ClassName(...)``.
    An argument whose own arguments are all plain values prints in full
    at any depth. The output stays short however large the graph behind
    the value.

    Args:
        value: The value to print
        depth: How many levels of arguments to print in full

    Returns:
        The class name and arguments, such as
        ``Product((Coefficient(...), TestFunction(...)))``.
    """
    try:
        arguments = value.init_args
    except TypeError:
        return object.__repr__(value)
    if depth <= 0 and not _has_plain_arguments(arguments):
        return f"{type(value).__name__}(...)"
    return (
        f"{type(value).__name__}({', '.join(_describe_argument(a, depth - 1) for a in arguments)})"
    )


def _describe_argument(value: Any, depth: int) -> str:
    """Print one argument of a value for :func:`describe`."""
    if isinstance(value, tuple):
        items = [_describe_argument(i, depth) for i in value]
        return f"({items[0]},)" if len(items) == 1 else f"({', '.join(items)})"
    if isinstance(value, list):
        return f"[{', '.join(_describe_argument(i, depth) for i in value)}]"
    if isinstance(value, Immutable) and type(value).__repr__ is Immutable.__repr__:
        return describe(value, depth)
    return repr(value)


def _has_plain_arguments(arguments: Any) -> bool:
    """Check that no argument has arguments of its own, looking one level down only."""
    if isinstance(arguments, tuple | list):
        return all(_has_plain_arguments(a) for a in arguments)
    if isinstance(arguments, Immutable):
        try:
            return arguments.init_args == ()
        except TypeError:
            return False
    return True


def _stored_argument(value: Any, name: str) -> Any:
    """Get the attribute ``_name`` a constructor stored its argument ``name`` in.

    Raises:
        TypeError: If ``_name`` is not set.
    """
    try:
        return value.__dict__[f"_{name}"]
    except KeyError:
        raise TypeError(
            f"{type(value).__name__} does not store its argument {name} as _{name}. "
            f"Store it there, or define init_args."
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
    graph.add_root_node(node)
    unexpanded = deque([node])
    while unexpanded:
        n = unexpanded.popleft()
        for successor in n.successors:
            if successor not in graph:
                unexpanded.append(successor)
            graph.add_edge(n, successor)
    return graph


def as_graph(node: GraphNode) -> Graph:
    """Convert to item to a graph."""
    if isinstance(node, RepresentedByGraph):
        return node.graph
    return generate_graph(node)
