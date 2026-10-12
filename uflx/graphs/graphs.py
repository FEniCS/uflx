"""Graphs."""

from __future__ import annotations

import functools
import inspect
import weakref
from abc import ABCMeta
from collections import deque
from collections.abc import Callable, Hashable, Iterable, Mapping
from collections.abc import Set as AbstractSet
from enum import Enum
from hashlib import blake2b
from typing import Any, ClassVar, Protocol, TypeVar, cast, runtime_checkable

import networkx as nx

_T = TypeVar("_T")


class _Interning(ABCMeta):
    """Return the value already built from the same class and arguments.

    The key is read after the outermost ``__init__`` returns. A
    constructor that normalises its arguments interns on the normalised
    ones. The table holds its values weakly and keeps nothing alive.
    """

    _values: weakref.WeakValueDictionary[Any, Any] = weakref.WeakValueDictionary()
    _rebuilt: ClassVar[set[type]] = set()
    check_every_value: ClassVar[bool] = False
    _checking: ClassVar[bool] = False

    def __call__(cls: type[_T], *args: Any, **kwargs: Any) -> _T:
        """Build a value or return the equal one already built.

        The first value of each class is rebuilt from its ``init_args``
        and must come back as the same object. This catches a class that
        stores or reads back its arguments wrongly, once per class. Every
        value is checked while ``check_every_value`` is set, as
        :mod:`uflx.testing` does.

        Raises:
            TypeError: If the first value of a class does not rebuild to
                itself.
        """
        value = super().__call__(*args, **kwargs)
        # Fingerprint first so that an argument that cannot be digested is refused with a reason.
        value.fingerprint
        key = (cls, _interning_key(value.init_args))
        interned = _Interning._values.setdefault(key, value)
        if not _Interning._checking and (
            _Interning.check_every_value or cls not in _Interning._rebuilt
        ):
            # A rebuild builds values too, and they need no check of their own.
            _Interning._checking = True
            try:
                _check_rebuilds(interned)
            except Exception:
                # A value that fails its check must not be returned by a later build.
                if interned is value:
                    del _Interning._values[key]
                raise
            finally:
                _Interning._checking = False
            _Interning._rebuilt.add(cls)
        return interned


def _check_rebuilds(value: Any) -> None:
    """Check that a value rebuilt from its ``init_args`` is the same object."""
    cls = type(value)
    try:
        rebuilt = cls(*value.init_args)
    except TypeError as error:
        raise TypeError(
            f"{cls.__qualname__} cannot be rebuilt from its init_args {describe(value)}. "
            f"Store every constructor argument as _name, or make init_args match the constructor."
        ) from error
    if rebuilt is not value:
        raise TypeError(
            f"{cls.__qualname__} rebuilt from its init_args is {describe(rebuilt)}, "
            f"not {describe(value)}. Store every constructor argument as _name in the form "
            f"it is compared in."
        )


def _interning_key(value: Any) -> Any:
    """Key an argument by what its fingerprint digests.

    A value is already interned and keys as itself. Any other argument
    keys as its type name and its repr, or its hash if it has no repr of
    its own. Interning and hashing then always agree. Arguments that are
    equal but print differently stay apart, such as ``1`` and ``1.0`` or
    ``0.0`` and ``-0.0``. Otherwise a value would print and generate code
    as whichever was built first, and its hash would depend on it too.
    """
    if isinstance(value, Immutable):
        return value
    if isinstance(value, tuple):
        return ("tuple", *map(_interning_key, value))
    return _leaf(value)


_DERIVED = ("__eq__", "__ne__", "__hash__", "__repr__", "__str__", "successors", "fingerprint")


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
    the attributes the constructor stores its arguments in.

    A subclass, in UFLx or in a consumer, gets these checks at runtime.
    Defining ``__eq__``, ``__ne__``, ``__hash__``, ``__repr__``,
    ``__str__``, ``successors`` or ``fingerprint`` is refused when the
    class is defined. An argument
    that is a list, a set, a dict or an object known only by its address
    is refused when a value is built. The first value of each class must
    rebuild from its ``init_args`` to itself. :mod:`uflx.testing` checks
    every value a test suite builds.
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
        """Refuse a class that defines its own identity, and wrap its ``__init__`` to freeze.

        Only an ``__init__`` defined in the class body is wrapped. An
        inherited one was already wrapped in its parent. Only the
        outermost ``__init__`` freezes, so an ``__init__`` that calls its
        parent's can still set attributes afterwards. The names of the
        constructor's arguments are recorded here for ``init_args``.

        Raises:
            TypeError: If the class body defines a method that would give
                the class an identity other than its ``init_args``.
        """
        super().__init_subclass__(**kwargs)
        for name in _DERIVED:
            # Defining __eq__ alone sets __hash__ to None in the class body.
            if name in cls.__dict__ and not (name == "__hash__" and cls.__dict__[name] is None):
                raise TypeError(
                    f"{cls.__qualname__} defines {name}. A value's identity and printing come "
                    f"from its init_args alone, so a class defines none of {', '.join(_DERIVED)}."
                )
        init = cls.__dict__.get("__init__")
        if init is not None:
            setattr(cls, "__init__", _freezing(init))
            cls._init_parameters = _parameters_of(cls, init)

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

    A value is digested by its class and arguments. A tuple is digested
    item by item in order. Any other object is digested by its class and
    its own repr, or its own hash if it has no repr of its own. Either
    must be the same in every process.

    Raises:
        TypeError: If an object is a list, a set or a mapping, or has
            neither its own repr nor its own hash.
    """
    if isinstance(value, Immutable):
        return value.fingerprint
    if isinstance(value, tuple):
        return digest("tuple", *map(fingerprint, value))
    return digest(*_leaf(value))


def _leaf(value: Any) -> tuple[str, str]:
    """Describe an argument that is neither a value nor a tuple, as fingerprint and interning do.

    Raises:
        TypeError: If the argument is a list, a set or a mapping, or has
            neither its own repr nor its own hash.
    """
    cls = type(value)
    if isinstance(value, list):
        raise TypeError(
            "A list can change after the value holding it is built. Store a tuple instead."
        )
    if isinstance(value, AbstractSet | Mapping):
        raise TypeError(
            f"{_qualified_name(cls)} prints in an order that differs between processes "
            f"or between equal values. Store a tuple instead."
        )
    if cls.__repr__ is not object.__repr__:
        return _qualified_name(cls), repr(value)
    if cls.__hash__ is not None and cls.__hash__ is not object.__hash__:
        return _qualified_name(cls), str(hash(value))
    raise TypeError(
        f"{_qualified_name(cls)} prints and hashes by its address, so it cannot be "
        f"digested the same way in every process. Give it a repr or a hash."
    )


def _values_in(argument: Any) -> Iterable[Immutable]:
    if isinstance(argument, Immutable):
        yield argument
    elif isinstance(argument, tuple):
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
    if isinstance(value, Immutable):
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


class NodeOrder(Enum):
    """A node order."""

    roots_first = 0
    leaves_first = 1


def print_node(graph: Graph, node: GraphNode, indentation: int = 0):
    """Print a graph using the node as the root node.

    A node with several parents prints in full the first time with a
    label such as ``[1]``, and as its label alone after that. The output
    is one line per edge however much the graph shares.
    """
    reachable = graph.descendants(node)
    shared = {n for n in reachable if sum(p in reachable for p in graph.predecessors(n)) > 1}
    labels: dict[GraphNode, str] = {}
    indent = "    " * indentation
    print(indent + repr(node))
    unprinted = [(child, indent, False) for child in reversed(list(graph.successors(node)))]
    if unprinted:
        unprinted[0] = (unprinted[0][0], indent, True)
    while unprinted:
        child, prefix, is_last = unprinted.pop()
        connector = "└── " if is_last else "├── "
        if child in labels:
            print(prefix + connector + labels[child])
            continue
        if child in shared:
            labels[child] = f"[{len(labels) + 1}]"
            print(f"{prefix}{connector}{labels[child]} {child!r}")
        else:
            print(prefix + connector + repr(child))
        extension = prefix + ("    " if is_last else "│   ")
        children = list(graph.successors(child))
        for i, grandchild in reversed(list(enumerate(children))):
            unprinted.append((grandchild, extension, i == len(children) - 1))


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
