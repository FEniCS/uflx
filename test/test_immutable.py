"""Test that everything in a graph is built once, never changed and hashed once."""

import copy
import os
import pickle
import subprocess
import sys
import textwrap

import pytest
from conftest import LagrangeElement, Triangle

from uflx import Coefficient, coordinate_element, dx, function_space
from uflx.expressions import Rational, RealScalar
from uflx.graphs.graphs import Immutable


@pytest.fixture
def values(lagrange_element):
    """One instance of each family that is hashed by value."""
    cell = Triangle()
    element = LagrangeElement(cell, 1)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    integral = Coefficient(space) * dx
    expression = RealScalar(2.0) * RealScalar(3.0)
    return [cell, element, domain, space, dx, integral, expression]


def test_every_family_is_a_value(values):
    """Each is immutable once built, whatever class defined its hash."""
    for value in values:
        assert isinstance(value, Immutable)
        with pytest.raises(AttributeError, match="immutable"):
            value.anything = 1
        with pytest.raises(AttributeError, match="immutable"):
            del value.anything


def test_every_familys_hash_is_kept(values):
    """The hash a class computes is computed once and kept on the instance.

    An integral is hashed by the graph its ``__init__`` generates. Its
    arguments are set by then, so it arrives with its hash kept.
    """
    for value in values:
        by_value = hasattr(value, "init_args") or any(
            "__hash__" in vars(cls)
            for cls in type(value).__mro__
            if cls is not object and cls is not Immutable
        )
        first = hash(value)
        assert hash(value) == first
        assert ("_hash" in value.__dict__) == by_value


def test_a_subclass_is_not_equal_to_its_parent():
    """Equal arguments make equal values only within one class."""

    class Pair(Immutable):
        def __init__(self, a, b):
            self._a = a
            self._b = b

    class Twin(Pair):
        pass

    assert Pair(1, 2) == Pair(1, 2)
    assert Pair(1, 2) != Twin(1, 2)


def test_an_argument_stored_under_another_name_is_refused():
    """The first hash names the argument and where it should be stored."""

    class Renamed(Immutable):
        def __init__(self, a):
            self._other = a

    with pytest.raises(TypeError, match="does not store its argument a as _a"):
        hash(Renamed(1))


def test_the_librarys_values_are_equal_by_what_they_were_built_from(lagrange_element):
    """Library values built the same way are equal without hand-written equality."""
    from uflx.domains import RD
    from uflx.maps import IdentityReferenceMap

    element = lagrange_element("triangle", 1, (2,))
    assert RD(2) == RD(2) and RD(2) != RD(3)
    assert IdentityReferenceMap() == IdentityReferenceMap()
    assert coordinate_element(element) == coordinate_element(element)
    assert function_space(coordinate_element(element), element) == function_space(
        coordinate_element(element), element
    )


def test_a_value_prints_as_what_it_was_built_from():
    """A repr shows the arguments and never an address."""
    from uflx.domains import RD
    from uflx.maps import IdentityReferenceMap

    assert repr(RD(2)) == "RD(2)"
    assert repr(IdentityReferenceMap()) == "IdentityReferenceMap()"
    assert repr(IdentityReferenceMap()) == repr(IdentityReferenceMap())


def test_a_repr_does_not_grow_with_its_graph():
    """Printing a deeply shared expression costs a bounded number of levels."""
    from uflx.geometry import SpatialCoordinate

    e = SpatialCoordinate(2).component(0)
    for _ in range(50):
        e = e / e - e
    assert repr(e) == "Sum((Product(...), Neg(...)))"


def test_no_class_in_uflx_prints_itself():
    """Guards against a handwritten repr in UFLx drifting from the value's init_args."""
    import importlib
    import pkgutil

    import uflx

    for module in pkgutil.walk_packages(uflx.__path__, "uflx."):
        importlib.import_module(module.name)

    def subclasses(cls):
        for sub in cls.__subclasses__():
            yield sub
            yield from subclasses(sub)

    for cls in subclasses(Immutable):
        if cls.__module__.startswith("uflx."):
            for method in ("__repr__", "__str__"):
                assert method not in vars(cls), (cls, method)


_HASHES = textwrap.dedent("""
    from uflx import dx
    from uflx.domains import RD
    from uflx.expressions import RealScalar
    from uflx.graphs.graphs import Immutable
    from uflx.maps import IdentityReferenceMap

    class Labelled(Immutable):
        def __init__(self, labels):
            self._labels = labels

    labels = Labelled(("left", "right", "top", "bottom"))
    values = (dx, RD(2), IdentityReferenceMap(), RealScalar(1.0) * RealScalar(2.0), labels)
    print([hash(v) for v in values])
""")
_ORDER = textwrap.dedent("""
    from uflx.geometry import SpatialCoordinate
    from uflx.graphs import as_graph

    x = SpatialCoordinate(2)
    e = x.component(0) * x.component(1) + x.component(1) / x.component(0)
    print(as_graph(e).ordered_nodes())
""")


@pytest.mark.parametrize("script", [_HASHES, _ORDER], ids=["hashes", "order"])
def test_every_process_agrees_whatever_its_hash_seed(script):
    """Guards against hashes or graph order differing between MPI ranks."""
    printed = {
        seed: subprocess.run(
            [sys.executable, "-c", script],
            env={**os.environ, "PYTHONHASHSEED": seed},
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        for seed in ("1", "2")
    }
    assert printed["1"] == printed["2"]


def test_a_deep_value_hashes_without_recursing():
    """Guards against the first hash of a deep graph overflowing the stack."""
    from uflx.geometry import SpatialCoordinate

    def deep():
        e = SpatialCoordinate(2).component(0)
        for _ in range(2000):
            e = e / e - e
        return e

    assert hash(deep()) == hash(deep())


class _Holder(Immutable):
    def __init__(self, a):
        self._a = a


class _Opaque:
    pass


@pytest.mark.parametrize(
    "argument, message",
    [
        (_Opaque(), "prints and hashes by its address"),
        (frozenset({"left", "right"}), "Store a tuple instead"),
        ({"left": 1, "right": 2}, "Store a tuple instead"),
        (["left", "right"], "Store a tuple instead"),
    ],
    ids=["address", "frozenset", "dict", "list"],
)
def test_an_argument_that_prints_differently_between_processes_is_refused(argument, message):
    """Guards against hashing by an address or by an order that differs between processes."""
    with pytest.raises(TypeError, match=message):
        hash(_Holder(argument))


def test_labels_in_a_tuple_are_hashed_in_order():
    """Guards against two differently ordered sequences of labels hashing alike."""
    assert hash(_Holder(("left", "right"))) != hash(_Holder(("right", "left")))


@pytest.mark.parametrize("name", ["__eq__", "__ne__", "__hash__", "successors", "fingerprint"])
def test_a_class_that_defines_its_own_identity_is_refused(name):
    """Guards against a consumer's class giving itself an identity other than its init_args."""
    with pytest.raises(TypeError, match=f"defines {name}"):
        type("Custom", (Immutable,), {name: lambda self: None})


def test_a_class_whose_init_args_lose_an_argument_is_refused():
    """Guards against a class whose values cannot be rebuilt from their init_args."""

    class Lossy(Immutable):
        def __init__(self, a, b):
            self._a = a
            self._b = b

        @property
        def init_args(self):
            return (self._a,)

    with pytest.raises(TypeError, match="cannot be rebuilt from its init_args"):
        Lossy(1, 2)


def test_a_class_that_alters_its_argument_is_refused():
    """Guards against a constructor whose stored argument rebuilds to a different value."""

    class Shifting(Immutable):
        def __init__(self, n):
            self._n = n + 1

    with pytest.raises(TypeError, match="rebuilt from its init_args is"):
        Shifting(1)


def _shared(levels):
    from uflx.geometry import SpatialCoordinate

    x = SpatialCoordinate(2)
    e = x.component(0)
    for _ in range(levels):
        e = e / e - e
    return x, e


@pytest.mark.parametrize(
    "build_a, build_b",
    [
        (lambda: RealScalar(2.0) * RealScalar(3.0), lambda: RealScalar(2.0) * RealScalar(3.0)),
        (lambda: Rational(2, 4), lambda: Rational(1, 2)),
        (lambda: _shared(1000)[1], lambda: _shared(1000)[1]),
    ],
)
def test_a_value_built_twice_is_one_object(build_a, build_b):
    """Guards against equal values built apart staying distinct objects."""
    assert build_a() is build_b()


def test_a_rebuilt_value_is_the_original():
    """Guards against a rebuild skipping the table or re-expanding shared nodes."""
    from uflx.algorithms import replace
    from uflx.geometry import SpatialCoordinate

    x, e = _shared(1000)
    y = SpatialCoordinate(3)
    assert replace(replace(e, {x: y}), {y: x}) is e


@pytest.mark.parametrize("first, second", [(1, 1.0), (True, 1), (0.0, -0.0), (-0.0, 0.0)])
def test_numbers_that_print_differently_stay_apart(first, second):
    """Guards against a number printing as an equal one that was built first."""
    kept = RealScalar(first)
    built = RealScalar(second)
    assert built is not kept
    assert repr(built) == f"RealScalar({second!r})"


def test_the_table_keeps_nothing_alive():
    """Guards against the table holding every value ever built."""
    import gc
    import weakref

    ref = weakref.ref(RealScalar(123.25))
    gc.collect()
    assert ref() is None


@pytest.mark.parametrize(
    "duplicate",
    [copy.copy, copy.deepcopy, lambda value: pickle.loads(pickle.dumps(value))],
)
def test_a_copy_is_the_value_itself(duplicate):
    """Guards against a copy that is unequal to its original."""
    value = RealScalar(2.0) * RealScalar(3.0)
    assert duplicate(value) is value


def test_a_shared_node_prints_once_then_by_its_label(capsys):
    """Guards against printing a node again under every parent it has."""
    from uflx.geometry import SpatialCoordinate
    from uflx.graphs import as_graph

    a = SpatialCoordinate(2).component(0)
    as_graph(a / a - a).print()
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 7
    assert sum(line.endswith("[1] SingleSpatialCoordinate(2, 0)") for line in lines) == 1
    assert sum(line.endswith("── [1]") for line in lines) == 2


def test_a_deeply_shared_graph_prints_one_line_per_edge(capsys):
    """Guards against printing a path count that doubles with every level of sharing."""
    from uflx.geometry import SpatialCoordinate
    from uflx.graphs import as_graph

    e = SpatialCoordinate(2).component(0)
    for _ in range(1000):
        e = e / e - e
    as_graph(e).print()
    assert len(capsys.readouterr().out.splitlines()) == 6 * 1000 + 1
