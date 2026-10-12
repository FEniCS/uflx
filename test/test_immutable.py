"""Test that everything in a graph is built once, never changed and hashed once."""

import collections
import copy
import enum
import os
import pickle
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from uflx.expressions import Rational, RealScalar
from uflx.graphs.graphs import Immutable, _Interning


def test_a_value_cannot_be_changed():
    """Guards against setting or deleting an attribute once a value is built."""
    value = _Holder(1)
    with pytest.raises(AttributeError, match="immutable"):
        value._a = 2
    with pytest.raises(AttributeError, match="immutable"):
        del value._a


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


def test_a_repr_does_not_grow_with_its_graph():
    """Printing a deeply shared expression costs a bounded number of levels."""
    from uflx.geometry import SpatialCoordinate

    e = SpatialCoordinate(2).component(0)
    for _ in range(50):
        e = e / e - e
    assert repr(e) == "Sum((Product(...), Neg(...)))"


_HASHES = textwrap.dedent("""
    import enum

    from uflx import dx
    from uflx.domains import RD
    from uflx.expressions import RealScalar
    from uflx.graphs.graphs import Immutable
    from uflx.maps import IdentityReferenceMap

    class Labelled(Immutable):
        def __init__(self, labels):
            self._labels = labels

    class Side(enum.Enum):
        LEFT = 1

    labels = Labelled(("left", "right", "top", "bottom", Side.LEFT))
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


class _Holder(Immutable):
    def __init__(self, a):
        self._a = a


class _Opaque:
    pass


class _Point(collections.namedtuple("_Point", "x y")):
    pass


class _Colour(enum.IntEnum):
    RED = 1


class _Shade(enum.IntEnum):
    RED = 1


@pytest.mark.parametrize(
    "argument",
    [
        _Opaque(),
        frozenset({"left", "right"}),
        ["left", "right"],
        np.float64(1.0),
        _Point(1, 2),
    ],
    ids=[
        "address",
        "frozenset",
        "list",
        "numpy-float",
        "namedtuple",
    ],
)
def test_an_argument_of_any_other_type_is_refused(argument):
    """Guards against accepting an argument whose description could differ between processes."""
    with pytest.raises(TypeError, match="cannot be an argument of a value"):
        _Holder(argument)


@pytest.mark.parametrize(
    "first, second",
    [
        (1, 1.0),
        (True, 1),
        (0.0, -0.0),
        (1.0, 1 + 0j),
        (0j, complex(-0.0, 0.0)),
        (_Colour.RED, _Shade.RED),
        (_Colour.RED, 1),
    ],
    ids=[
        "int-float",
        "bool-int",
        "signed-zero",
        "float-complex",
        "complex-signed-zero",
        "enums",
        "enum-int",
    ],
)
def test_accepted_arguments_that_print_differently_stay_apart(first, second):
    """Guards against interning merging accepted arguments that are equal but not the same."""
    kept = _Holder(first)
    built = _Holder(second)
    assert built is not kept
    assert built.init_args == (second,)


@pytest.mark.parametrize(
    "name",
    [
        "__eq__",
        "__ne__",
        "__hash__",
        "__repr__",
        "__str__",
        "__copy__",
        "__deepcopy__",
        "successors",
        "fingerprint",
    ],
)
def test_a_class_that_defines_its_own_identity_is_refused(name):
    """Guards against a consumer's class giving itself an identity or printing of its own."""
    with pytest.raises(TypeError, match=f"defines {name}"):
        type("Custom", (Immutable,), {name: lambda self: None})


@pytest.mark.parametrize("mixin_first", [True, False], ids=["mixin-first", "mixin-last"])
def test_an_identity_method_from_a_mixin_is_refused(mixin_first):
    """Guards against a mixin supplying an identity method that the class body does not."""

    class Loose:
        def __eq__(self, other):
            return True

    bases = (Loose, Immutable) if mixin_first else (Immutable, Loose)
    with pytest.raises(TypeError, match=r"inherits from .*Loose __eq__"):
        type("Custom", bases, {})


def test_a_class_whose_init_args_lose_an_argument_is_refused(monkeypatch):
    """Guards against a class whose values cannot be rebuilt from their init_args.

    It is refused every time, not only the first. Outside a test only the
    first value of a class is checked, so the test checks that way too.
    """
    monkeypatch.setattr(_Interning, "check_every_value", False)

    class Lossy(Immutable):
        def __init__(self, a, b):
            self._a = a
            self._b = b

        @property
        def init_args(self):
            return (self._a,)

    for _ in range(2):
        with pytest.raises(TypeError, match="cannot be rebuilt from its init_args"):
            Lossy(1, 2)


def test_a_class_that_alters_its_argument_is_refused(monkeypatch):
    """Guards against a constructor whose stored argument rebuilds to a different value.

    It is refused every time, and the refused value is never returned.
    """
    monkeypatch.setattr(_Interning, "check_every_value", False)

    class Shifting(Immutable):
        def __init__(self, n):
            self._n = n + 1

    for _ in range(2):
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


def test_a_deeply_shared_graph_prints_one_line_per_edge(capsys):
    """Guards against printing a path count that doubles with every level of sharing."""
    from uflx.geometry import SpatialCoordinate
    from uflx.graphs import as_graph

    e = SpatialCoordinate(2).component(0)
    for _ in range(1000):
        e = e / e - e
    as_graph(e).print()
    assert len(capsys.readouterr().out.splitlines()) == 6 * 1000 + 1
