"""Test that everything in a graph is built once, never changed and hashed once."""

import os
import subprocess
import sys

import pytest
from conftest import LagrangeElement, Triangle

from uflx import Coefficient, coordinate_element, dx, function_space
from uflx.expressions import RealScalar
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


def test_no_class_in_uflx_defines_its_own_identity():
    """Every value in UFLx is identified by its init_args alone."""
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
            for method in ("__eq__", "__hash__"):
                assert method not in vars(cls), (cls, method)


_HASHES = (
    "from uflx import dx\n"
    "from uflx.domains import RD\n"
    "from uflx.expressions import RealScalar\n"
    "from uflx.maps import IdentityReferenceMap\n"
    "values = (dx, RD(2), IdentityReferenceMap(), RealScalar(1.0) * RealScalar(2.0))\n"
    "print([hash(v) for v in values])"
)
_ORDER = (
    "from uflx.geometry import SpatialCoordinate\n"
    "from uflx.graphs import as_graph\n"
    "x = SpatialCoordinate(2)\n"
    "e = x.component(0) * x.component(1) + x.component(1) / x.component(0)\n"
    "print(as_graph(e).ordered_nodes())"
)


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


def test_an_argument_known_only_by_its_address_is_refused():
    """Guards against hashing by an address that differs between processes."""

    class Opaque:
        pass

    class Holder(Immutable):
        def __init__(self, a):
            self._a = a

    with pytest.raises(TypeError, match="prints and hashes by its address"):
        hash(Holder(Opaque()))
