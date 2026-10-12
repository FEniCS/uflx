"""Test that everything in a graph is built once, never changed and hashed once."""

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

    A value hashed by identity, such as a function space, keeps nothing.
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


def test_a_value_that_says_what_it_was_built_from_needs_no_equality_of_its_own():
    """Identity is the class and the arguments, so saying the arguments is enough."""

    class Pair(Immutable):
        def __init__(self, a, b):
            self.a, self.b = a, b

        @property
        def init_args(self):
            return (self.a, self.b)

    class Twin(Pair):
        pass

    assert Pair(1, 2) == Pair(1, 2)
    assert hash(Pair(1, 2)) == hash(Pair(1, 2))
    assert Pair(1, 2) != Pair(2, 1)
    assert Pair(1, 2) != Twin(1, 2)
    assert len({Pair(1, 2), Pair(1, 2), Twin(1, 2)}) == 2


def test_a_value_that_says_nothing_of_itself_has_no_identity():
    """Neither arguments nor its own equality: refused at the first hash, with what to do."""

    class Mute(Immutable):
        pass

    with pytest.raises(TypeError, match="init_args, or writes its own"):
        hash(Mute())
    with pytest.raises(TypeError, match="init_args"):
        Mute() == Mute()


def test_the_librarys_values_are_equal_by_what_they_were_built_from(lagrange_element):
    """Two of anything built the same way are the same, with no equality written by hand."""
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
    """So that nothing that hashes or compares by repr sees an address."""
    from uflx.domains import RD
    from uflx.maps import IdentityReferenceMap

    assert repr(RD(2)) == "RD(2)"
    assert repr(IdentityReferenceMap()) == "IdentityReferenceMap()"
    assert repr(IdentityReferenceMap()) == repr(IdentityReferenceMap())
