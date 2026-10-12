"""Test that everything in a graph is built once, never changed and hashed once.

:mod:`test_expressions` tests this for expressions. Every other family
follows the same rule through the same mixin, so each is tested here.
"""

import pytest
from conftest import LagrangeElement, Triangle

from uflx import Coefficient, coordinate_element, dx, function_space
from uflx.graphs.graphs import Immutable


@pytest.fixture
def values(lagrange_element):
    """One instance of each family that is hashed by value."""
    cell = Triangle()
    element = LagrangeElement(cell, 1)
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, element)
    integral = Coefficient(space) * dx
    return [cell, element, domain, space, dx, integral]


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
        by_value = any(
            "__hash__" in vars(cls)
            for cls in type(value).__mro__
            if cls is not object and cls is not Immutable
        )
        first = hash(value)
        assert hash(value) == first
        assert ("_hash" in value.__dict__) == by_value


def test_a_derived_value_may_still_be_cached():
    """Any derived value can be cached through object.__setattr__."""
    cell = Triangle()
    object.__setattr__(cell, "_memo", 42)
    assert cell.__dict__["_memo"] == 42
    with pytest.raises(AttributeError, match="immutable"):
        object.__getattribute__(cell, "__setattr__")("_memo", 43)


def test_only_a_classs_own_methods_are_wrapped():
    """A subclass that inherits its parent's methods is frozen and hashed by them.

    Wrapping happens where a method is defined, so it is not done twice.
    A subclass with no methods of its own still gets both.
    """

    class Plain(Triangle):
        pass

    assert "__hash__" not in Plain.__dict__
    assert "__init__" not in Plain.__dict__
    plain = Plain()
    assert hash(plain) == hash(plain)
    assert "_hash" in plain.__dict__
    with pytest.raises(AttributeError, match="immutable"):
        plain.anything = 1
