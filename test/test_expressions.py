"""Test what every expression node promises: to be what it was built as.

A node's arguments are its identity in the graph and its hash, so a
node that changed after being built would be found under the wrong
hash and compare unequal to itself. The base class refuses the change
rather than trusting every subclass to.
"""

import pytest

from uflx.expressions import AbstractExpression, Integer, RealScalar
from uflx.tensors import Tensor


def test_a_node_cannot_be_changed_once_built():
    """Setting or deleting an attribute after construction is refused."""
    product = RealScalar(2.0) * RealScalar(3.0)

    with pytest.raises(AttributeError, match="immutable"):
        product._items = ()
    with pytest.raises(AttributeError, match="immutable"):
        product.anything = 1
    with pytest.raises(AttributeError, match="immutable"):
        del product._items


def test_a_subclass_that_calls_its_parents_init_still_builds():
    """Only the outermost ``__init__`` freezes, so cooperative construction works."""

    class Doubled(Tensor):
        def __init__(self, entries):
            super().__init__(entries)
            self.doubled = [e * Integer(2) for e in entries]

    node = Doubled([RealScalar(1.0), RealScalar(2.0)])
    assert len(node.doubled) == 2
    assert isinstance(node, AbstractExpression)
    with pytest.raises(AttributeError, match="immutable"):
        node.doubled = []


def test_hashing_is_not_a_change():
    """The kept hash is written after construction, and is the one allowed write."""
    node = RealScalar(1.0) + RealScalar(2.0)
    assert hash(node) == hash(node)
    assert node == RealScalar(1.0) + RealScalar(2.0)
