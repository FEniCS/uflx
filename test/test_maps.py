"""Test maps."""

from uflx.expressions import Integer
from uflx.maps import (
    BlockedValueMap,
    IdentityValueMap,
)
from uflx.tensors import Vector


def test_identity_map(lagrange_element):
    """Test an identity map."""
    id_map = IdentityValueMap()

    five = Integer(5)
    assert id_map.push_forward(five) == five
    assert id_map.pull_back(five) == five


def test_blocked_map(lagrange_element):
    """Test a blocked map."""
    id_map = IdentityValueMap()
    b_map = BlockedValueMap(id_map, (4,))

    v = Vector([Integer(i) for i in range(4)])
    assert b_map.push_forward(v) == v
    assert b_map.pull_back(v) == v
