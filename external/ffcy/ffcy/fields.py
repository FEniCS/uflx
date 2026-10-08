# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

"""Arithmetic on scalar fields at quadrature points, folding known constants.

A field is either a tensor value of shape [cells, nq, ..., nq] or a float that is
the same at every point. Building ops needs an active ImplicitBuilder.
"""

from collections.abc import Sequence

from xdsl.dialects import arith, math
from xdsl.ir import SSAValue

Field = float | SSAValue


def is_constant(a: Field, value: float) -> bool:
    return isinstance(a, float) and a == value


def field(a: Field) -> SSAValue:
    if isinstance(a, SSAValue):
        return a
    raise NotImplementedError(f"Constant factors other than 0 and ±1, here {a}")


def add(a: Field, b: Field) -> Field:
    if isinstance(a, float) and isinstance(b, float):
        return a + b
    if is_constant(a, 0):
        return b
    if is_constant(b, 0):
        return a
    return arith.AddfOp(field(a), field(b)).result


def neg(a: Field) -> Field:
    if isinstance(a, float):
        return -a
    return arith.NegfOp(a).result


def sub(a: Field, b: Field) -> Field:
    if isinstance(a, float) and isinstance(b, float):
        return a - b
    if is_constant(b, 0):
        return a
    if is_constant(a, 0):
        return neg(b)
    return arith.SubfOp(field(a), field(b)).result


def mul(a: Field, b: Field) -> Field:
    if isinstance(a, float) and isinstance(b, float):
        return a * b
    if isinstance(b, float):
        a, b = b, a
    if is_constant(a, 0):
        return 0.0
    if is_constant(a, 1):
        return b
    if is_constant(a, -1):
        return neg(b)
    return arith.MulfOp(field(a), field(b)).result


def div(a: Field, b: Field) -> Field:
    if isinstance(b, float):
        return mul(a, 1 / b)
    if is_constant(a, 0):
        return 0.0
    return arith.DivfOp(field(a), b).result


def absolute(a: Field) -> Field:
    if isinstance(a, float):
        return abs(a)
    return math.AbsFOp(a).result


def total(values: Sequence[Field]) -> Field:
    result: Field = 0.0
    for v in values:
        result = add(result, v)
    return result


def product(values: Sequence[Field]) -> Field:
    result: Field = 1.0
    for v in values:
        result = mul(result, v)
    return result
