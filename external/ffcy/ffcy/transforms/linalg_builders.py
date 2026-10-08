# Builders for linalg, tensor and arith ops on batched tensors with a leading cell axis.

from collections.abc import Callable, Sequence
from functools import reduce

import numpy as np
import numpy.typing as npt
from xdsl.builder import ImplicitBuilder
from xdsl.dialects import arith, tensor
from xdsl.dialects.builtin import (
    DYNAMIC_INDEX,
    AffineMapAttr,
    ArrayAttr,
    DenseIntOrFPElementsAttr,
    FloatAttr,
    IndexType,
    IntegerAttr,
    TensorType,
    UnitAttr,
    f64,
    i64,
)
from xdsl.dialects.linalg import ops as linalg
from xdsl.dialects.linalg.attrs import IteratorType, IteratorTypeAttr
from xdsl.dialects.utils.fast_math import FastMathFlag
from xdsl.ir import Block, Region, SSAValue
from xdsl.ir.affine import AffineExpr, AffineMap

TABLE = "ffcy.table"
# The axis each thread loops over in a pencil schedule. Contractions loop over the
# contracted axis. Pointwise ops loop over the last axis by default, where evaluation
# ends and integration starts, so a thread can pass on the pencil it computed. Fills
# are left untagged so that a schedule can fuse them into the op that consumes them.
PENCIL = "ffcy.pencil"

# Lets LLVM fuse multiplies and adds into FMAs.
CONTRACT = arith.FastMathFlagsAttr([FastMathFlag.ALLOW_CONTRACT])


def constant(array: npt.NDArray[np.float64]) -> SSAValue:
    tensor_type = TensorType(f64, array.shape)
    values = DenseIntOrFPElementsAttr.from_list(tensor_type, array.ravel().tolist())
    op = arith.ConstantOp(values)
    # Lets a schedule find the tables and move them into GPU kernels.
    op.attributes[TABLE] = UnitAttr()
    return op.result


def num_cells_of(x: SSAValue) -> SSAValue:
    # A collapse keeps the cell axis, so its source is asked instead.
    if isinstance(x.owner, tensor.CollapseShapeOp):
        x = x.owner.src
    zero = arith.ConstantOp(IntegerAttr(0, IndexType()))
    return tensor.DimOp(x, zero).result


def empty(num_cells: SSAValue, shape: Sequence[int]) -> SSAValue:
    """An output for pointwise ops, which never read it."""
    return tensor.EmptyOp([num_cells], TensorType(f64, [DYNAMIC_INDEX, *shape])).tensor


def zeros(num_cells: SSAValue, shape: Sequence[int]) -> SSAValue:
    empty = tensor.EmptyOp([num_cells], TensorType(f64, [DYNAMIC_INDEX, *shape]))
    zero = arith.ConstantOp(FloatAttr(0.0, f64))
    return linalg.FillOp([zero.result], [empty.tensor]).res[0]


def dims(*positions: int) -> tuple[AffineExpr, ...]:
    return tuple(AffineExpr.dimension(p) for p in positions)


def generic(
    inputs: Sequence[SSAValue],
    maps: Sequence[AffineMap],
    output: SSAValue,
    body: Callable[..., SSAValue],
    num_reductions: int = 0,
    pencil: int | None = None,
) -> SSAValue:
    block = Block(arg_types=[v.type.element_type for v in (*inputs, output)])
    with ImplicitBuilder(block):
        linalg.YieldOp(body(*block.args))
    num_parallel = maps[-1].num_dims - num_reductions
    iterator_types = [IteratorType.PARALLEL] * num_parallel + [IteratorType.REDUCTION] * (
        num_reductions
    )
    op = linalg.GenericOp(
        inputs,
        [output],
        Region(block),
        [AffineMapAttr(m) for m in maps],
        [IteratorTypeAttr(t) for t in iterator_types],
        [output.type],
    )
    if pencil is None:
        pencil = len(maps[-1].results) - 1
    op.attributes[PENCIL] = IntegerAttr(pencil, i64)
    return op.res[0]


def contract(
    num_cells: SSAValue, table: npt.NDArray[np.float64], x: SSAValue, axis: int
) -> SSAValue:
    """Contract `axis` of x with the second index of table."""
    shape = list(x.type.get_shape())
    rank = len(shape)
    shape[axis] = table.shape[0]
    k = rank
    output_dims = list(range(rank))
    x_dims = [k if d == axis else d for d in output_dims]
    maps = [
        AffineMap(rank + 1, 0, dims(axis, k)),
        AffineMap(rank + 1, 0, dims(*x_dims)),
        AffineMap(rank + 1, 0, dims(*output_dims)),
    ]
    return generic(
        [constant(table), x],
        maps,
        zeros(num_cells, shape[1:]),
        lambda a, b, acc: arith.AddfOp(acc, arith.MulfOp(a, b, CONTRACT), CONTRACT).result,
        num_reductions=1,
        pencil=axis,
    )


def contract_axes(
    num_cells: SSAValue, tables: Sequence[npt.NDArray[np.float64]], x: SSAValue
) -> SSAValue:
    """Contract axis i + 1 of x with tables[i] for each i."""
    for axis, table in enumerate(tables, start=1):
        x = contract(num_cells, table, x, axis)
    return x


def reassociation_attr(groups: Sequence[Sequence[int]]) -> ArrayAttr:
    return ArrayAttr(ArrayAttr(IntegerAttr(i, i64) for i in group) for group in groups)


def expand(num_cells: SSAValue, x: SSAValue, groups: Sequence[Sequence[int]]) -> SSAValue:
    """Split each non-cell axis of x into the sizes in groups."""
    shape = [DYNAMIC_INDEX, *(n for group in groups for n in group)]
    # Undoing a collapse returns its source, so ops stay directly connected.
    if isinstance(x.owner, tensor.CollapseShapeOp) and x.owner.src.type.get_shape() == tuple(shape):
        return x.owner.src
    reassociation, start = [[0]], 1
    for group in groups:
        reassociation.append(list(range(start, start + len(group))))
        start += len(group)
    return tensor.ExpandShapeOp(
        x,
        [num_cells],
        reassociation_attr(reassociation),
        shape,
        TensorType(x.type.element_type, shape),
    ).result


def collapse(x: SSAValue) -> SSAValue:
    """Merge all non-cell axes of x."""
    shape = x.type.get_shape()
    return tensor.CollapseShapeOp.build(
        operands=[x],
        result_types=[TensorType(f64, [DYNAMIC_INDEX, int(np.prod(shape[1:]))])],
        properties={"reassociation": reassociation_attr([[0], list(range(1, len(shape)))])},
    ).result


def multiply_all(values: Sequence[SSAValue]) -> SSAValue:
    return reduce(lambda a, b: arith.MulfOp(a, b, CONTRACT).result, values)


def add_all(values: Sequence[SSAValue]) -> SSAValue:
    return reduce(lambda a, b: arith.AddfOp(a, b, CONTRACT).result, values)
