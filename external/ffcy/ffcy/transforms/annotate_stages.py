# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

from dataclasses import dataclass
from itertools import count

from xdsl.context import Context
from xdsl.dialects import func, tensor
from xdsl.dialects.builtin import IntegerAttr, ModuleOp, i64
from xdsl.dialects.linalg import ops as linalg
from xdsl.ir import Operation
from xdsl.ir.affine import AffineDimExpr
from xdsl.passes import ModulePass

from ffcy.transforms.linalg_builders import PENCIL

# A stage is one phase of a pencil schedule, ending with a barrier. Ops that end a
# stage, whose results later stages read, are tagged ffcy.sink with the stage. The
# other ops of a stage, and the fills they write into, are tagged ffcy.stage and
# fused into the threads of its sinks.
STAGE = "ffcy.stage"
SINK = "ffcy.sink"
BUFFER = "ffcy.buffer"


def pencil(op: linalg.GenericOp) -> int:
    return op.attributes[PENCIL].value.data


def reads_own_pencil(producer: linalg.GenericOp, consumer: linalg.GenericOp) -> bool:
    """Whether each thread of the consumer only reads the pencil it computed itself.

    Both ops loop over the same axis, and the consumer reads the producer's result
    at its own indices along every other axis.
    """
    axis = pencil(producer)
    if axis != pencil(consumer):
        return False
    for operand, indexing_map in zip(consumer.operands, consumer.indexing_maps, strict=False):
        if operand.owner is producer:
            for i, expr in enumerate(indexing_map.data.results):
                if i != axis and not (isinstance(expr, AffineDimExpr) and expr.position == i):
                    return False
    return True


def grid(op: linalg.GenericOp) -> tuple[int, ...]:
    """The shape of op's result without its pencil axis, the pencils its threads take."""
    shape = list(op.results[0].type.get_shape())
    del shape[pencil(op)]
    return tuple(shape)


def users(op: Operation) -> list[Operation]:
    """Ops that read op's results, ignoring size queries."""
    return [
        use.operation
        for result in op.results
        for use in result.uses
        if not isinstance(use.operation, tensor.DimOp)
    ]


def fusible(op: linalg.GenericOp) -> bool:
    """Whether every user of op reads only the pencils its own threads computed."""
    return bool(users(op)) and all(
        isinstance(user, linalg.GenericOp)
        and any(operand.owner is op for operand in user.inputs)
        and reads_own_pencil(op, user)
        for user in users(op)
    )


def stages(
    generics: list[linalg.GenericOp], is_fusible: dict[linalg.GenericOp, bool]
) -> dict[Operation, int]:
    """Stages as soon as possible, with fused producers in the stage of their users.

    The sinks of a stage take the same grid of pencils, so sinks as soon as possible
    whose grids differ, such as those of geometry of another degree, start stages of
    their own in turn.
    """
    depth: dict[Operation, int] = {}
    for op in generics:
        depth[op] = max(
            (
                depth[p] + (0 if is_fusible[p] else 1)
                for p in (operand.owner for operand in op.inputs)
                if isinstance(p, linalg.GenericOp)
            ),
            default=0,
        )
    key = {op: (depth[op], grid(op)) for op in generics}
    for op in reversed(generics):
        if is_fusible[op]:
            key[op] = max(key[user] for user in users(op))
    order = {
        k: i for i, k in enumerate(sorted(set(key.values()), key=lambda k: (k[0], first(key, k))))
    }
    return {op: order[key[op]] for op in generics}


def first(key: dict[Operation, tuple[int, tuple[int, ...]]], k: tuple) -> int:
    """The position of the first op with key `k`."""
    return next(i for i, op in enumerate(key) if key[op] == k)


def annotate_stages(function: func.FuncOp) -> None:
    generics = [op for op in function.walk() if isinstance(op, linalg.GenericOp)]
    is_fusible = {op: fusible(op) for op in generics}

    # Producers that cannot join their users' stage end their own, so repeat until no
    # more producers are ruled out.
    while True:
        stage = stages(generics, is_fusible)
        ruled_out = [
            op
            for op in generics
            if is_fusible[op] and len({stage[user] for user in users(op)}) != 1
        ]
        if not ruled_out:
            break
        for op in ruled_out:
            is_fusible[op] = False

    for op in generics:
        tag = STAGE if is_fusible[op] else SINK
        op.attributes[tag] = IntegerAttr(stage[op], i64)
        # Fills shared by several ops are cloned, so each stage can fuse its own.
        for i, value in enumerate(op.outputs, start=len(op.inputs)):
            if isinstance(fill := value.owner, linalg.FillOp):
                clone = fill.clone()
                clone.attributes[STAGE] = IntegerAttr(stage[op], i64)
                fill.parent_block().insert_op_before(clone, fill)
                op.operands[i] = clone.results[0]
    for fill in [op for op in function.walk() if isinstance(op, linalg.FillOp)]:
        if not any(fill.results[0].uses):
            fill.parent_block().erase_op(fill)

    # Each op gets an empty tensor of its own, tagged so that CSE keeps them apart.
    # Otherwise ops of a stage would share a destination, and fusion would write one
    # op's results into another's.
    buffers = count()
    for empty in [op for op in function.walk() if isinstance(op, tensor.EmptyOp)]:
        for use in list(empty.tensor.uses):
            clone = empty.clone()
            clone.attributes[BUFFER] = IntegerAttr(next(buffers), i64)
            empty.parent_block().insert_op_before(clone, empty)
            use.operation.operands[use.index] = clone.tensor
        empty.parent_block().erase_op(empty)


# Groups the ops of each function into the stages of a pencil schedule.
@dataclass(frozen=True)
class AnnotateStages(ModulePass):
    name = "ffcy-annotate-stages"

    def apply(self, ctx: Context, op: ModuleOp) -> None:
        for function in [f for f in op.ops if isinstance(f, func.FuncOp)]:
            annotate_stages(function)
