from dataclasses import dataclass

from xdsl.context import Context
from xdsl.dialects import func, tensor
from xdsl.dialects.builtin import (
    DYNAMIC_INDEX,
    AffineMapAttr,
    ArrayAttr,
    DenseArrayBase,
    ModuleOp,
    TensorType,
    i64,
)
from xdsl.dialects.linalg import ops as linalg
from xdsl.dialects.linalg.attrs import IteratorType, IteratorTypeAttr
from xdsl.ir.affine import AffineDimExpr, AffineExpr, AffineMap
from xdsl.passes import ModulePass
from xdsl.rewriter import Rewriter

from ffcy.transforms.linalg_builders import reassociation_attr
from ffcy.transforms.split_cells import is_batched


def interleave_type(t: TensorType, lanes: int) -> TensorType:
    return TensorType(t.element_type, [DYNAMIC_INDEX, *t.get_shape()[1:], lanes])


def interleave_map(m: AffineMap) -> AffineMap:
    """Add a last loop over lanes, which follows the cell loop d0 to the end."""
    results = list(m.results)
    if results and isinstance(results[0], AffineDimExpr) and results[0].position == 0:
        results.append(AffineExpr.dimension(m.num_dims))
    return AffineMap(m.num_dims + 1, m.num_symbols, tuple(results))


def interleave_reassociation(groups: ArrayAttr, rank: int) -> ArrayAttr:
    """The lanes are a last group of their own, after the `rank` axes of the groups."""
    indices = [[i.value.data for i in group] for group in groups]
    return reassociation_attr([*indices, [rank]])


def interleave_cells(module: ModuleOp, lanes: int) -> None:
    values = [
        value
        for op in module.walk()
        for value in (*op.results, *(a for r in op.regions for b in r.blocks for a in b.args))
        if is_batched(value.type)
    ]
    for value in values:
        Rewriter.replace_value_with_new_type(value, interleave_type(value.type, lanes))

    for op in module.walk():
        if isinstance(op, linalg.GenericOp):
            op.properties["indexing_maps"] = ArrayAttr(
                AffineMapAttr(interleave_map(m.data)) for m in op.indexing_maps
            )
            op.properties["iterator_types"] = ArrayAttr(
                [*op.iterator_types, IteratorTypeAttr(IteratorType.PARALLEL)]
            )
        if isinstance(op, tensor.CollapseShapeOp) and is_batched(op.src.type):
            rank = len(op.src.type.get_shape()) - 1
            op.properties["reassociation"] = interleave_reassociation(op.reassociation, rank)
        if isinstance(op, tensor.ExpandShapeOp) and is_batched(op.result.type):
            rank = len(op.result.type.get_shape()) - 1
            op.properties["reassociation"] = interleave_reassociation(op.reassociation, rank)
            shape = list(op.static_output_shape.get_values())
            op.properties["static_output_shape"] = DenseArrayBase.from_list(i64, [*shape, lanes])
        if isinstance(op, func.FuncOp):
            op.update_function_type()


# Interleaves the cells of every tensor on cells in blocks of `lanes`, so that a
# tensor of shape [cells, ...] becomes [cells / lanes, ..., lanes]. Each lane holds
# a different cell, and every loop over a block runs over the lanes innermost, which
# a CPU vectorises across cells. The number of cells must be a multiple of `lanes`.
@dataclass(frozen=True)
class InterleaveCells(ModulePass):
    name = "ffcy-interleave-cells"

    lanes: int

    def apply(self, ctx: Context, op: ModuleOp) -> None:
        interleave_cells(op, self.lanes)
