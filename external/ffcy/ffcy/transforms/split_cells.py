from dataclasses import dataclass

from xdsl.context import Context
from xdsl.dialects import func, tensor
from xdsl.dialects.builtin import (
    DYNAMIC_INDEX,
    AffineMapAttr,
    ArrayAttr,
    DenseArrayBase,
    IntegerAttr,
    ModuleOp,
    TensorType,
    i64,
)
from xdsl.dialects.linalg import ops as linalg
from xdsl.dialects.linalg.attrs import IteratorType, IteratorTypeAttr
from xdsl.ir import Attribute
from xdsl.ir.affine import AffineDimExpr, AffineExpr, AffineMap
from xdsl.passes import ModulePass
from xdsl.rewriter import Rewriter

from ffcy.transforms.linalg_builders import PENCIL, reassociation_attr


def is_batched(t: Attribute) -> bool:
    """Values on cells have a leading cell axis, unlike assembled vectors of rank 1."""
    shape = t.get_shape() if isinstance(t, TensorType) else ()
    return len(shape) > 1 and shape[0] == DYNAMIC_INDEX


def split_type(t: TensorType, cells: int) -> TensorType:
    return TensorType(t.element_type, [DYNAMIC_INDEX, cells, *t.get_shape()[1:]])


def split_map(m: AffineMap) -> AffineMap:
    """Add a leading loop, replacing each use of the cell loop d0 by (d0, d1)."""
    shifted = [AffineExpr.dimension(i + 1) for i in range(m.num_dims)]
    results = []
    for expr in m.results:
        if isinstance(expr, AffineDimExpr) and expr.position == 0:
            results += [AffineExpr.dimension(0), AffineExpr.dimension(1)]
        else:
            results.append(expr.replace_dims_and_symbols(shifted, []))
    return AffineMap(m.num_dims + 1, m.num_symbols, tuple(results))


def split_reassociation(groups: ArrayAttr) -> ArrayAttr:
    """The leading cell group [0] becomes [0], [1], and later axes shift by one."""
    indices = [[i.value.data for i in group] for group in groups]
    return reassociation_attr([[0], [1], *[[i + 1 for i in group] for group in indices[1:]]])


def split_cells(module: ModuleOp, cells: int) -> None:
    values = [
        value
        for op in module.walk()
        for value in (*op.results, *(a for r in op.regions for b in r.blocks for a in b.args))
        if is_batched(value.type)
    ]
    for value in values:
        Rewriter.replace_value_with_new_type(value, split_type(value.type, cells))

    for op in module.walk():
        if isinstance(op, linalg.GenericOp):
            op.properties["indexing_maps"] = ArrayAttr(
                AffineMapAttr(split_map(m.data)) for m in op.indexing_maps
            )
            op.properties["iterator_types"] = ArrayAttr(
                [IteratorTypeAttr(IteratorType.PARALLEL), *op.iterator_types]
            )
        if PENCIL in op.attributes:
            op.attributes[PENCIL] = IntegerAttr(op.attributes[PENCIL].value.data + 1, i64)
        if isinstance(op, tensor.ExpandShapeOp | tensor.CollapseShapeOp):
            op.properties["reassociation"] = split_reassociation(op.reassociation)
        if isinstance(op, tensor.ConcatOp):
            op.properties["dim"] = IntegerAttr(op.dim.value.data + 1, i64)
        if isinstance(op, tensor.ExpandShapeOp):
            shape = list(op.static_output_shape.get_values())
            op.properties["static_output_shape"] = DenseArrayBase.from_list(
                i64, [shape[0], cells, *shape[1:]]
            )
        if isinstance(op, func.FuncOp):
            op.update_function_type()


# Splits the leading cell axis of every tensor on cells into blocks of `cells` cells,
# so that a GPU schedule sees a static number of cells per block. The number of
# cells must be a multiple of `cells`, and the memory layout is unchanged.
@dataclass(frozen=True)
class SplitCells(ModulePass):
    name = "ffcy-split-cells"

    cells: int

    def apply(self, ctx: Context, op: ModuleOp) -> None:
        split_cells(op, self.cells)
