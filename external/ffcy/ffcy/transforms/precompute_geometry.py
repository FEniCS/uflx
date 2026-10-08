from collections.abc import Sequence
from dataclasses import dataclass

from xdsl.builder import Builder, ImplicitBuilder
from xdsl.context import Context
from xdsl.dialects import func, tensor
from xdsl.dialects.builtin import DYNAMIC_INDEX, ModuleOp, TensorType, f64
from xdsl.ir import BlockArgument, Operation, SSAValue
from xdsl.ir.affine import AffineExpr, AffineMap
from xdsl.passes import ModulePass
from xdsl.rewriter import InsertPoint

from ffcy.dialects import uflx
from ffcy.transforms.linalg_builders import dims, empty, expand, generic, num_cells_of


def depends_on(ops: list[Operation], values: set[SSAValue]) -> set[Operation]:
    """The ops whose results depend on any of `values`, in program order."""
    dependent = set()
    for op in ops:
        if any(v in values for v in op.operands):
            dependent.add(op)
            values |= set(op.results)
    return dependent


def qdata_type(values: Sequence[SSAValue]) -> TensorType:
    """Fields at the quadrature points packed as [cells, fields, nq, ..., nq]."""
    types = {v.type for v in values}
    if len(types) != 1:
        raise NotImplementedError(f"Precomputed fields must share one type, not {types}")
    (field_type,) = types
    return TensorType(f64, [DYNAMIC_INDEX, len(values), *field_type.get_shape()[1:]])


def pack(values: Sequence[SSAValue]) -> SSAValue:
    num_cells = num_cells_of(values[0])
    first, *rest = values[0].type.get_shape()[1:]
    fields = [expand(num_cells, v, [[1, first], *([n] for n in rest)]) for v in values]
    return tensor.ConcatOp(fields, 1, qdata_type(values)).result


def unpack(qdata: SSAValue, i: int) -> SSAValue:
    shape = qdata.type.get_shape()[2:]
    rank = len(shape) + 1
    field = AffineMap(
        rank, 0, (AffineExpr.dimension(0), AffineExpr.constant(i), *dims(*range(1, rank)))
    )
    return generic(
        [qdata],
        [field, AffineMap.identity(rank)],
        empty(num_cells_of(qdata), shape),
        lambda value, _: value,
    )


def precompute_geometry(module: ModuleOp, function: func.FuncOp) -> None:
    block = function.body.block
    ops = list(block.ops)
    applied = {op.input for op in ops if isinstance(op, uflx.CoefficientOp)}
    data = set(block.args) - applied
    # Ops that only need the mesh data move to setup. Constants are copied as needed.
    hoisted = depends_on(ops, set(data)) - depends_on(ops, set(applied))
    hoisted = {op for op in hoisted if not isinstance(op, func.ReturnOp)}
    if not hoisted:
        return
    cut = list(
        dict.fromkeys(
            result
            for op in ops
            if op in hoisted
            for result in op.results
            if any(use.operation not in hoisted for use in result.uses)
        )
    )
    inputs = [arg for arg in block.args if any(u.operation in hoisted for u in arg.uses)]

    packed_type = qdata_type(cut)
    setup = func.FuncOp(
        f"{function.sym_name.data}_setup", ([v.type for v in inputs], [packed_type])
    )
    mapping: dict[SSAValue, SSAValue] = dict(zip(inputs, setup.args, strict=True))
    for value, arg in zip(inputs, setup.args, strict=True):
        arg.name_hint = value.name_hint
    setup_block = setup.body.block
    for op in ops:
        if op not in hoisted:
            continue
        for operand in op.operands:
            if operand not in mapping and not operand.owner.operands:
                constant = operand.owner.clone()
                setup_block.add_op(constant)
                mapping.update(zip(operand.owner.results, constant.results, strict=True))
        clone = op.clone(value_mapper=mapping)
        setup_block.add_op(clone)
        mapping.update(zip(op.results, clone.results, strict=True))
    with ImplicitBuilder(setup_block):
        func.ReturnOp(pack([mapping[v] for v in cut]))
    module.body.block.insert_op_before(setup, function)

    qdata = block.insert_arg(packed_type, 0)
    qdata.name_hint = "qdata"
    with ImplicitBuilder(Builder(InsertPoint.at_start(block))):
        fields = [unpack(qdata, i) for i in range(len(cut))]
    for value, field in zip(cut, fields, strict=True):
        value.replace_all_uses_with(field)
    for op in reversed(ops):
        if op in hoisted:
            block.erase_op(op)
    for arg in list(block.args):
        if isinstance(arg, BlockArgument) and not arg.uses:
            block.erase_arg(arg)
    function.update_function_type()


# Moves everything that does not depend on the applied E-vectors, the inputs of
# uflx.coefficient ops, into a setup function. Its results are packed into one
# qdata tensor, so every form has the signatures setup(mesh data) -> qdata and
# action(qdata, E-vectors).
@dataclass(frozen=True)
class PrecomputeGeometry(ModulePass):
    name = "ffcy-precompute-geometry"

    def apply(self, ctx: Context, op: ModuleOp) -> None:
        for function in [f for f in op.ops if isinstance(f, func.FuncOp)]:
            precompute_geometry(op, function)
