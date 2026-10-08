from dataclasses import dataclass

from xdsl.builder import Builder, ImplicitBuilder
from xdsl.context import Context
from xdsl.dialects import func
from xdsl.dialects.builtin import DYNAMIC_INDEX, IntegerAttr, ModuleOp, TensorType, f64, i32, i64
from xdsl.passes import ModulePass
from xdsl.rewriter import InsertPoint

from ffcy.dialects import uflx

SCATTER_DOFMAP = "ffcy.scatter_dofmap"


def assemble(function: func.FuncOp, scatter: str) -> None:
    block = function.body.block
    coefficients = [op for op in block.ops if isinstance(op, uflx.CoefficientOp)]
    (u,) = {op.input for op in coefficients}
    ndofs = u.type.get_shape()[1]
    # A dof is shared by at most 2^d cells of a box mesh.
    width = 2 ** uflx.tdim(coefficients[0].element)

    vector = block.insert_arg(TensorType(f64, [DYNAMIC_INDEX]), u.index)
    vector.name_hint = "x"
    dofmap = block.insert_arg(TensorType(i32, [DYNAMIC_INDEX, ndofs]), len(block.args))
    dofmap.name_hint = "dofmap"

    with ImplicitBuilder(Builder(InsertPoint.at_start(block))):
        gathered = uflx.GatherOp(vector, dofmap, coefficients[0].element).output
        u.replace_all_uses_with(gathered)
    block.erase_arg(u)

    if scatter == "atomic":
        function.attributes[SCATTER_DOFMAP] = IntegerAttr(dofmap.index, i64)
    elif scatter == "transpose":
        transpose = block.insert_arg(TensorType(i32, [DYNAMIC_INDEX, width]), len(block.args))
        transpose.name_hint = "transpose"
        ret = block.last_op
        assert isinstance(ret, func.ReturnOp)
        with ImplicitBuilder(Builder(InsertPoint.before(ret))):
            ret.operands = [uflx.ScatterAddOp(v, transpose).output for v in ret.operands]
    function.update_function_type()


# Turns the action on E-vectors into one on assembled vectors, y = G^T A G x, where
# G gathers the dofs of each cell through the dofmap. The scatter sums each dof's
# entries through the transposed dofmap, so the result does not depend on the order
# in which cells finish. Otherwise the action returns its values on cells. With
# scatter "none" the caller sums them, and with "atomic" ffcy.scatter_dofmap names the
# dofmap argument through which code generation adds them atomically.
@dataclass(frozen=True)
class Assemble(ModulePass):
    name = "ffcy-assemble"

    scatter: str = "transpose"

    def apply(self, ctx: Context, op: ModuleOp) -> None:
        for function in op.ops:
            if isinstance(function, func.FuncOp) and function.sym_name.data == "action":
                assemble(function, self.scatter)
