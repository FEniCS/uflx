from dataclasses import dataclass

from xdsl.context import Context
from xdsl.dialects import func, tensor
from xdsl.dialects.builtin import ModuleOp
from xdsl.passes import ModulePass


def expand_results(function: func.FuncOp) -> None:
    (ret,) = [op for op in function.walk() if isinstance(op, func.ReturnOp)]
    for i, value in enumerate(ret.operands):
        if isinstance(value.owner, tensor.CollapseShapeOp):
            ret.operands[i] = value.owner.src
    function.update_function_type()


# A flat E-vector and its tensor product shape share one memory layout, so callers
# can pass either view. Returning the uncollapsed tensor lets bufferization write
# the result straight into a caller provided buffer.
@dataclass(frozen=True)
class ExpandResults(ModulePass):
    name = "ffcy-expand-results"

    def apply(self, ctx: Context, op: ModuleOp) -> None:
        for function in [f for f in op.ops if isinstance(f, func.FuncOp)]:
            expand_results(function)
