from dataclasses import dataclass

from xdsl.context import Context
from xdsl.dialects import arith
from xdsl.dialects.builtin import ModuleOp
from xdsl.dialects.utils.fast_math import FastMathFlag
from xdsl.passes import ModulePass


# Allows each division to be computed as a multiplication by the reciprocal of its
# divisor, which can differ in the last bit. Divisions by the same value that end up
# in one block then share one reciprocal, which LLVM computes once.
@dataclass(frozen=True)
class AllowReciprocal(ModulePass):
    name = "ffcy-allow-reciprocal"

    def apply(self, ctx: Context, op: ModuleOp) -> None:
        for division in op.walk():
            if isinstance(division, arith.DivfOp):
                flags = division.fastmath.data | {FastMathFlag.ALLOW_RECIP}
                division.properties["fastmath"] = arith.FastMathFlagsAttr(flags)
