from basix import CellType
from xdsl.context import Context
from xdsl.dialects import arith
from xdsl.dialects.utils.fast_math import FastMathFlag

from ffcy.pipeline import compile_form
from ffcy.transforms.allow_reciprocal import AllowReciprocal
from forms import FORMS


def test_allow_reciprocal():
    """Every division allows a reciprocal and keeps the flags it had."""
    module = compile_form(FORMS["laplace"][0](1, CellType.hexahedron), 3)
    AllowReciprocal().apply(Context(), module)
    divisions = [op for op in module.walk() if isinstance(op, arith.DivfOp)]
    assert divisions
    for division in divisions:
        assert division.fastmath.data == {FastMathFlag.ALLOW_CONTRACT, FastMathFlag.ALLOW_RECIP}
