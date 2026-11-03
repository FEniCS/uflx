"""UFL #492 analogues: a single form with unit weight.

Run with the spec interpreter's src directory on PYTHONPATH:
    python -m pytest -q -o addopts='' reproductions/test_ufl_492.py
"""

import pytest
from uflx_compiler.algebra import _normalizer
from uflx_oracle import evaluate
from uflx_spec_interpreter import dual_space, form, function, num
from uflx_spec_interpreter._src.arithmetic import _total
from uflx_spec_interpreter._src.nodes import Product
from uflx_toy.spaces import space


@pytest.mark.parametrize("opaque", [False, True])
def test_single_form_with_unit_weight(opaque):
    """A unit-weight form agrees with its unweighted operand."""
    V = space("pr492", scalar_type="real")
    v = function(V, "v")
    f = function(dual_space(V), "f") if opaque else form(num(3) * v, (v,))
    env = {f: lambda x: 3 * x} if opaque else {}

    # The sum helper adds no wrapper for a single term.
    assert _total((f,)) is f

    for weighted in (num(1) * f, f * num(1)):
        # Arithmetic construction is lazy; collection is a separate operation.
        assert isinstance(weighted, Product)
        assert weighted is not f
        simplified = _normalizer()(weighted)
        assert simplified is f
        assert simplified.argument_spaces == (V,)
        assert evaluate(weighted, env)((2.0,)) == 6.0
        assert evaluate(simplified, env)((2.0,)) == 6.0
