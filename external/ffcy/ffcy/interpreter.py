# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

from math import prod

import numpy as np
import numpy.typing as npt
from xdsl.context import Context
from xdsl.dialects import arith, math, tensor
from xdsl.dialects.builtin import DYNAMIC_INDEX, DenseIntOrFPElementsAttr, ModuleOp
from xdsl.dialects.linalg import ops as linalg
from xdsl.interpreter import (
    Interpreter,
    InterpreterFunctions,
    PythonValues,
    impl,
    register_impls,
)
from xdsl.interpreters import register_implementations
from xdsl.interpreters.shaped_array import ShapedArray
from xdsl.interpreters.utils.ptr import TypedPtr


def shaped_array(data: list[float], shape: list[int]) -> ShapedArray[float]:
    return ShapedArray(TypedPtr.new_float64(data), shape)


# Implementations of ops the upstream xDSL interpreter lacks or gets wrong.
@register_impls
class MissingFunctions(InterpreterFunctions):
    @impl(arith.ConstantOp)
    def run_constant(
        self, interpreter: Interpreter, op: arith.ConstantOp, args: PythonValues
    ) -> PythonValues:
        if isinstance(op.value, DenseIntOrFPElementsAttr):
            shape = list(op.value.get_type().get_shape())
            return (shaped_array(list(op.value.get_values()), shape),)
        return (op.value.value.data,)

    @impl(tensor.ExpandShapeOp)
    def run_expand_shape(
        self, interpreter: Interpreter, op: tensor.ExpandShapeOp, args: PythonValues
    ) -> PythonValues:
        source, *dynamic_sizes = args
        sizes = iter(dynamic_sizes)
        shape = [
            next(sizes) if n == DYNAMIC_INDEX else n for n in op.static_output_shape.get_values()
        ]
        return (ShapedArray(source.data_ptr.copy(), shape),)

    @impl(tensor.CollapseShapeOp)
    def run_collapse_shape(
        self, interpreter: Interpreter, op: tensor.CollapseShapeOp, args: PythonValues
    ) -> PythonValues:
        (source,) = args
        shape = [prod(source.shape[i.value.data] for i in group) for group in op.reassociation]
        return (ShapedArray(source.data_ptr.copy(), shape),)

    @impl(tensor.ConcatOp)
    def run_concat(
        self, interpreter: Interpreter, op: tensor.ConcatOp, args: PythonValues
    ) -> PythonValues:
        arrays = [np.array(a.data).reshape(a.shape) for a in args]
        result = np.concatenate(arrays, axis=op.dim.value.data)
        return (shaped_array(result.ravel().tolist(), list(result.shape)),)

    @impl(linalg.FillOp)
    def run_fill(
        self, interpreter: Interpreter, op: linalg.FillOp, args: PythonValues
    ) -> PythonValues:
        value, output = args
        return (shaped_array([value] * len(output.data), list(output.shape)),)

    @impl(arith.NegfOp)
    def run_negf(
        self, interpreter: Interpreter, op: arith.NegfOp, args: PythonValues
    ) -> PythonValues:
        (value,) = args
        return (-value,)

    @impl(arith.DivfOp)
    def run_divf(
        self, interpreter: Interpreter, op: arith.DivfOp, args: PythonValues
    ) -> PythonValues:
        a, b = args
        return (a / b,)

    @impl(math.AbsFOp)
    def run_absf(
        self, interpreter: Interpreter, op: math.AbsFOp, args: PythonValues
    ) -> PythonValues:
        (value,) = args
        return (abs(value),)


def run(
    module: ModuleOp, name: str, *arrays: npt.NDArray[np.float64]
) -> tuple[npt.NDArray[np.float64], ...]:
    interpreter = Interpreter(module)
    register_implementations(interpreter, Context())
    interpreter.register_implementations(MissingFunctions(), override=True)
    args = tuple(
        ShapedArray(TypedPtr.new_int32(a.ravel().tolist()), list(a.shape))
        if a.dtype == np.int32
        else shaped_array(a.ravel().tolist(), list(a.shape))
        for a in arrays
    )
    results = interpreter.call_op(name, args)
    return tuple(np.array(r.data).reshape(r.shape) for r in results)
