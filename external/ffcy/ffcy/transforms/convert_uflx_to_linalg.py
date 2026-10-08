# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from functools import reduce
from itertools import product

import numpy as np
import numpy.typing as npt
from xdsl.builder import ImplicitBuilder
from xdsl.context import Context
from xdsl.dialects import arith, math, tensor
from xdsl.dialects.builtin import (
    DYNAMIC_INDEX,
    FloatAttr,
    IndexType,
    IntegerAttr,
    ModuleOp,
    TensorType,
    f64,
    i32,
)
from xdsl.ir import Operation, SSAValue
from xdsl.ir.affine import AffineExpr, AffineMap
from xdsl.passes import ModulePass
from xdsl.pattern_rewriter import (
    GreedyRewritePatternApplier,
    PatternRewriter,
    PatternRewriteWalker,
    RewritePattern,
    op_type_rewrite_pattern,
)

from ffcy.dialects import uflx
from ffcy.transforms.linalg_builders import (
    CONTRACT,
    add_all,
    collapse,
    constant,
    contract,
    contract_axes,
    dims,
    empty,
    expand,
    generic,
    multiply_all,
    num_cells_of,
    reassociation_attr,
)


def tables(
    op: uflx.TensorProductOp, order: int
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Derivatives up to `order` of the 1D basis at the 1D points, and the 1D weights.

    Table n has shape [points, dofs] and holds the nth derivatives.
    """
    factor = op.element.tensor_factors()[0]
    points, weights = op.quadrature.tensor_factor()
    return factor.tabulate(order, points[:, None])[:, :, :, 0], weights


def lower_coefficient(op: uflx.CoefficientOp) -> SSAValue:
    derivative = op.derivative.get_values()
    T, _ = tables(op, max(derivative))
    num_cells = num_cells_of(op.input)
    u = expand(num_cells, op.input, [[T.shape[2]] * len(derivative)])
    return contract_axes(num_cells, [T[a] for a in derivative], u)


def lower_jacobian_sum_factorised(op: uflx.JacobianOp) -> SSAValue:
    g, r = op.component.get_values()
    T, _ = tables(op, 1)
    d = uflx.tdim(op.element)
    num_cells = num_cells_of(op.input)
    x = expand(num_cells, op.input, [[T.shape[2]] * d, [d]])
    # The derivative along r of every coordinate, with shape [cells, nq, ..., nq, gdim].
    derivative = contract_axes(num_cells, [T[int(axis == r)] for axis in range(d)], x)
    return generic(
        [derivative],
        [
            AffineMap(d + 1, 0, (*dims(*range(d + 1)), AffineExpr.constant(g))),
            AffineMap.identity(d + 1),
        ],
        empty(num_cells, [T.shape[1]] * d),
        lambda value, _: value,
    )


def lower_jacobian_multilinear(op: uflx.JacobianOp) -> SSAValue:
    """Evaluate a component of J pointwise from the cell's vertices, for degree 1 geometry."""
    g, r = op.component.get_values()
    T, _ = tables(op, 1)
    B, D = T
    d = uflx.tdim(op.element)
    num_cells = num_cells_of(op.input)
    x = expand(num_cells, op.input, [[2] * d, [d]])
    vertices = list(product((0, 1), repeat=d))
    # Basix tabulates the derivatives of the 1D Q1 basis as constants in the point.
    dB = D[0]

    def vertex_map(v: tuple[int, ...]) -> AffineMap:
        constants = (AffineExpr.constant(i) for i in (*v, g))
        return AffineMap(d + 1, 0, (AffineExpr.dimension(0), *constants))

    def table_map(s: int, a: int) -> AffineMap:
        return AffineMap(d + 1, 0, (AffineExpr.dimension(s + 1), AffineExpr.constant(a)))

    def body(*args: SSAValue) -> SSAValue:
        corners, values = args[: len(vertices)], args[len(vertices) :]
        X = dict(zip(vertices, corners, strict=True))
        L = [values[2 * s : 2 * s + 2] for s in range(d)]
        dl = [arith.ConstantOp(FloatAttr(float(c), f64)).result for c in dB]
        others = [s for s in range(d) if s != r]
        terms = []
        # Along axis r the coordinate is linear, so J interpolates the edge differences.
        for rest in product((0, 1), repeat=d - 1):
            weight = [L[s][a] for s, a in zip(others, rest, strict=True)]
            v0, v1 = (tuple([*rest[:r], a, *rest[r:]]) for a in (0, 1))
            edge = add_all(
                [arith.MulfOp(dl[a], X[v], CONTRACT).result for a, v in enumerate((v0, v1))]
            )
            terms.append(multiply_all([*weight, edge]))
        return add_all(terms)

    return generic(
        [x] * len(vertices) + [constant(B)] * (2 * d),
        [vertex_map(v) for v in vertices]
        + [table_map(s, a) for s in range(d) for a in (0, 1)]
        + [AffineMap.identity(d + 1)],
        empty(num_cells, [B.shape[0]] * d),
        lambda *args: body(*args[:-1]),
    )


def lower_jacobian(op: uflx.JacobianOp, sum_factorise: bool = False) -> SSAValue:
    if op.element.degree.data == 1 and not sum_factorise:
        return lower_jacobian_multilinear(op)
    return lower_jacobian_sum_factorised(op)


def integrate(
    num_cells: SSAValue, f: SSAValue, T: npt.NDArray, w: npt.NDArray, derivative: Sequence[int]
) -> SSAValue:
    """One term of an integral, with shape [cells, n, ..., n]."""
    d = len(derivative)
    weighted = generic(
        [f, constant(reduce(np.multiply.outer, [w] * d))],
        [
            AffineMap.identity(d + 1),
            AffineMap(d + 1, 0, dims(*range(1, d + 1))),
            AffineMap.identity(d + 1),
        ],
        empty(num_cells, [T.shape[1]] * d),
        lambda fq, wq, _: arith.MulfOp(fq, wq, CONTRACT).result,
    )
    # Starting from the last axis, where evaluation ends, lets pencils line up.
    for axis in reversed(range(1, d + 1)):
        weighted = contract(num_cells, T[derivative[axis - 1]].T, weighted, axis)
    return weighted


def lower_integral(op: uflx.IntegralOp) -> SSAValue:
    derivatives = [list(a.get_values()) for a in op.derivatives]
    T, w = tables(op, max(max(a) for a in derivatives))
    num_cells = num_cells_of(op.integrands[0])
    terms = [
        integrate(num_cells, f, T, w, a) for f, a in zip(op.integrands, derivatives, strict=True)
    ]
    rank = len(terms[0].type.get_shape())
    total = reduce(
        lambda a, b: generic(
            [a, b],
            [AffineMap.identity(rank)] * 3,
            empty(num_cells, a.type.get_shape()[1:]),
            lambda x, y, _: arith.AddfOp(x, y, CONTRACT).result,
            # Integration ends on the first axis.
            pencil=1,
        ),
        terms,
    )
    return collapse(total)


def lower_gather(op: uflx.GatherOp) -> SSAValue:
    """Gather in tensor product shape, along the axis the first contraction loops over."""
    n, d = op.element.degree.data + 1, uflx.tdim(op.element)
    num_cells = num_cells_of(op.dofmap)
    dofmap = expand(num_cells, op.dofmap, [[n] * d])

    def body(dof: SSAValue, _: SSAValue) -> SSAValue:
        index = arith.IndexCastOp(dof, IndexType()).result
        return tensor.ExtractOp(op.vector, [index], f64).result

    values = generic(
        [dofmap], [AffineMap.identity(d + 1)] * 2, empty(num_cells, [n] * d), body, pencil=1
    )
    return collapse(values)


def lower_scatter_add(op: uflx.ScatterAddOp) -> SSAValue:
    """Sum each dof's entries of the flattened input, skipping the -1 padding.

    The number of entries is fixed, so the sum is pointwise over dofs, which needs no
    zeroed output.
    """
    flat = tensor.CollapseShapeOp.build(
        operands=[op.input],
        result_types=[TensorType(f64, [DYNAMIC_INDEX])],
        properties={"reassociation": reassociation_attr([[0, 1]])},
    ).result
    width = op.transpose.type.get_shape()[1]

    def entry(e: SSAValue) -> SSAValue:
        zero = arith.ConstantOp(IntegerAttr(0, i32)).result
        valid = arith.CmpiOp(e, zero, "sge").result
        index = arith.IndexCastOp(arith.SelectOp(valid, e, zero), IndexType()).result
        value = tensor.ExtractOp(flat, [index], f64).result
        return arith.SelectOp(valid, value, arith.ConstantOp(FloatAttr(0.0, f64))).result

    return generic(
        [op.transpose] * width,
        [AffineMap(1, 0, (AffineExpr.dimension(0), AffineExpr.constant(k))) for k in range(width)]
        + [AffineMap.identity(1)],
        empty(num_cells_of(op.transpose), []),
        lambda *args: add_all([entry(e) for e in args[:-1]]),
    )


LOWERINGS: dict[type[uflx.TensorProductOp], Callable[..., SSAValue]] = {
    uflx.CoefficientOp: lower_coefficient,
    uflx.IntegralOp: lower_integral,
}


class LowerRestriction(RewritePattern):
    """Lower gathers onto cells and scatters back to assembled vectors."""

    def match_and_rewrite(self, op: Operation, rewriter: PatternRewriter, /):
        lowering = {uflx.GatherOp: lower_gather, uflx.ScatterAddOp: lower_scatter_add}.get(type(op))
        if lowering is None:
            return
        with ImplicitBuilder(rewriter):
            result = lowering(op)
        rewriter.replace(op, [], [result])


@dataclass(frozen=True)
class LowerTensorProductOp(RewritePattern):
    sum_factorise_geometry: bool = False

    @op_type_rewrite_pattern
    def match_and_rewrite(self, op: uflx.TensorProductOp, rewriter: PatternRewriter, /):
        with ImplicitBuilder(rewriter):
            if isinstance(op, uflx.JacobianOp):
                result = lower_jacobian(op, self.sum_factorise_geometry)
            else:
                result = LOWERINGS[type(op)](op)
        rewriter.replace(op, [], [result])


ELEMENTWISE: dict[type[Operation], Callable[..., SSAValue]] = {
    math.AbsFOp: lambda x: math.AbsFOp(x).result,
    arith.NegfOp: lambda x: arith.NegfOp(x).result,
    arith.AddfOp: lambda a, b: arith.AddfOp(a, b, CONTRACT).result,
    arith.SubfOp: lambda a, b: arith.SubfOp(a, b, CONTRACT).result,
    arith.MulfOp: lambda a, b: arith.MulfOp(a, b, CONTRACT).result,
    arith.DivfOp: lambda a, b: arith.DivfOp(a, b, CONTRACT).result,
}


class LowerElementwise(RewritePattern):
    """Lower elementwise ops on tensors to linalg.generic with a fresh output."""

    def match_and_rewrite(self, op: Operation, rewriter: PatternRewriter, /):
        body = ELEMENTWISE.get(type(op))
        if body is None or not isinstance(op.results[0].type, TensorType):
            return
        shape = op.results[0].type.get_shape()
        with ImplicitBuilder(rewriter):
            result = generic(
                op.operands,
                [AffineMap.identity(len(shape))] * (len(op.operands) + 1),
                empty(num_cells_of(op.operands[0]), shape[1:]),
                lambda *args: body(*args[:-1]),
            )
        rewriter.replace(op, [], [result])


@dataclass(frozen=True)
class ConvertUFLxToLinalg(ModulePass):
    name = "convert-uflx-to-linalg"

    # Use contractions for degree 1 geometry too, rather than the pointwise form.
    sum_factorise_geometry: bool = False

    def apply(self, ctx: Context, op: ModuleOp) -> None:
        lower = LowerTensorProductOp(self.sum_factorise_geometry)
        PatternRewriteWalker(
            GreedyRewritePatternApplier([lower, LowerRestriction(), LowerElementwise()]),
            apply_recursively=False,
        ).rewrite_module(op)
