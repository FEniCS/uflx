# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

import subprocess
from collections.abc import Sequence
from pathlib import Path

from uflx.integrals import AbstractIntegral
from xdsl.context import Context
from xdsl.dialects.builtin import ModuleOp
from xdsl.passes import ModulePass
from xdsl.transforms.common_subexpression_elimination import CommonSubexpressionElimination
from xdsl.transforms.dead_code_elimination import DeadCodeElimination
from xdsl.transforms.mlir_opt import DEFAULT_MLIR_OPT_EXECUTABLE

from ffcy.frontend import from_uflx
from ffcy.transforms.assemble import Assemble
from ffcy.transforms.convert_uflx_to_linalg import ConvertUFLxToLinalg
from ffcy.transforms.expand_results import ExpandResults
from ffcy.transforms.precompute_geometry import PrecomputeGeometry


def apply(module: ModuleOp, *passes: ModulePass) -> ModuleOp:
    ctx = Context()
    for p in passes:
        p.apply(ctx, module)
    module.verify()
    return module


# Returns text because xDSL cannot parse every dialect a pipeline lowers to.
def run_mlir_opt(
    module: ModuleOp | str,
    pipeline: str,
    plugins: Sequence[str | Path] = (),
    executable: str = DEFAULT_MLIR_OPT_EXECUTABLE,
) -> str:
    command = [executable, *(f"--load-pass-plugin={p}" for p in plugins)]
    command.append(f"--pass-pipeline=builtin.module({pipeline})")
    result = subprocess.run(command, input=str(module), capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"mlir-opt failed with {pipeline}\n{result.stderr}")
    return result.stdout


# Lowers a form's integral to linalg on tensors with a leading cell axis, which a
# target's layout pass then blocks, such as SplitCells or InterleaveCells.
def compile_form(
    integral: AbstractIntegral,
    quadrature_degree: int,
    precompute_geometry: bool = False,
    sum_factorise_geometry: bool = False,
    expand_results: bool = False,
    assemble: bool = False,
    scatter: str = "transpose",
) -> ModuleOp:
    module = from_uflx(integral, quadrature_degree)
    passes: list[ModulePass] = [PrecomputeGeometry()] if precompute_geometry else []
    if assemble:
        passes += [Assemble(scatter)]
    passes += [ConvertUFLxToLinalg(sum_factorise_geometry), CommonSubexpressionElimination()]
    if expand_results:
        passes += [ExpandResults(), DeadCodeElimination()]
    return apply(module, *passes)
