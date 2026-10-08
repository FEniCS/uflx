"""Emit a form's action as MLIR for a target's pipeline.

The action's cells are split into blocks for GPUs or interleaved for CPUs, its ops
are grouped into the stages of a pencil schedule, and its divisions may multiply by
a reciprocal, as both targets' pipelines expect. Precomputed geometry is evaluated
by bench/cases.py, so only the action is emitted.
"""

import argparse
from pathlib import Path

from xdsl.dialects import func

from bench.action import add_arguments, compile_action, parse
from ffcy.pipeline import apply
from ffcy.transforms.allow_reciprocal import AllowReciprocal
from ffcy.transforms.annotate_stages import AnnotateStages
from ffcy.transforms.interleave_cells import InterleaveCells
from ffcy.transforms.split_cells import SplitCells


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    layout = parser.add_mutually_exclusive_group()
    layout.add_argument(
        "--cells-per-block", type=int, help="split the cell axis into blocks of this size"
    )
    layout.add_argument(
        "--lanes", type=int, help="interleave cells in blocks of this size, for CPUs"
    )
    parser.add_argument("-o", "--output", type=Path, help="write the action as MLIR")
    args = parse(parser)

    module = compile_action(args)
    # The setup is evaluated for cases by bench/cases.py.
    for f in list(module.ops):
        if isinstance(f, func.FuncOp) and f.sym_name.data.endswith("_setup"):
            module.body.block.erase_op(f)
    layout = []
    if args.cells_per_block:
        layout = [SplitCells(args.cells_per_block)]
    elif args.lanes:
        layout = [InterleaveCells(args.lanes)]
    apply(module, *layout, AnnotateStages(), AllowReciprocal())
    text = str(module) + "\n"
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
