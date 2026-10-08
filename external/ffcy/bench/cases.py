"""Write a case that checks and times a form's action.

A case holds the action's arguments and expected result for a few cells, which a
driver repeats to any number of cells. Precomputed geometry is evaluated by the
interpreter, so the case's arguments are those of the action alone.
"""

import argparse
from pathlib import Path

import numpy as np
from basix import CellType
from basix.cell import topology
from xdsl.dialects.builtin import ModuleOp

from bench.action import add_arguments, compile_action, function, parse
from ffcy.interpreter import run
from ffcy.reference import (
    assembled_action,
    box_mesh,
    raise_geometry,
    random_cells,
    tp_element,
    transpose_dofmap,
)
from forms import FORMS

NUM_CELLS = 4
# A box of 64 cells splits into blocks of up to 64.
BOX_CELLS = 64


def trailing_shape(t) -> str:
    return " ".join(map(str, t.get_shape()[1:]))


def write_case(prefix: Path, module: ModuleOp, args: argparse.Namespace) -> None:
    """Write the action's arguments and expected result for NUM_CELLS cells.

    prefix.case lists the shape of each argument after the cell axis, and
    prefix.bin holds the arrays in the same order.
    """
    cell = CellType[args.cell]
    rng = np.random.default_rng(0)
    x = raise_geometry(cell, random_cells(cell, NUM_CELLS, rng), args.geometry_degree)
    u = rng.uniform(-1, 1, (NUM_CELLS, tp_element(cell, args.degree).dim))
    inputs = [*run(module, "action_setup", x), u] if args.precompute_geometry else [x, u]
    y = FORMS[args.form][1](cell, x, u, args.degree, args.quadrature_degree)

    signature = function(module, "action").function_type
    lines = [f"cells {NUM_CELLS}"]
    for kind, types in (("input", signature.inputs), ("output", signature.outputs)):
        lines += [f"{kind} {trailing_shape(t)}" for t in types]
    prefix.with_suffix(".case").write_text("\n".join(lines) + "\n")
    with prefix.with_suffix(".bin").open("wb") as out:
        for a in (*inputs, y):
            np.ascontiguousarray(a, dtype=np.float64).tofile(out)


def write_assembled_case(prefix: Path, module: ModuleOp, args: argparse.Namespace) -> None:
    """Write the action on assembled vectors over a box of BOX_CELLS cells.

    prefix.case gives the number of cells and dofs, the position of each of a cell's
    dofs on its lattice, from which a driver numbers larger boxes, and the shapes of
    the geometry and of the action's values on each cell. prefix.bin holds the
    geometry, the vector, the dofmap, its transpose and the expected result.
    """
    cell = CellType[args.cell]
    d = len(topology(cell)) - 1
    rng = np.random.default_rng(0)
    n = round(BOX_CELLS ** (1 / d))
    x, dofmap, num_dofs = box_mesh(cell, args.degree, n, rng)
    x = raise_geometry(cell, x, args.geometry_degree)
    transpose = transpose_dofmap(dofmap, num_dofs, 2**d)
    v = rng.uniform(-1, 1, num_dofs)
    geometry = run(module, "action_setup", x)[0] if args.precompute_geometry else x
    y = assembled_action(
        FORMS[args.form][1], cell, x, dofmap, num_dofs, v, args.degree, args.quadrature_degree
    )

    signature = function(module, "action").function_type
    lattice = box_mesh(cell, args.degree, 1, rng)[1][0]
    size = args.degree + 1
    positions = np.stack(np.unravel_index(lattice, (size,) * d), axis=1)
    lines = [
        f"cells {len(dofmap)}",
        f"dofs {num_dofs}",
        f"lattice {' '.join(map(str, positions.ravel()))}",
        f"geometry {trailing_shape(signature.inputs.data[0])}",
        f"values {trailing_shape(signature.outputs.data[0])}",
    ]
    prefix.with_suffix(".case").write_text("\n".join(lines) + "\n")
    with prefix.with_suffix(".bin").open("wb") as out:
        for a, dtype in (
            (geometry, np.float64),
            (v, np.float64),
            (dofmap, np.int32),
            (transpose, np.int32),
            (y, np.float64),
        ):
            np.ascontiguousarray(a, dtype=dtype).tofile(out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    parser.add_argument("prefix", type=Path, help="write prefix.case and prefix.bin")
    args = parse(parser)
    write = write_assembled_case if args.assemble else write_case
    write(args.prefix, compile_action(args), args)


if __name__ == "__main__":
    main()
