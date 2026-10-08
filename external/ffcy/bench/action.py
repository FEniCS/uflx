# The action of a form chosen on the command line, for bench/emit.py and bench/cases.py.

import argparse

from basix import CellType
from xdsl.dialects import func
from xdsl.dialects.builtin import ModuleOp

from ffcy.pipeline import compile_form
from forms import FORMS


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--form", choices=sorted(FORMS), default="mass")
    parser.add_argument("--cell", choices=["quadrilateral", "hexahedron"], default="hexahedron")
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument(
        "--quadrature-degree", type=int, help="defaults to 2k + 1, as many points as dofs"
    )
    parser.add_argument(
        "--geometry-degree",
        type=int,
        default=1,
        help="describe the cells by geometry of this degree, with straight sides",
    )
    parser.add_argument("--sum-factorise-geometry", action="store_true")
    parser.add_argument("--precompute-geometry", action="store_true")
    parser.add_argument(
        "--assemble",
        action="store_true",
        help="act on assembled vectors through a dofmap, leaving the scatter to the driver",
    )
    parser.add_argument(
        "--atomic-scatter",
        action="store_true",
        help="with --assemble, scatter in the action with atomic adds",
    )


def parse(parser: argparse.ArgumentParser) -> argparse.Namespace:
    args = parser.parse_args()
    args.quadrature_degree = args.quadrature_degree or 2 * args.degree + 1
    return args


# Results keep their tensor product shape, so they can be written in place.
def compile_action(args: argparse.Namespace) -> ModuleOp:
    form = FORMS[args.form][0]
    return compile_form(
        form(args.degree, CellType[args.cell], geometry_degree=args.geometry_degree),
        args.quadrature_degree,
        args.precompute_geometry,
        sum_factorise_geometry=args.sum_factorise_geometry,
        expand_results=True,
        assemble=args.assemble,
        scatter="atomic" if args.atomic_scatter else "none",
    )


def function(module: ModuleOp, name: str) -> func.FuncOp:
    return next(f for f in module.ops if isinstance(f, func.FuncOp) and f.sym_name.data == name)
