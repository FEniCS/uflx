# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

"""Check and time every benchmark built in a CMake build directory."""

import argparse
import re
import subprocess
from pathlib import Path

NAME = re.compile(
    r"lib(?P<case>(?P<form>[a-z]+)_[a-z]+_q(?P<degree>\d+)_(?P<geometry>[a-z0-9_]+))"
    r"_c(?P<cells>\d+)\.so$"
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path, help="CMake build directory")
    parser.add_argument("--filter", default="", help="only names containing this")
    parser.add_argument("--seconds", type=float, default=2.0)
    parser.add_argument("--num-dofs", type=float, default=4e6)
    args = parser.parse_args()

    runs = []
    for library in args.build.iterdir():
        match = NAME.match(library.name)
        if match and args.filter in library.name:
            key = (match["form"], int(match["degree"]), match["geometry"], int(match["cells"]))
            runs.append((key, library.name, match["case"]))
    for (form, degree, geometry, cells), name, case in sorted(runs):
        command = [args.build / "ffcy_driver", args.build / name, args.build / "cases" / case]
        command += ["--cells-per-block", str(cells), "--seconds", str(args.seconds)]
        command += ["--num-dofs", str(args.num_dofs)]
        # A fused action scatters its own values.
        if geometry.endswith("_fused"):
            command += ["--scatter", "fused"]
        run = subprocess.run(command, capture_output=True, text=True)
        if run.returncode:
            result = "failed: " + run.stderr.strip().splitlines()[-1]
        else:
            error, timing = run.stdout.strip().split(", ", 1)
            # A kernel that fails to launch leaves NaN in y, which only the check notices.
            result = timing if float(error.split()[1]) < 1e-12 else f"failed check, {error}"
        print(f"{form:9s} Q{degree} {geometry:21s} c{cells:<3d} {result}", flush=True)


if __name__ == "__main__":
    main()
