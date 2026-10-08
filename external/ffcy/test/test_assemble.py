# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

import numpy as np
import pytest
from basix import CellType
from basix.cell import topology

from ffcy.interpreter import run
from ffcy.pipeline import compile_form
from ffcy.reference import assembled_action, box_mesh, transpose_dofmap
from forms import FORMS


@pytest.mark.parametrize(("form", "action"), FORMS.values(), ids=FORMS.keys())
@pytest.mark.parametrize("precompute", [False, True])
@pytest.mark.parametrize("cell", [CellType.quadrilateral, CellType.hexahedron])
@pytest.mark.parametrize("degree", [1, 2])
@pytest.mark.parametrize("scatter", ["transpose", "none", "atomic"])
def test_assembled_action(scatter, degree, cell, precompute, form, action):
    """y = G^T A G x on a box mesh, with G gathering each cell's dofs."""
    rng = np.random.default_rng(degree)
    d = len(topology(cell)) - 1
    x, dofmap, num_dofs = box_mesh(cell, degree, 2, rng)
    transpose = transpose_dofmap(dofmap, num_dofs, 2**d)
    v = rng.uniform(-1, 1, num_dofs)
    quadrature_degree = 2 * degree + 1

    module = compile_form(
        form(degree, cell), quadrature_degree, precompute, assemble=True, scatter=scatter
    )
    geometry = run(module, "action_setup", x)[0] if precompute else x
    if scatter == "transpose":
        (y,) = run(module, "action", geometry, v, dofmap, transpose)
    else:
        (values,) = run(module, "action", geometry, v, dofmap)
        y = np.zeros(num_dofs)
        np.add.at(y, dofmap, values.reshape(dofmap.shape))
    expected = assembled_action(action, cell, x, dofmap, num_dofs, v, degree, quadrature_degree)
    assert np.allclose(y, expected)
