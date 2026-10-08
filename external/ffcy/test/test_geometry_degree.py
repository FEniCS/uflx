# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

import numpy as np
import pytest
from basix import CellType

from ffcy.interpreter import run
from ffcy.pipeline import compile_form
from ffcy.reference import raise_geometry, random_cells, tp_element
from forms import FORMS


@pytest.mark.parametrize(("form", "action"), FORMS.values(), ids=FORMS.keys())
@pytest.mark.parametrize("precompute", [False, True])
@pytest.mark.parametrize("cell", [CellType.quadrilateral, CellType.hexahedron])
@pytest.mark.parametrize("geometry_degree", [2, 3])
def test_geometry_degree(geometry_degree, cell, precompute, form, action):
    """Cells described by higher degree geometry have the action of their Q1 cells."""
    degree, num_cells = 2, 2
    rng = np.random.default_rng(0)
    x = random_cells(cell, num_cells, rng)
    xk = raise_geometry(cell, x, geometry_degree)
    u = rng.uniform(-1, 1, (num_cells, tp_element(cell, degree).dim))
    quadrature_degree = 2 * degree + 1

    module = compile_form(
        form(degree, cell, geometry_degree=geometry_degree), quadrature_degree, precompute
    )
    geometry = run(module, "action_setup", xk)[0] if precompute else xk
    (y,) = run(module, "action", geometry, u)
    expected = action(cell, x, u, degree, quadrature_degree)
    assert np.allclose(y.reshape(expected.shape), expected)
    assert np.allclose(action(cell, xk, u, degree, quadrature_degree), expected)
