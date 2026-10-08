import numpy as np
import pytest
from basix import CellType

from ffcy.interpreter import run
from ffcy.pipeline import apply, compile_form
from ffcy.reference import random_cells, tp_element
from ffcy.transforms.interleave_cells import InterleaveCells
from forms import FORMS


@pytest.mark.parametrize(("form", "action"), FORMS.values(), ids=FORMS.keys())
@pytest.mark.parametrize("precompute", [False, True])
@pytest.mark.parametrize("cell", [CellType.quadrilateral, CellType.hexahedron])
def test_interleave_cells(form, action, precompute, cell):
    """Interleaving cells in blocks leaves the action unchanged."""
    degree, lanes, num_blocks = 2, 4, 2
    rng = np.random.default_rng(0)
    x = random_cells(cell, lanes * num_blocks, rng)
    u = rng.uniform(-1, 1, (len(x), tp_element(cell, degree).dim))
    quadrature_degree = 2 * degree + 2

    module = apply(
        compile_form(form(degree, cell), quadrature_degree, precompute, expand_results=True),
        InterleaveCells(lanes),
    )

    def interleaved(a):
        return np.moveaxis(a.reshape(num_blocks, lanes, *a.shape[1:]), 1, -1).copy()

    if precompute:
        (qdata,) = run(module, "action_setup", interleaved(x))
        (y,) = run(module, "action", qdata, interleaved(u))
    else:
        (y,) = run(module, "action", interleaved(x), interleaved(u))
    y = np.moveaxis(y, -1, 1).reshape(u.shape)
    assert np.allclose(y, action(cell, x, u, degree, quadrature_degree))
