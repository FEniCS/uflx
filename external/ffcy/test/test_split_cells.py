import numpy as np
import pytest
from basix import CellType

from ffcy.interpreter import run
from ffcy.pipeline import apply, compile_form
from ffcy.reference import random_cells, tp_element
from ffcy.transforms.split_cells import SplitCells
from forms import FORMS


@pytest.mark.parametrize(("form", "action"), FORMS.values(), ids=FORMS.keys())
@pytest.mark.parametrize("precompute", [False, True])
def test_split_cells(form, action, precompute):
    """Splitting the cell axis into blocks leaves the action unchanged."""
    cell, degree, cells_per_block, num_blocks = CellType.hexahedron, 2, 2, 2
    rng = np.random.default_rng(0)
    x = random_cells(cell, cells_per_block * num_blocks, rng)
    u = rng.uniform(-1, 1, (len(x), tp_element(cell, degree).dim))
    quadrature_degree = 2 * degree + 2

    module = apply(
        compile_form(form(degree, cell), quadrature_degree, precompute, expand_results=True),
        SplitCells(cells_per_block),
    )

    def blocked(a):
        return a.reshape(num_blocks, cells_per_block, *a.shape[1:])

    if precompute:
        (qdata,) = run(module, "action_setup", blocked(x))
        (y,) = run(module, "action", qdata, blocked(u))
    else:
        (y,) = run(module, "action", blocked(x), blocked(u))
    assert np.allclose(y.reshape(u.shape), action(cell, x, u, degree, quadrature_degree))
