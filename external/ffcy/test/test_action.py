import numpy as np
import pytest
from basix import CellType

from ffcy.dialects import uflx
from ffcy.interpreter import run
from ffcy.pipeline import compile_form
from ffcy.reference import jacobian_determinant, random_cells, tp_element
from forms import FORMS
from forms.mass import mass_form


@pytest.mark.parametrize(("form", "action"), FORMS.values(), ids=FORMS.keys())
@pytest.mark.parametrize(
    ("precompute", "sum_factorise"), [(False, False), (False, True), (True, False)]
)
@pytest.mark.parametrize("cell", [CellType.quadrilateral, CellType.hexahedron])
@pytest.mark.parametrize(
    "degree", [1, 2, 3, *(pytest.param(k, marks=pytest.mark.slow) for k in range(4, 8))]
)
def test_action(degree, cell, precompute, sum_factorise, form, action):
    rng = np.random.default_rng(degree)
    num_cells = 2
    x = random_cells(cell, num_cells, rng)
    u = rng.uniform(-1, 1, (num_cells, tp_element(cell, degree).dim))
    quadrature_degree = 2 * degree + x.shape[2] - 1

    module = compile_form(
        form(degree, cell), quadrature_degree, precompute, sum_factorise_geometry=sum_factorise
    )
    assert not any(op.dialect_name() == uflx.UFLx.name for op in module.walk())

    if precompute:
        (qdata,) = run(module, "action_setup", x)
        (y,) = run(module, "action", qdata, u)
    else:
        (y,) = run(module, "action", x, u)
    assert np.allclose(y, action(cell, x, u, degree, quadrature_degree))


@pytest.mark.parametrize("cell", [CellType.quadrilateral, CellType.hexahedron])
def test_mass_qdata(cell):
    """The mass action precomputes |det J| at each quadrature point."""
    degree, num_cells, quadrature_degree = 2, 2, 5
    x = random_cells(cell, num_cells, np.random.default_rng(0))
    module = compile_form(mass_form(degree, cell), quadrature_degree, precompute_geometry=True)
    (qdata,) = run(module, "action_setup", x)
    expected = np.abs(jacobian_determinant(cell, x, quadrature_degree))
    assert np.allclose(qdata.reshape(num_cells, -1), expected)
