from basix import CellType
from basix.cell import topology
from uflx import (
    Coefficient,
    TestFunction,
    coordinate_element,
    dx,
    function_space,
    grad,
    inner,
)

from forms.mass import tp_lagrange


def laplace_form(degree: int, cell: CellType = CellType.quadrilateral, geometry_degree: int = 1):
    gdim = len(topology(cell)) - 1
    coordinates = coordinate_element(tp_lagrange(cell, geometry_degree, shape=(gdim,)))
    V = function_space(coordinates, tp_lagrange(cell, degree))

    v = TestFunction(V)
    u = Coefficient(V)

    return inner(grad(u), grad(v)) * dx
