# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

from basix import CellType, ElementFamily, LagrangeVariant
from basix.cell import topology
from basix.finite_element import tp_dof_ordering
from basix_uflx import element
from uflx import (
    Coefficient,
    TestFunction,
    coordinate_element,
    dx,
    function_space,
    inner,
)


def tp_lagrange(cell: CellType, degree: int, shape: tuple[int, ...] | None = None):
    family, variant = ElementFamily.P, LagrangeVariant.gll_warped
    dof_ordering = tp_dof_ordering(family, cell, degree, variant)
    return element(family, cell, degree, variant, shape=shape, dof_ordering=dof_ordering)


def mass_form(degree: int, cell: CellType = CellType.quadrilateral, geometry_degree: int = 1):
    gdim = len(topology(cell)) - 1
    coordinates = coordinate_element(tp_lagrange(cell, geometry_degree, shape=(gdim,)))
    V = function_space(coordinates, tp_lagrange(cell, degree))

    v = TestFunction(V)
    u = Coefficient(V)

    return inner(u, v) * dx
