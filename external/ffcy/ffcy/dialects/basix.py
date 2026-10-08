# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

from collections.abc import Sequence

import basix
import numpy as np
import numpy.typing as npt
from xdsl.dialects.builtin import IntAttr
from xdsl.dialects.utils import EnumAttribute
from xdsl.ir import Attribute, Dialect, ParametrizedAttribute
from xdsl.irdl import irdl_attr_definition
from xdsl.parser import AttrParser
from xdsl.printer import Printer
from xdsl.utils.str_enum import StrEnum


# Members mirror the Basix enums of the same name, listing only what is supported.
class CellType(StrEnum):
    interval = "interval"
    quadrilateral = "quadrilateral"
    hexahedron = "hexahedron"


def topological_dimension(cell_type: CellType) -> int:
    return len(basix.cell.topology(basix.CellType[cell_type])) - 1


class ElementFamily(StrEnum):
    P = "P"


class LagrangeVariant(StrEnum):
    gll_warped = "gll_warped"


class QuadratureType(StrEnum):
    gauss_jacobi = "gauss_jacobi"
    gll = "gll"


@irdl_attr_definition
class CellTypeAttr(EnumAttribute[CellType]):
    name = "basix.cell_type"


@irdl_attr_definition
class ElementFamilyAttr(EnumAttribute[ElementFamily]):
    name = "basix.element_family"


@irdl_attr_definition
class LagrangeVariantAttr(EnumAttribute[LagrangeVariant]):
    name = "basix.lagrange_variant"


@irdl_attr_definition
class QuadratureTypeAttr(EnumAttribute[QuadratureType]):
    name = "basix.quadrature_type"


@irdl_attr_definition
class ElementAttr(ParametrizedAttribute):
    name = "basix.element"

    family: ElementFamilyAttr
    cell_type: CellTypeAttr
    degree: IntAttr
    lagrange_variant: LagrangeVariantAttr

    def __init__(
        self,
        family: ElementFamily,
        cell_type: CellType,
        degree: int,
        lagrange_variant: LagrangeVariant,
    ):
        super().__init__(
            ElementFamilyAttr(family),
            CellTypeAttr(cell_type),
            IntAttr(degree),
            LagrangeVariantAttr(lagrange_variant),
        )

    @classmethod
    def from_basix(cls, element: basix.finite_element.FiniteElement) -> "ElementAttr":
        return cls(
            ElementFamily(element.family.name),
            CellType(element.cell_type.name),
            element.degree,
            LagrangeVariant(element.lagrange_variant.name),
        )

    def tensor_factors(self) -> list[basix.finite_element.FiniteElement]:
        element = basix.create_tp_element(
            basix.ElementFamily[self.family.data],
            basix.CellType[self.cell_type.data],
            self.degree.data,
            basix.LagrangeVariant[self.lagrange_variant.data],
        )
        return element.get_tensor_product_representation()[0]

    @classmethod
    def parse_parameters(cls, parser: AttrParser) -> Sequence[Attribute]:
        with parser.in_angle_brackets():
            family = parser.parse_str_enum(ElementFamily)
            parser.parse_punctuation(",")
            cell_type = parser.parse_str_enum(CellType)
            parser.parse_punctuation(",")
            degree = parser.parse_integer()
            parser.parse_punctuation(",")
            lagrange_variant = parser.parse_str_enum(LagrangeVariant)
        return cls(family, cell_type, degree, lagrange_variant).parameters

    def print_parameters(self, printer: Printer) -> None:
        printer.print_string(
            f"<{self.family.data}, {self.cell_type.data}, {self.degree.data}, "
            f"{self.lagrange_variant.data}>"
        )


@irdl_attr_definition
class QuadratureAttr(ParametrizedAttribute):
    name = "basix.quadrature"

    cell_type: CellTypeAttr
    degree: IntAttr
    type: QuadratureTypeAttr

    def __init__(
        self,
        cell_type: CellType,
        degree: int,
        type: QuadratureType = QuadratureType.gauss_jacobi,
    ):
        super().__init__(CellTypeAttr(cell_type), IntAttr(degree), QuadratureTypeAttr(type))

    def tensor_factor(self) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        points, weights = basix.make_quadrature(
            basix.CellType.interval,
            self.degree.data,
            rule=basix.QuadratureType[self.type.data],
        )
        return points[:, 0], weights

    @classmethod
    def parse_parameters(cls, parser: AttrParser) -> Sequence[Attribute]:
        with parser.in_angle_brackets():
            cell_type = parser.parse_str_enum(CellType)
            parser.parse_punctuation(",")
            degree = parser.parse_integer()
            parser.parse_punctuation(",")
            type = parser.parse_str_enum(QuadratureType)
        return cls(cell_type, degree, type).parameters

    def print_parameters(self, printer: Printer) -> None:
        printer.print_string(f"<{self.cell_type.data}, {self.degree.data}, {self.type.data}>")


Basix = Dialect(
    "basix",
    [],
    [
        CellTypeAttr,
        ElementFamilyAttr,
        LagrangeVariantAttr,
        QuadratureTypeAttr,
        ElementAttr,
        QuadratureAttr,
    ],
)
