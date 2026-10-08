"""UFLx operations on a batch of quadrilateral or hexahedral cells.

E-vectors are in Basix tensor product DOF ordering and values at quadrature
points have shape [cells, nq, ..., nq] in the order of `basix.make_quadrature`.
Derivatives are multi-indices of reference derivative orders, one per axis, as
in Basix.
"""

from collections.abc import Sequence

from xdsl.dialects.builtin import (
    ArrayAttr,
    DenseArrayBase,
    Float64Type,
    IntegerType,
    TensorType,
    f64,
    i64,
)
from xdsl.ir import Dialect, Operation, SSAValue
from xdsl.irdl import (
    IRDLOperation,
    ParsePropInAttrDict,
    irdl_op_definition,
    operand_def,
    prop_def,
    result_def,
    traits_def,
    var_operand_def,
)
from xdsl.traits import Pure
from xdsl.utils.exceptions import VerifyException

from ffcy.dialects.basix import CellType, ElementAttr, QuadratureAttr, topological_dimension


def tdim(element: ElementAttr) -> int:
    return topological_dimension(element.cell_type.data)


def num_dofs(element: ElementAttr) -> int:
    return (element.degree.data + 1) ** tdim(element)


def points_shape(quadrature: QuadratureAttr) -> tuple[int, ...]:
    _, weights = quadrature.tensor_factor()
    return (len(weights),) * topological_dimension(quadrature.cell_type.data)


def values_type(like: SSAValue | Operation, quadrature: QuadratureAttr) -> TensorType:
    num_cells = SSAValue.get(like).type.get_shape()[0]
    return TensorType(f64, [num_cells, *points_shape(quadrature)])


def check_shape(op: IRDLOperation, value: SSAValue, shape: Sequence[int], name: str) -> None:
    if tuple(value.type.get_shape()[1:]) != tuple(shape):
        raise VerifyException(f"{op.name} expects {name} of shape [cells, {list(shape)}]")


class TensorProductOp(IRDLOperation):
    element = prop_def(ElementAttr)
    quadrature = prop_def(QuadratureAttr)

    irdl_options = (ParsePropInAttrDict(),)
    traits = traits_def(Pure())

    def verify_(self) -> None:
        cell_type = self.element.cell_type.data
        if cell_type not in (CellType.quadrilateral, CellType.hexahedron):
            raise VerifyException(f"{self.name} does not support {cell_type} cells")
        if self.quadrature.cell_type.data != cell_type:
            raise VerifyException(f"{self.name} element and quadrature cells differ")
        cells = {v.type.get_shape()[0] for v in (*self.operands, *self.results)}
        if len(cells) != 1:
            raise VerifyException(f"{self.name} must preserve the number of cells")


@irdl_op_definition
class CoefficientOp(TensorProductOp):
    """A reference derivative of a coefficient at the quadrature points."""

    name = "uflx.coefficient"

    input = operand_def(TensorType[Float64Type])
    output = result_def(TensorType[Float64Type])
    derivative = prop_def(DenseArrayBase)

    assembly_format = "$input attr-dict `:` type($input) `->` type($output)"

    def __init__(
        self,
        input: SSAValue | Operation,
        element: ElementAttr,
        quadrature: QuadratureAttr,
        derivative: Sequence[int],
    ):
        super().__init__(
            operands=[input],
            result_types=[values_type(input, quadrature)],
            properties={
                "element": element,
                "quadrature": quadrature,
                "derivative": DenseArrayBase.from_list(i64, derivative),
            },
        )

    def verify_(self) -> None:
        super().verify_()
        if len(self.derivative.get_values()) != tdim(self.element):
            raise VerifyException(f"{self.name} needs one derivative order per axis")
        check_shape(self, self.input, [num_dofs(self.element)], "dofs")
        check_shape(self, self.output, points_shape(self.quadrature), "values")


@irdl_op_definition
class JacobianOp(TensorProductOp):
    """Component [g, r] of the Jacobian, the derivative of x_g along reference axis r."""

    name = "uflx.jacobian"

    input = operand_def(TensorType[Float64Type])
    output = result_def(TensorType[Float64Type])
    component = prop_def(DenseArrayBase)

    assembly_format = "$input attr-dict `:` type($input) `->` type($output)"

    def __init__(
        self,
        coordinate_dofs: SSAValue | Operation,
        element: ElementAttr,
        quadrature: QuadratureAttr,
        component: Sequence[int],
    ):
        super().__init__(
            operands=[coordinate_dofs],
            result_types=[values_type(coordinate_dofs, quadrature)],
            properties={
                "element": element,
                "quadrature": quadrature,
                "component": DenseArrayBase.from_list(i64, component),
            },
        )

    def verify_(self) -> None:
        super().verify_()
        d = tdim(self.element)
        if not all(0 <= i < d for i in self.component.get_values()):
            raise VerifyException(f"{self.name} component out of range")
        check_shape(self, self.input, [num_dofs(self.element), d], "coordinate dofs")
        check_shape(self, self.output, points_shape(self.quadrature), "values")


@irdl_op_definition
class IntegralOp(TensorProductOp):
    """Integrate terms against reference derivatives of the test functions.

    dofs[c, i] = sum_t sum_q w_q D_t phi_i(X_q) integrands[t][c, q]
    """

    name = "uflx.integral"

    integrands = var_operand_def(TensorType[Float64Type])
    dofs = result_def(TensorType[Float64Type])
    derivatives = prop_def(ArrayAttr[DenseArrayBase])

    assembly_format = "$integrands attr-dict `:` type($integrands) `->` type($dofs)"

    def __init__(
        self,
        integrands: Sequence[SSAValue],
        element: ElementAttr,
        quadrature: QuadratureAttr,
        derivatives: Sequence[Sequence[int]],
    ):
        num_cells = integrands[0].type.get_shape()[0]
        super().__init__(
            operands=[integrands],
            result_types=[TensorType(f64, [num_cells, num_dofs(element)])],
            properties={
                "element": element,
                "quadrature": quadrature,
                "derivatives": ArrayAttr(DenseArrayBase.from_list(i64, d) for d in derivatives),
            },
        )

    def verify_(self) -> None:
        super().verify_()
        if len(self.derivatives) != len(self.integrands):
            raise VerifyException(f"{self.name} needs one derivative per integrand")
        for integrand in self.integrands:
            check_shape(self, integrand, points_shape(self.quadrature), "integrands")
        check_shape(self, self.dofs, [num_dofs(self.element)], "dofs")


@irdl_op_definition
class GatherOp(IRDLOperation):
    """The values of an assembled vector on each cell, output[c, i] = vector[dofmap[c, i]]."""

    name = "uflx.gather"

    vector = operand_def(TensorType[Float64Type])
    dofmap = operand_def(TensorType[IntegerType])
    output = result_def(TensorType[Float64Type])
    element = prop_def(ElementAttr)

    irdl_options = (ParsePropInAttrDict(),)
    traits = traits_def(Pure())

    assembly_format = (
        "$vector `,` $dofmap attr-dict `:` type($vector) `,` type($dofmap) `->` type($output)"
    )

    def __init__(
        self, vector: SSAValue | Operation, dofmap: SSAValue | Operation, element: ElementAttr
    ):
        shape = SSAValue.get(dofmap).type.get_shape()
        super().__init__(
            operands=[vector, dofmap],
            result_types=[TensorType(f64, shape)],
            properties={"element": element},
        )

    def verify_(self) -> None:
        if len(self.vector.type.get_shape()) != 1:
            raise VerifyException(f"{self.name} expects a vector of rank 1")
        check_shape(self, self.dofmap, [num_dofs(self.element)], "dofmap")
        check_shape(self, self.output, [num_dofs(self.element)], "output")


@irdl_op_definition
class ScatterAddOp(IRDLOperation):
    """Sum values on each cell into an assembled vector.

    output[j] = sum of input over the flat E-vector indices transpose[j, :], where an
    index c * ndofs + i refers to input[c, i] and -1 pads the rows. Each dof sums a
    fixed number of entries in a fixed order, so the result is deterministic.
    """

    name = "uflx.scatter_add"

    input = operand_def(TensorType[Float64Type])
    transpose = operand_def(TensorType[IntegerType])
    output = result_def(TensorType[Float64Type])

    traits = traits_def(Pure())

    assembly_format = (
        "$input `,` $transpose attr-dict `:` type($input) `,` type($transpose) `->` type($output)"
    )

    def __init__(self, input: SSAValue | Operation, transpose: SSAValue | Operation):
        num_dofs = SSAValue.get(transpose).type.get_shape()[0]
        super().__init__(operands=[input, transpose], result_types=[TensorType(f64, [num_dofs])])

    def verify_(self) -> None:
        if len(self.input.type.get_shape()) != 2 or len(self.transpose.type.get_shape()) != 2:
            raise VerifyException(f"{self.name} expects an input and transpose of rank 2")


UFLx = Dialect("uflx", [CoefficientOp, JacobianOp, IntegralOp, GatherOp, ScatterAddOp])
