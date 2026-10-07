"""Generation of MLIR.

The generated module uses the func, scf, arith, math and memref dialects. Tables are module-level
``memref.global`` constants. The kernel takes its arrays as dynamically sized memrefs and is
marked ``llvm.emit_c_interface``, so after lowering to the LLVM dialect it can be called through
``_mlir_ciface_tabulate_tensor_f64`` with pointers to memref descriptors.
"""

import math
from typing import Any

import numpy as np
from uflx.expressions import (
    Abs,
    AbstractExpression,
    Div,
    Integer,
    Neg,
    Product,
    RealScalar,
    Reciprocal,
    Subtract,
    Sum,
)
from uflx.geometry import CoordinateDofComponent
from uflx.points import PointComponent

from uflx_codegeneration import symbols
from uflx_codegeneration.generator import Generator, handles
from uflx_codegeneration.nodes import AddToLocalTensor, ArrayEntry, FunctionCall, Loop, Variable
from uflx_codegeneration.quadrature import QuadratureLoop

# MLIR types of the C types used for kernel and geometry function arguments
_TYPES = {"double*": "memref<?xf64>", "const double*": "memref<?xf64>", "int": "index"}


def _float(value: float) -> str:
    """Format a float as an MLIR floating point literal."""
    if not math.isfinite(value):
        raise ValueError(f"Cannot write the non-finite value {value} as an MLIR literal")
    return format(value, ".17e")


def _dense(table: np.ndarray) -> str:
    """Format a table as the body of an MLIR dense elements attribute."""
    if table.ndim == 1:
        return "[" + ", ".join(_float(float(v)) for v in table) + "]"
    return "[" + ", ".join(_dense(t) for t in table) + "]"


def _memref_type(shape: tuple[int, ...]) -> str:
    """The MLIR type of a statically shaped f64 memref."""
    return "memref<" + "".join(f"{n}x" for n in shape) + "f64>"


class MLIRGenerator(Generator):
    """Generator of MLIR for lowered UFLx expressions and code-structure nodes.

    MLIR is in SSA form: every operation defines a value, and operations cannot be nested. The
    expression handlers therefore append operations to the current block and return the name of
    the SSA value they define. Values are memoised per region (function body or loop body), so
    an expression that is used more than once is computed once, without a separate
    common-subexpression elimination pass.
    """

    target = "MLIR"

    def __init__(self, tables: dict[str, np.ndarray]):
        """Initialise.

        Args:
            tables: The tables used by the kernel and its geometry functions.
        """
        self.tables = tables
        self._lines: list[str] = []
        self._indent = 0
        self._regions: list[dict[Any, str]] = []
        self._next = 0

    # Building blocks

    def _append(self, line: str):
        """Append an operation to the current block."""
        self._lines.append("  " * self._indent + line)

    def _define(self, operation: str) -> str:
        """Append an operation that defines a value, and return the value's name."""
        name = f"%v{self._next}"
        self._next += 1
        self._append(f"{name} = {operation}")
        return name

    def _lookup(self, key: Any) -> str | None:
        """Find a value computed in the current region or a region enclosing it."""
        for region in reversed(self._regions):
            if key in region:
                return region[key]
        return None

    def value(self, expression: Any) -> str:
        """Get the SSA value of an expression, computing it if it is not yet available."""
        name = self._lookup(expression)
        if name is None:
            name = self._handler(type(expression))(self, expression)
            self._regions[-1][expression] = name
        return name

    def _index(self, i: int | str) -> str:
        """Get an index value: a constant, or the induction variable of an enclosing loop."""
        if isinstance(i, str) and not i.isdigit():
            name = self._lookup(("index variable", i))
            if name is None:
                raise NotImplementedError(f"Cannot generate MLIR for the index expression {i!r}")
            return name
        key = ("index constant", int(i))
        name = self._lookup(key)
        if name is None:
            name = self._define(f"arith.constant {int(i)} : index")
            self._regions[-1][key] = name
        return name

    def _type(self, value: Any) -> str:
        """The MLIR type of an argument of a function call."""
        if isinstance(value, Variable):
            return _TYPES[value._dtype]
        if isinstance(value, AbstractExpression):
            return "f64"
        return "index"

    def _region(self, header: str, body: Any, bindings: dict[Any, str]):
        """Generate a region (eg a loop body) with its own scope of values."""
        self._append(header + " {")
        self._indent += 1
        self._regions.append(dict(bindings))
        try:
            self.code(body)
        finally:
            self._regions.pop()
            self._indent -= 1
        self._append("}")

    def function(
        self,
        name: str,
        arguments: list[tuple[str, str]],
        tables: dict[str, np.ndarray],
        body: Any,
        returns: bool,
        attributes: str = "",
    ) -> str:
        """Generate a function.

        Args:
            name: The name of the function.
            arguments: The name and MLIR type of each argument.
            tables: The tables the function reads; each is fetched once on entry.
            body: The function body: a statement, or (if returns is True) the expression whose
                value is returned.
            returns: Whether the function returns the f64 value of body.
            attributes: Attributes of the function, eg "attributes {llvm.emit_c_interface}".

        Returns:
            The function.
        """
        self._lines, self._indent, self._next = [], 1, 0
        bindings: dict[Any, str] = {}
        for argument, mlir_type in arguments:
            key = ("index variable", argument) if mlir_type == "index" else ("argument", argument)
            bindings[key] = f"%{argument}"
        self._regions = [bindings]
        for table_name, table in tables.items():
            bindings[("table", table_name)] = self._define(
                f"memref.get_global @{table_name} : {_memref_type(table.shape)}"
            )
        if returns:
            self._append(f"return {self.value(body)} : f64")
        else:
            self.code(body)
            self._append("return")
        signature = ", ".join(f"%{a}: {t}" for a, t in arguments)
        result = " -> f64" if returns else ""
        visibility = "private " if returns else ""
        return (
            f"func.func {visibility}@{name}({signature}){result} {attributes}".rstrip()
            + " {\n"
            + "\n".join(self._lines)
            + "\n}"
        )

    def module(
        self,
        geometry_functions: dict[str, tuple[str, list[Variable], dict[str, np.ndarray], Any]],
        kernel_tables: dict[str, np.ndarray],
        kernel: Any,
    ) -> str:
        """Generate the module for a kernel: tables, geometry functions and the kernel.

        Args:
            geometry_functions: For each geometry function, its return type, arguments,
                tables and the expression it returns.
            kernel_tables: The tables used by the kernel.
            kernel: The kernel body.

        Returns:
            The MLIR module.
        """
        parts = [
            f'memref.global "private" constant @{name} : {_memref_type(table.shape)} = '
            f"dense<{_dense(table)}>"
            for name, table in self.tables.items()
        ]
        for name, (dtype, inputs, tables, expression) in geometry_functions.items():
            if dtype != "double":
                raise NotImplementedError(f"Cannot generate MLIR for a function returning {dtype}")
            arguments = [(i._variable, _TYPES[i._dtype]) for i in inputs]
            parts.append(self.function(name, arguments, tables, expression, True))
        arguments = [
            (symbols.local_tensor, "memref<?xf64>"),
            (symbols.coefficients, "memref<?xf64>"),
            (symbols.constants, "memref<?xf64>"),
            (symbols.coordinate_dofs, "memref<?xf64>"),
            (symbols.entity_local_index, "memref<?xi32>"),
            (symbols.quadrature_permutation, "memref<?xi8>"),
            (symbols.custom_data, "!llvm.ptr"),
        ]
        parts.append(
            self.function(
                "tabulate_tensor_f64",
                arguments,
                kernel_tables,
                kernel,
                False,
                "attributes {llvm.emit_c_interface}",
            )
        )
        return "module {\n" + "\n\n".join(parts) + "\n}\n"

    # Expressions: each handler returns the name of the SSA value it defines

    @handles(RealScalar, Integer)
    def _scalar(self, node: RealScalar | Integer) -> str:
        return self._define(f"arith.constant {_float(float(node.value))} : f64")

    @handles(Sum)
    def _sum(self, node: Sum) -> str:
        values = [self.value(i) for i in node._items]
        result = values[0]
        for v in values[1:]:
            result = self._define(f"arith.addf {result}, {v} : f64")
        return result

    @handles(Product)
    def _product(self, node: Product) -> str:
        if node.value_shape != ():
            raise NotImplementedError("Cannot generate code for multiplication of non-scalars")
        values = [self.value(i) for i in node._items]
        result = values[0]
        for v in values[1:]:
            result = self._define(f"arith.mulf {result}, {v} : f64")
        return result

    @handles(Div)
    def _div(self, node: Div) -> str:
        first, second = self.value(node.first), self.value(node.second)
        return self._define(f"arith.divf {first}, {second} : f64")

    @handles(Subtract)
    def _subtract(self, node: Subtract) -> str:
        first, second = self.value(node.first), self.value(node.second)
        return self._define(f"arith.subf {first}, {second} : f64")

    @handles(Reciprocal)
    def _reciprocal(self, node: Reciprocal) -> str:
        argument = self.value(node.argument)
        one = self.value(RealScalar(1.0))
        return self._define(f"arith.divf {one}, {argument} : f64")

    @handles(Abs)
    def _abs(self, node: Abs) -> str:
        return self._define(f"math.absf {self.value(node.argument)} : f64")

    @handles(Neg)
    def _neg(self, node: Neg) -> str:
        return self._define(f"arith.negf {self.value(node.argument)} : f64")

    @handles(PointComponent)
    def _point_component(self, node: PointComponent) -> str:
        if not isinstance(node.component_index, int):
            raise NotImplementedError("Cannot generate code for a symbolic point component")
        component = node.point.component(node.component_index)
        if isinstance(component, PointComponent):
            raise NotImplementedError("Cannot generate code for a symbolic point component")
        return self.value(component)

    @handles(CoordinateDofComponent)
    def _coordinate_dof_component(self, node: CoordinateDofComponent) -> str:
        coordinate_dofs = self._lookup(("argument", symbols.coordinate_dofs))
        index = self._index(node._tdim * node._point + node._component)
        return self._define(f"memref.load {coordinate_dofs}[{index}] : memref<?xf64>")

    @handles(ArrayEntry)
    def _array_entry(self, node: ArrayEntry) -> str:
        table = self._lookup(("table", node.array))
        if table is None:
            raise NotImplementedError(f"Cannot generate MLIR for an entry of {node.array}")
        indices = ", ".join(self._index(i) for i in node.index)
        shape = self.tables[node.array].shape
        return self._define(f"memref.load {table}[{indices}] : {_memref_type(shape)}")

    @handles(FunctionCall)
    def _function_call(self, node: FunctionCall) -> str:
        arguments, types = [], []
        for i in node.inputs:
            if isinstance(i, AbstractExpression):
                arguments.append(self.value(i))
            else:
                arguments.append(self._index(i))
            types.append(self._type(i))
        return self._define(
            f"func.call @{node.function}({', '.join(arguments)}) : ({', '.join(types)}) -> f64"
        )

    @handles(Variable)
    def _variable(self, node: Variable) -> str:
        key = ("index variable" if _TYPES[node._dtype] == "index" else "argument", node._variable)
        name = self._lookup(key)
        if name is None:
            raise NotImplementedError(f"Variable {node._variable} is not defined")
        return name

    # Statements: each handler appends operations and returns ""

    @handles(AddToLocalTensor)
    def _add_to_local_tensor(self, node: AddToLocalTensor) -> str:
        value = self.value(node.body)
        index = self._index(node.component[0])
        for component, size in zip(node.component[1:], node.shape[1:]):
            index = self._define(f"arith.muli {index}, {self._index(size)} : index")
            index = self._define(f"arith.addi {index}, {self._index(component)} : index")
        tensor = self._lookup(("argument", symbols.local_tensor))
        old = self._define(f"memref.load {tensor}[{index}] : memref<?xf64>")
        new = self._define(f"arith.addf {old}, {value} : f64")
        self._append(f"memref.store {new}, {tensor}[{index}] : memref<?xf64>")
        return ""

    def _loop(self, variable: str, start: int | str, end: int | str, body: Any):
        """Generate an scf.for loop."""
        lower, upper, step = self._index(start), self._index(end), self._index(1)
        induction = f"%{variable}"
        self._region(
            f"scf.for {induction} = {lower} to {upper} step {step}",
            body,
            {("index variable", variable): induction},
        )

    @handles(Loop)
    def _for(self, node: Loop) -> str:
        self._loop(node.variable, node.start, node.end, node.body)
        return ""

    @handles(QuadratureLoop)
    def _quadrature_loop(self, node: QuadratureLoop) -> str:
        self._loop(node.variable, 0, node.rule.npoints, node.body)
        return ""
