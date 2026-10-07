"""Generation of C code."""

import re
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
from uflx_codegeneration.nodes import (
    AddToLocalTensor,
    ArrayEntry,
    FunctionCall,
    Loop,
    Variable,
    flatten_component,
)
from uflx_codegeneration.quadrature import QuadratureLoop
from uflx_codegeneration.utils import indented

# Code that is a literal, a name or an array entry with simple indices: never worth binding to
# a temporary.
_ATOMIC = re.compile(r"-?[\w.]+(\[[^\[\]\x00]*\])*")
# Placeholder for an operand in the counting pass of common-subexpression elimination.
_PLACEHOLDER = re.compile(r"\x00(\d+)\x00")


class _CSEScope:
    """State of common-subexpression elimination for one statement.

    In the counting pass, CGenerator.emit records how often each expression is used and returns
    a placeholder instead of the expression's code, so no code is duplicated even for heavily
    shared expressions. In the emission pass, the expressions chosen as temporaries are declared
    on first use and referred to by name afterwards.
    """

    def __init__(self, counting: bool, temporaries: set[Any]):
        """Initialise."""
        self.counting = counting
        self.temporaries = temporaries
        self.uses: dict[Any, int] = {}
        self.code: dict[Any, str] = {}
        self.ids: dict[Any, int] = {}
        self.by_id: list[Any] = []
        self.alias: dict[Any, Any] = {}
        self.names: dict[Any, str] = {}
        self.declarations: list[str] = []


class CGenerator(Generator):
    """Generator of C code for lowered UFLx expressions and code-structure nodes.

    Expressions are generated as nested C expressions. With common-subexpression elimination,
    each statement first binds the subexpressions it uses more than once to temporaries.
    """

    target = "C"

    def __init__(self, cse: bool):
        """Initialise.

        Args:
            cse: Evaluate each common subexpression of a statement once, into a temporary.
        """
        self.cse = cse
        self._scope: _CSEScope | None = None

    def emit(self, expression: Any) -> str:
        """Generate code for an expression that is the operand of another expression.

        Inside a statement with common-subexpression elimination, this counts the uses of the
        expression (counting pass), or returns the name of its temporary (emission pass).
        """
        scope = self._scope
        if scope is None:
            return self._handler(type(expression))(self, expression)

        if scope.counting:
            if expression in scope.uses:
                # A node whose code is exactly its operand's code (eg a PointComponent) is an
                # alias: reusing it also reuses the node it forwards to.
                node: Any = expression
                while node is not None:
                    scope.uses[node] += 1
                    node = scope.alias.get(node)
            else:
                scope.uses[expression] = 1
                code = self._handler(type(expression))(self, expression)
                scope.code[expression] = code
                if match := _PLACEHOLDER.fullmatch(code):
                    scope.alias[expression] = scope.by_id[int(match[1])]
                scope.ids[expression] = len(scope.by_id)
                scope.by_id.append(expression)
            return f"\x00{scope.ids[expression]}\x00"

        if expression in scope.names:
            return scope.names[expression]
        code = self._handler(type(expression))(self, expression)
        if expression in scope.temporaries:
            name = symbols.global_variable_namer.temporary()
            scope.declarations.append(f"const double {name} = {code};")
            scope.names[expression] = name
            return name
        return code

    def statement(self, expression: Any) -> tuple[list[str], str]:
        """Generate code for the scalar expression of one statement.

        With common-subexpression elimination, every subexpression that is used more than once
        and is not atomic (a literal, a name or an array entry) is evaluated once, into a
        ``const double`` temporary. Each temporary is declared on first use, after the
        temporaries it depends on.

        Args:
            expression: The expression.

        Returns:
            The declarations of the temporaries, and the code for the expression.
        """
        if not self.cse:
            return [], self.code(expression)

        outer = self._scope
        try:
            self._scope = counting = _CSEScope(True, set())
            self.emit(expression)
            self._scope = emission = _CSEScope(
                False,
                {
                    node
                    for node, uses in counting.uses.items()
                    if uses > 1
                    and node not in counting.alias
                    and not _ATOMIC.fullmatch(counting.code[node])
                },
            )
            code = self.emit(expression)
        finally:
            self._scope = outer
        return emission.declarations, code

    # Expressions

    @handles(Product)
    def _product(self, node: Product) -> str:
        if node.value_shape != ():
            raise NotImplementedError("Cannot generate code for multiplication of non-scalars")
        return "(" + " * ".join([self.emit(i) for i in node._items]) + ")"

    @handles(Sum)
    def _sum(self, node: Sum) -> str:
        return "(" + " + ".join([self.emit(i) for i in node._items]) + ")"

    @handles(Div)
    def _div(self, node: Div) -> str:
        return f"({self.emit(node.first)} / {self.emit(node.second)})"

    @handles(Subtract)
    def _subtract(self, node: Subtract) -> str:
        return f"({self.emit(node.first)} - {self.emit(node.second)})"

    @handles(Reciprocal)
    def _reciprocal(self, node: Reciprocal) -> str:
        return f"(1.0 / {self.emit(node.argument)})"

    @handles(Abs)
    def _abs(self, node: Abs) -> str:
        return f"fabs({self.emit(node.argument)})"

    @handles(Neg)
    def _neg(self, node: Neg) -> str:
        return f"-{self.emit(node.argument)}"

    @handles(PointComponent)
    def _point_component(self, node: PointComponent) -> str:
        if not isinstance(node.component_index, int):
            raise NotImplementedError("Cannot generate code for a symbolic point component")
        component = node.point.component(node.component_index)
        if isinstance(component, PointComponent):
            raise NotImplementedError("Cannot generate code for a symbolic point component")
        return self.emit(component)

    @handles(CoordinateDofComponent)
    def _coordinate_dof_component(self, node: CoordinateDofComponent) -> str:
        return f"{symbols.coordinate_dofs}[{node._tdim * node._point + node._component}]"

    @handles(RealScalar, Integer)
    def _scalar(self, node: RealScalar | Integer) -> str:
        return f"{node.value}"

    @handles(ArrayEntry)
    def _array_entry(self, node: ArrayEntry) -> str:
        return f"{node.array}[" + "][".join(f"{i}" for i in node.index) + "]"

    @handles(FunctionCall)
    def _function_call(self, node: FunctionCall) -> str:
        inputs = [
            self.emit(i) if isinstance(i, AbstractExpression) else f"{i}" for i in node.inputs
        ]
        return f"{node.function}(" + ", ".join(inputs) + ")"

    @handles(Variable)
    def _variable(self, node: Variable) -> str:
        return node._variable

    # Statements

    @handles(AddToLocalTensor)
    def _add_to_local_tensor(self, node: AddToLocalTensor) -> str:
        declarations, body = self.statement(node.body)
        entry = f"{symbols.local_tensor}[{flatten_component(node.component, node.shape)}]"
        return "\n".join([*declarations, f"{entry} += {body};"])

    @handles(Loop)
    def _loop(self, node: Loop) -> str:
        return (
            f"for (int {node.variable}={node.start}; {node.variable}!={node.end}; "
            f"++{node.variable})\n"
            "{\n" + indented(self.code(node.body), 2) + "\n}"
        )

    @handles(QuadratureLoop)
    def _quadrature_loop(self, node: QuadratureLoop) -> str:
        return (
            f"for (int {node.variable}=0; {node.variable}!={node.rule.npoints}; "
            f"++{node.variable})\n"
            "{\n" + indented(self.code(node.body), 2) + "\n}"
        )


def c_table(table: np.ndarray) -> str:
    """Convert a numpy array to C."""
    if len(table.shape) == 1:
        return "{" + ", ".join(f"{i}" for i in table) + "}"
    return "{" + ", ".join(c_table(i) for i in table) + "}"


def tables_to_c(tables: dict[str, np.ndarray]) -> str:
    """Convert tables of values to a string of code."""
    return "\n".join(
        f"static const double {variable}["
        + "][".join(f"{i}" for i in table.shape)
        + "] = "
        + c_table(table)
        + ";"
        for variable, table in tables.items()
    )
