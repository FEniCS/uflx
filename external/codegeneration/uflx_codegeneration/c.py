"""Generation of C code."""

import re
from collections.abc import Callable
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
# Stand-in for the code of an operand that has not been generated yet.
_MISSING = "\x01"


class CGenerator(Generator):
    """Generator of C code for lowered UFLx expressions and code-structure nodes.

    Expressions are generated as nested C expressions. With common-subexpression elimination,
    each statement first binds the subexpressions it uses more than once to temporaries.

    Expressions are generated without recursion (see _evaluate), so the depth of an expression
    is not limited by Python's recursion limit.
    """

    target = "C"

    def __init__(self, cse: bool):
        """Initialise.

        Args:
            cse: Evaluate each common subexpression of a statement once, into a temporary.
        """
        self.cse = cse
        # While _evaluate runs a handler: gives the code of each operand the handler asks for.
        self._operand: Callable[[Any], str] | None = None

    def emit(self, expression: Any) -> str:
        """Generate code for an expression that is the operand of another expression."""
        if self._operand is not None:
            return self._operand(expression)
        return self._expand(expression)

    def _evaluate(
        self, root: Any, operand_code: Callable[[Any, dict[Any, tuple[str, list[Any]]]], str]
    ) -> dict[Any, tuple[str, list[Any]]]:
        """Run the handler of root and of every expression its code needs, without recursion.

        A handler gets the code of each operand from operand_code. If an operand's handler has
        not run yet, the handler is run again once it has: handlers only build strings, so
        running one twice is harmless.

        Args:
            root: The expression.
            operand_code: The code to use for an operand, given the results so far.

        Returns:
            For each expression, its code and the operands its handler asked for, in order.
        """
        results: dict[Any, tuple[str, list[Any]]] = {}
        operands: list[Any] = []
        missing: list[Any] = []

        def operand(expression: Any) -> str:
            operands.append(expression)
            if expression in results:
                return operand_code(expression, results)
            missing.append(expression)
            return _MISSING

        outer, self._operand = self._operand, operand
        try:
            stack = [root]
            while stack:
                node = stack[-1]
                if node in results:
                    stack.pop()
                    continue
                operands.clear()
                missing.clear()
                code = self._handler(type(node))(self, node)
                if missing:
                    stack.extend(missing)
                else:
                    results[node] = (code, list(operands))
                    stack.pop()
        finally:
            self._operand = outer
        return results

    def _expand(self, expression: Any) -> str:
        """Generate code for an expression, with every subexpression written out in place."""
        return self._evaluate(expression, lambda e, results: results[e][0])[expression][0]

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
            return [], self._expand(expression)

        # Counting pass: each expression's code with placeholders for its operands, so no code
        # is duplicated even for heavily shared expressions.
        ids: dict[Any, int] = {}

        def placeholder(e: Any, results: Any) -> str:
            return f"\x00{ids.setdefault(e, len(ids))}\x00"

        counting = self._evaluate(expression, placeholder)
        # An expression whose code is exactly its operand's code (eg a PointComponent) is an
        # alias: reusing it also reuses the expression it forwards to.
        alias = {
            node: operands[0]
            for node, (code, operands) in counting.items()
            if len(operands) == 1 and code == placeholder(operands[0], counting)
        }

        # Count the uses of each expression, visiting operands in the order the handlers ask
        # for them and expanding each expression on its first use.
        uses = {expression: 1}
        iterators = [iter(counting[expression][1])]
        while iterators:
            operand = next(iterators[-1], None)
            if operand is None:
                iterators.pop()
            elif operand in uses:
                node: Any = operand
                while node is not None:
                    uses[node] += 1
                    node = alias.get(node)
            else:
                uses[operand] = 1
                iterators.append(iter(counting[operand][1]))
        temporaries = {
            node
            for node, n in uses.items()
            if n > 1 and node not in alias and not _ATOMIC.fullmatch(counting[node][0])
        }

        # Name the temporaries in the order of their declarations: on first use, after the
        # temporaries they depend on.
        names: dict[Any, str] = {}
        visited = {expression}
        stack = [(expression, iter(counting[expression][1]))]
        while stack:
            node, operands = stack[-1]
            operand = next(operands, None)
            if operand is None:
                stack.pop()
                if node in temporaries:
                    names[node] = symbols.global_variable_namer.temporary()
            elif operand not in visited:
                visited.add(operand)
                stack.append((operand, iter(counting[operand][1])))

        # Emission pass: refer to temporaries by name.
        emission = self._evaluate(
            expression, lambda e, results: names[e] if e in names else results[e][0]
        )
        declarations = [
            f"const double {name} = {emission[node][0]};" for node, name in names.items()
        ]
        return declarations, emission[expression][0]

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
