from itertools import product
from typing import Any

import basix
from uflx import dx
from uflx.algorithms import pull_back_to_reference, replace, simplify
from uflx.expressions import (
    Abs,
    AbstractExpression,
    Div,
    Integer,
    Neg,
    Product,
    Rational,
    RealScalar,
    Reciprocal,
    ScalarMult,
    Subtract,
    Sum,
)
from uflx.functions import Coefficient, TestFunction
from uflx.geometry import (
    Jacobian,
    JacobianDeterminant,
    JacobianInverse,
    JacobianInverseTranspose,
    JacobianTranspose,
)
from uflx.graphs import GraphNode, as_graph
from uflx.integrals import AbstractIntegral
from uflx.maps import apply_push_forwards
from uflx.operators import Inner, ReferenceGrad
from uflx.tensors import Matrix, Vector
from xdsl.builder import ImplicitBuilder
from xdsl.dialects import func
from xdsl.dialects.builtin import DYNAMIC_INDEX, ModuleOp, TensorType, f64
from xdsl.ir import SSAValue

from ffcy import fields
from ffcy.dialects import uflx
from ffcy.dialects.basix import CellType, ElementAttr, QuadratureAttr
from ffcy.fields import Field


def basix_element(element) -> basix.finite_element.FiniteElement:
    # basix_uflx has no public accessor for the wrapped Basix element.
    e = getattr(element, "sub_element", element)._element
    if not e.has_tensor_product_factorisation:
        raise ValueError(f"{e} must use Basix tensor product DOF ordering")
    return e


def e_vector_type(*shape: int) -> TensorType:
    return TensorType(f64, [DYNAMIC_INDEX, *shape])


class JacobianComponent(AbstractExpression):
    """The derivative of x_g along reference axis r, at the quadrature points."""

    def __init__(self, g: int, r: int):
        self.g = g
        self.r = r

    @property
    def value_shape(self) -> tuple[int, ...]:
        return ()

    @property
    def is_real_valued(self) -> bool:
        return True

    @property
    def successors(self) -> set[GraphNode]:
        return set()

    @property
    def init_args(self) -> tuple[Any, ...]:
        return self.g, self.r

    def component(self, *indices: int) -> AbstractExpression:
        raise ValueError("Cannot get a component of a scalar expression")


def expand_geometry(nodes, gdim: int, tdim: int) -> dict:
    """Write the Jacobian and quantities derived from it in terms of its components."""
    J = Matrix([[JacobianComponent(g, r) for r in range(tdim)] for g in range(gdim)])
    expansions = {
        Jacobian: lambda: J,
        JacobianTranspose: J.transpose,
        JacobianDeterminant: J.compute_determinant,
        JacobianInverse: J.compute_inverse,
        JacobianInverseTranspose: lambda: J.compute_inverse().transpose(),
    }
    return {n: expansions[type(n)]() for n in nodes if type(n) in expansions}


def indices(shape: tuple[int, ...]):
    return product(*(range(n) for n in shape))


class Translator:
    """Translates scalar expressions into fields at the quadrature points."""

    def __init__(self, coordinate_dofs: SSAValue, element: ElementAttr, quadrature: QuadratureAttr):
        self.coordinate_dofs = coordinate_dofs
        self.element = element
        self.quadrature = quadrature
        self.translated: dict[AbstractExpression, Field] = {}

    def __call__(self, expr: AbstractExpression) -> Field:
        if expr not in self.translated:
            self.translated[expr] = self.translate(expr)
        return self.translated[expr]

    def translate(self, expr: AbstractExpression) -> Field:
        match expr:
            case Integer() | RealScalar():
                return float(expr.value)
            case Rational():
                return expr.numerator / expr.denominator
            case JacobianComponent():
                component = [expr.g, expr.r]
                return uflx.JacobianOp(
                    self.coordinate_dofs, self.element, self.quadrature, component
                ).output
            case Abs():
                return fields.absolute(self(expr.argument))
            case Neg():
                return fields.neg(self(expr.argument))
            case Sum():
                items = expr.init_args[0]
                result = fields.total([self(i) for i in items if not isinstance(i, Neg)])
                for i in items:
                    if isinstance(i, Neg):
                        result = fields.sub(result, self(i.argument))
                return result
            case Product():
                items = expr.init_args[0]
                numerator = fields.product(
                    [self(i) for i in items if not isinstance(i, Reciprocal)]
                )
                denominator = fields.product(
                    [self(i.argument) for i in items if isinstance(i, Reciprocal)]
                )
                return fields.div(numerator, denominator)
            case Subtract():
                return fields.sub(self(expr.first), self(expr.second))
            case ScalarMult():
                return fields.mul(self(expr.first), self(expr.second))
            case Div():
                return fields.div(self(expr.first), self(expr.second))
            case Inner():
                # Expressions are real-valued, so the inner product needs no conjugate.
                a, b = expr.first, expr.second
                return fields.total(
                    [
                        fields.mul(self(a.component(*i)), self(b.component(*i)))
                        for i in indices(a.value_shape)
                    ]
                )
            case _:
                raise NotImplementedError(f"Cannot translate {type(expr).__name__}")


def derivatives(d: int) -> list[tuple[int, ...]]:
    """The value, then the first reference derivative along each axis."""
    return [(0,) * d] + [tuple(int(a == r) for a in range(d)) for r in range(d)]


def from_uflx(integral: AbstractIntegral, quadrature_degree: int) -> ModuleOp:
    """Translate the action of a bilinear form, linear in one coefficient u.

    The integrand is written as sum_kj (D_k v) C_kj (D_j u), where D are values and
    first reference derivatives. Each C_kj is found by substituting values and
    unit vectors for the evaluations of u and v.
    """
    if integral.measure != dx:
        raise NotImplementedError("Only cell integrals are supported")
    integrand = simplify(apply_push_forwards(pull_back_to_reference(integral))).integrand
    nodes = as_graph(integrand).nodes
    tests = {n for n in nodes if isinstance(n, TestFunction)}
    coefficients = {n for n in nodes if isinstance(n, Coefficient)}
    if len(tests) != 1 or len(coefficients) != 1:
        raise NotImplementedError("The integrand must have one test function and one coefficient")
    ((test,), (u,)) = (tests, coefficients)
    gradients = {n.argument: n for n in nodes if isinstance(n, ReferenceGrad)}

    test_element = basix_element(test.function_space.elements[0])
    u_element = basix_element(u.function_space.elements[0])
    domain = test.function_space.domain
    coordinate_element = basix_element(domain.elements[0])
    quadrature = QuadratureAttr(CellType(coordinate_element.cell_type.name), quadrature_degree)
    d = len(test_element.get_tensor_product_representation()[0])
    modes = derivatives(d)
    integrand = simplify(replace(integrand, expand_geometry(nodes, domain.geometric_dimension, d)))

    def substitution(f, derivative: tuple[int, ...]) -> dict:
        values = {f: Integer(int(not any(derivative)))}
        if f in gradients:
            values[gradients[f]] = Vector([Integer(a) for a in derivative])
        return values

    C = {
        (k, j): c
        for k in modes
        for j in modes
        if not (
            c := simplify(replace(integrand, substitution(test, k) | substitution(u, j)))
        ).is_zero
    }

    inputs = [
        e_vector_type(coordinate_element.dim, domain.geometric_dimension),
        e_vector_type(u_element.dim),
    ]
    function = func.FuncOp("action", (inputs, [e_vector_type(test_element.dim)]))
    coordinate_dofs, u_dofs = function.args
    coordinate_dofs.name_hint = "coordinate_dofs"
    u_dofs.name_hint = u.label.replace("-", "_")

    with ImplicitBuilder(function.body):
        translate = Translator(
            coordinate_dofs, ElementAttr.from_basix(coordinate_element), quadrature
        )
        u_attr = ElementAttr.from_basix(u_element)
        evaluations = {
            j: uflx.CoefficientOp(u_dofs, u_attr, quadrature, j).output
            for j in modes
            if any(jj == j for _, jj in C)
        }
        fluxes = {
            k: fields.total(
                [fields.mul(translate(c), evaluations[j]) for (kk, j), c in C.items() if kk == k]
            )
            for k in modes
            if any(kk == k for kk, _ in C)
        }
        terms = [(k, f) for k, f in fluxes.items() if not fields.is_constant(f, 0)]
        integral_op = uflx.IntegralOp(
            [fields.field(f) for _, f in terms],
            ElementAttr.from_basix(test_element),
            quadrature,
            [k for k, _ in terms],
        )
        func.ReturnOp(integral_op)

    return ModuleOp([function])
