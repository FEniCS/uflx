"""Test common-subexpression elimination in the generated C."""

import re
from typing import Any

import numpy as np
import pytest
from cffi import FFI
from uflx import (
    SpatialCoordinate,
    TestFunction,
    TrialFunction,
    coordinate_element,
    dx,
    function_space,
    grad,
    inner,
)
from uflx.expressions import AbstractExpression, Sum
from uflx.geometry import JacobianDeterminant
from uflx.points import Point, PointComponent

import uflx_codegeneration
from uflx_codegeneration import symbols
from uflx_codegeneration.c import CGenerator
from uflx_codegeneration.nodes import ArrayEntry


@pytest.fixture(autouse=True)
def reset_names():
    """Number temporaries from cse0 in every test."""
    symbols.global_variable_namer.reset()


def _evaluate(declarations: list[str], code: str, w: list[float]) -> float:
    """Evaluate generated scalar C code (temporaries plus an expression) in Python."""
    names: dict[str, Any] = {"w": w, "fabs": abs}
    for declaration in declarations:
        match = re.fullmatch(r"const double (\w+) = (.*);", declaration)
        assert match is not None
        names[match[1]] = eval(match[2], names)
    return eval(code, names)


def _build(leaf: AbstractExpression, n: int) -> AbstractExpression:
    """Build e_n = e_{n-1} / e_{n-1} - e_{n-1}, which uses e_{n-1} three times."""
    e = leaf
    for _ in range(n):
        e = e / e - e
    return e


def test_repeated_subexpression_is_bound_once():
    """A non-atomic subexpression used twice is evaluated once."""
    s = ArrayEntry("w", (0,)) + ArrayEntry("w", (1,))
    assert CGenerator(cse=True).statement(s * s) == (
        ["const double cse0 = (w[0] + w[1]);"],
        "(cse0 * cse0)",
    )


def test_atomic_subexpressions_are_not_bound():
    """Array entries, names and literals are not bound to temporaries."""
    a = ArrayEntry("w", (0,))
    assert CGenerator(cse=True).statement(a * a) == ([], "(w[0] * w[0])")


def test_uses_through_an_alias_are_counted():
    """A PointComponent forwards its component's code, so its uses count for the component."""
    s = ArrayEntry("w", (0,)) + ArrayEntry("w", (1,))
    component = PointComponent(Point([s, ArrayEntry("w", (2,))]), 0)
    declarations, code = CGenerator(cse=True).statement(component * s)
    assert declarations == ["const double cse0 = (w[0] + w[1]);"]
    assert code == "(cse0 * cse0)"


def test_cse_matches_the_expanded_expression():
    """The code with temporaries evaluates to the same value as the fully expanded code."""
    e = _build(ArrayEntry("w", (0,)), 6)
    w = [0.3]
    declarations, code = CGenerator(cse=True).statement(e)
    assert len(declarations) > 0
    assert _evaluate(declarations, code, w) == _evaluate([], CGenerator(cse=False).code(e), w)


def test_cse_code_is_linear_in_the_number_of_nodes():
    """Fully expanded, this expression has ~3^200 terms; with CSE its code is linear in size."""
    n = 200
    declarations, code = CGenerator(cse=True).statement(_build(ArrayEntry("w", (0,)), n))
    assert len(declarations) <= 4 * n
    assert len(code) + sum(len(d) for d in declarations) < 100 * n


def test_generator_does_not_patch_expression_classes():
    """C generation lives on the generator, not as methods added to UFLx's classes."""
    assert not hasattr(Sum, "generate_c")
    assert not hasattr(ArrayEntry, "generate_c")


def test_unsupported_node_raises():
    """A node without a C handler raises a NotImplementedError naming its type."""
    with pytest.raises(NotImplementedError, match="SingleSpatialCoordinate"):
        CGenerator(cse=True).code(SpatialCoordinate(2)[0])


def test_generators_are_independent():
    """Generators with different options can be used side by side."""
    s = ArrayEntry("w", (0,)) + ArrayEntry("w", (1,))
    with_cse, without_cse = CGenerator(cse=True), CGenerator(cse=False)
    assert without_cse.statement(s * s) == ([], "((w[0] + w[1]) * (w[0] + w[1]))")
    assert with_cse.statement(s * s) == (["const double cse0 = (w[0] + w[1]);"], "(cse0 * cse0)")
    assert without_cse.statement(s * s) == ([], "((w[0] + w[1]) * (w[0] + w[1]))")


def _tabulate(code: str, signature: str, name: str, code_dir: str, shape: tuple[int, ...]):
    """Compile a kernel and tabulate it on one tilted cell."""
    ffi = FFI()
    ffi.cdef(signature)
    ffi.set_source(name, code)
    lib: Any = ffi.dlopen(ffi.compile(code_dir))
    coords = np.array([[0.1, 0.0], [1.3, 0.4], [0.2, 0.9]])
    empty = np.zeros(0)
    result = np.zeros(shape)
    lib.tabulate_tensor_f64(
        ffi.cast("double*", result.ctypes.data),
        ffi.cast("double*", empty.ctypes.data),
        ffi.cast("double*", empty.ctypes.data),
        ffi.cast("double*", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    return result


@pytest.mark.parametrize("degree", [1, 2, 3])
@pytest.mark.parametrize("kernel", ["mass", "stiffness", "load"])
def test_kernels_with_and_without_cse_agree(lagrange_element, code_dir, kernel, degree):
    """CSE changes the code, not the element tensor."""
    element = lagrange_element("triangle", degree)
    space = function_space(coordinate_element(lagrange_element("triangle", 1, (2,))), element)
    u, v = TrialFunction(space), TestFunction(space)
    x = SpatialCoordinate(2)
    form = {
        "mass": inner(u, v) * dx,
        "stiffness": inner(grad(u), grad(v)) * dx,
        "load": x[0] * v * dx,
    }[kernel]
    shape = (element.dim,) if kernel == "load" else (element.dim, element.dim)

    code, signature = uflx_codegeneration.generate(form)
    expanded_code, _ = uflx_codegeneration.generate(form, cse=False)
    assert "cse" not in expanded_code
    assert len(code) <= len(expanded_code)

    result = _tabulate(code, signature, f"test_cse_{kernel}_{degree}", code_dir, shape)
    expected = _tabulate(
        expanded_code, signature, f"test_no_cse_{kernel}_{degree}", code_dir, shape
    )
    np.testing.assert_allclose(result, expected, rtol=1e-13, atol=1e-15)


def _p1_space(lagrange_element):
    """A P1 space on a triangle with P1 geometry."""
    return function_space(
        coordinate_element(lagrange_element("triangle", 1, (2,))), lagrange_element("triangle", 1)
    )


def _nested_load(lagrange_element, n: int):
    """Load vector of e_n * v, with e_k = e_{k-1} * e_{k-1} - e_{k-1} and e_0 = x[0].

    Each level uses the previous one three times, so the integrand written out without
    common subexpressions has about 3^n copies of x[0].
    """
    e = SpatialCoordinate(2)[0]
    for _ in range(n):
        e = e * e - e
    return e * TestFunction(_p1_space(lagrange_element)) * dx


def test_nested_kernel_code_grows_linearly_with_cse(lagrange_element):
    """Without CSE the nested kernel's code grows ~9x per two levels; with CSE it barely grows."""
    sizes = {
        (n, cse): len(uflx_codegeneration.generate(_nested_load(lagrange_element, n), cse=cse)[0])
        for n in (0, 4, 6)
        for cse in (False, True)
    }
    for cse in (False, True):
        growth_4 = sizes[(4, cse)] - sizes[(0, cse)]
        growth_6 = sizes[(6, cse)] - sizes[(0, cse)]
        if cse:
            assert growth_6 - growth_4 < 500
        else:
            assert growth_6 > 4 * growth_4


@pytest.mark.parametrize("n", [4, 6])
def test_nested_kernel_with_and_without_cse_agree(lagrange_element, code_dir, n):
    """The nested kernel computes the same load vector with and without CSE."""
    form = _nested_load(lagrange_element, n)
    code, signature = uflx_codegeneration.generate(form)
    expanded_code, _ = uflx_codegeneration.generate(form, cse=False)

    result = _tabulate(code, signature, f"test_cse_nested{n}", code_dir, (3,))
    expected = _tabulate(expanded_code, signature, f"test_no_cse_nested{n}", code_dir, (3,))
    # The C compiler may contract a * b + c into a fused multiply-add differently once the
    # expression is split into statements, so allow for rounding.
    np.testing.assert_allclose(result, expected, rtol=1e-12)


def test_repeated_geometry_function_call_is_evaluated_once(lagrange_element, code_dir):
    """A geometry function used twice is called once.

    The C compiler cannot merge the two calls itself: in a shared library built with -fPIC
    the geometry functions are interposable, so the compiler may not assume they are pure.
    """
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    space = function_space(domain, lagrange_element("triangle", 1))
    u, v = TrialFunction(space), TestFunction(space)
    # |det J| written in the integrand, on top of the |det J| of the change of variables
    form = abs(JacobianDeterminant(domain)) * u * v * dx

    code, signature = uflx_codegeneration.generate(form)
    expanded_code, _ = uflx_codegeneration.generate(form, cse=False)

    def calls_in_kernel(c: str) -> int:
        return len(re.findall(r"\bgeo\d+\(", c.split("void tabulate_tensor_f64(")[1]))

    assert calls_in_kernel(code) == 1
    assert calls_in_kernel(expanded_code) == 2

    result = _tabulate(code, signature, "test_cse_det", code_dir, (3, 3))
    expected = _tabulate(expanded_code, signature, "test_no_cse_det", code_dir, (3, 3))
    np.testing.assert_allclose(result, expected, rtol=1e-13)
