"""Test code generation for integrals over equal measures."""

from typing import Any

import numpy as np
from cffi import FFI
from uflx import TestFunction, TrialFunction, coordinate_element, dx, function_space, grad, inner
from uflx.integrals import Measure

import uflx_codegeneration


def tabulate(form, name, code_dir, coords, shape):
    """Generate, compile and run the kernel of a form on one cell."""
    code, signature = uflx_codegeneration.generate(form)
    ffi = FFI()
    ffi.cdef(signature)
    ffi.set_source(name, code)
    lib: Any = ffi.dlopen(ffi.compile(code_dir))
    result = np.zeros(shape)
    empty = np.zeros(0)
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


def test_equal_measure(lagrange_element, code_dir):
    """A form over Measure(codim=0) has the same element matrix as over dx."""
    element = lagrange_element("triangle", 1)
    space = function_space(coordinate_element(lagrange_element("triangle", 1, (2,))), element)
    u = TrialFunction(space)
    v = TestFunction(space)
    integrand = inner(grad(u), grad(v)) + u * v
    coords = np.array([[0.0, 0.0], [1.0, 0.2], [0.3, 1.0]])

    expected = tabulate(integrand * dx, "test_measure_dx", code_dir, coords, (3, 3))
    result = tabulate(integrand * Measure(codim=0), "test_measure_codim0", code_dir, coords, (3, 3))
    assert np.allclose(result, expected, rtol=1e-14, atol=0)
