"""Test generation of MLIR.

The generated MLIR is lowered to the LLVM dialect with mlir-opt, translated to LLVM IR with
mlir-translate and compiled with llc, and the element tensors it computes are compared with
those of the C kernels. The tests are skipped if these tools are not found: set UFLX_MLIR_BIN
to the directory containing them, or put them on the PATH.
"""

import ctypes
import glob
import os
import re
import shutil
import subprocess
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
from uflx.geometry import JacobianDeterminant

import uflx_codegeneration
from uflx_codegeneration import symbols
from uflx_codegeneration.mlir import MLIRGenerator
from uflx_codegeneration.nodes import ArrayEntry

_LOWERING = [
    "--convert-scf-to-cf",
    "--convert-math-to-llvm",
    "--convert-arith-to-llvm",
    "--expand-strided-metadata",
    "--finalize-memref-to-llvm",
    "--convert-func-to-llvm",
    "--convert-cf-to-llvm",
    "--reconcile-unrealized-casts",
]


def _find_tool(name: str) -> str | None:
    """Find an MLIR or LLVM tool."""
    directories = [os.environ.get("UFLX_MLIR_BIN", ""), *sorted(glob.glob("/usr/lib/llvm-*/bin"))]
    for directory in directories:
        if directory and os.path.isfile(os.path.join(directory, name)):
            return os.path.join(directory, name)
    return shutil.which(name)


_TOOLS = {name: _find_tool(name) for name in ("mlir-opt", "mlir-translate", "llc")}
requires_mlir = pytest.mark.skipif(
    any(path is None for path in _TOOLS.values()), reason="MLIR tools not found"
)


def _tool(name: str) -> str:
    """The path of a tool that was found."""
    path = _TOOLS[name]
    assert path is not None
    return path


@pytest.fixture(autouse=True)
def reset_names():
    """Number generated names from 0 in every test."""
    symbols.global_variable_namer.reset()


class _MemRef1D(ctypes.Structure):
    """The descriptor of a one-dimensional memref, as passed through the MLIR C interface."""

    _fields_ = [
        ("allocated", ctypes.c_void_p),
        ("aligned", ctypes.c_void_p),
        ("offset", ctypes.c_int64),
        ("sizes", ctypes.c_int64 * 1),
        ("strides", ctypes.c_int64 * 1),
    ]


def _memref(array: np.ndarray) -> _MemRef1D:
    """Describe a contiguous array as a one-dimensional memref."""
    pointer = array.ctypes.data
    return _MemRef1D(pointer, pointer, 0, (ctypes.c_int64 * 1)(array.size), (ctypes.c_int64 * 1)(1))


_COORDS = np.array([[0.1, 0.0], [1.3, 0.4], [0.2, 0.9]])


def _tabulate_mlir(module: str, directory: str, shape: tuple[int, ...]) -> np.ndarray:
    """Compile a generated MLIR module and tabulate its kernel on one tilted cell."""
    source = os.path.join(directory, "kernel.mlir")
    with open(source, "w") as f:
        f.write(module)
    llvm_dialect = os.path.join(directory, "kernel.llvm.mlir")
    llvm_ir = os.path.join(directory, "kernel.ll")
    obj = os.path.join(directory, "kernel.o")
    library = os.path.join(directory, "libkernel.so")
    subprocess.run([_tool("mlir-opt"), source, *_LOWERING, "-o", llvm_dialect], check=True)
    subprocess.run(
        [_tool("mlir-translate"), "--mlir-to-llvmir", llvm_dialect, "-o", llvm_ir], check=True
    )
    subprocess.run(
        [_tool("llc"), "-O2", "-relocation-model=pic", "-filetype=obj", llvm_ir, "-o", obj],
        check=True,
    )
    subprocess.run(["cc", "-shared", obj, "-o", library], check=True)

    lib = ctypes.CDLL(library)
    result = np.zeros(shape)
    coordinate_dofs = _COORDS.flatten()
    empty_f64, empty_i32, empty_i8 = np.zeros(0), np.zeros(0, np.int32), np.zeros(0, np.int8)
    descriptors = [
        _memref(a) for a in (result, empty_f64, empty_f64, coordinate_dofs, empty_i32, empty_i8)
    ]
    lib._mlir_ciface_tabulate_tensor_f64(*(ctypes.byref(d) for d in descriptors), None)
    return result


def _tabulate_c(code: str, signature: str, name: str, code_dir: str, shape: tuple[int, ...]):
    """Compile a generated C kernel and tabulate it on one tilted cell."""
    ffi = FFI()
    ffi.cdef(signature)
    ffi.set_source(name, code)
    lib: Any = ffi.dlopen(ffi.compile(code_dir))
    result = np.zeros(shape)
    coordinate_dofs = _COORDS.flatten()
    empty = np.zeros(0)
    lib.tabulate_tensor_f64(
        ffi.cast("double*", result.ctypes.data),
        ffi.cast("double*", empty.ctypes.data),
        ffi.cast("double*", empty.ctypes.data),
        ffi.cast("double*", coordinate_dofs.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    return result


def _form(lagrange_element, kernel: str, degree: int):
    """A form to compile, and the shape of its element tensor."""
    domain = coordinate_element(lagrange_element("triangle", 1, (2,)))
    element = lagrange_element("triangle", degree)
    space = function_space(domain, element)
    u, v = TrialFunction(space), TestFunction(space)
    n = element.dim
    if kernel == "mass":
        return inner(u, v) * dx, (n, n)
    if kernel == "stiffness":
        return inner(grad(u), grad(v)) * dx, (n, n)
    if kernel == "load":
        return SpatialCoordinate(2)[0] * v * dx, (n,)
    if kernel == "det":
        return abs(JacobianDeterminant(domain)) * u * v * dx, (n, n)
    if kernel == "nested":
        e = SpatialCoordinate(2)[0]
        for _ in range(4):
            e = e * e - e
        return e * v * dx, (n,)
    raise ValueError(kernel)


def test_shared_subexpression_is_computed_once():
    """In SSA form, an expression used twice is computed once."""
    s = ArrayEntry("T", (0,)) + ArrayEntry("T", (1,))
    generator = MLIRGenerator({"T": np.array([1.0, 2.0])})
    function = generator.function("f", [], {"T": np.array([1.0, 2.0])}, s * s, True)
    assert function.count("arith.addf") == 1
    assert function.count("arith.mulf") == 1


def test_unsupported_node_raises():
    """A node without an MLIR handler raises a NotImplementedError naming its type."""
    generator = MLIRGenerator({})
    with pytest.raises(NotImplementedError, match="MLIR generation is not implemented"):
        generator.function("f", [], {}, SpatialCoordinate(2)[0], True)


def test_geometry_function_is_called_once(lagrange_element):
    """|det J| in the integrand and in the change of variables is one call."""
    form, _ = _form(lagrange_element, "det", 1)
    module, _ = uflx_codegeneration.generate(form, language="MLIR")
    kernel = module.split("@tabulate_tensor_f64(")[1]
    assert len(re.findall(r"func\.call @geo\d+\(", kernel)) == 1


@requires_mlir
@pytest.mark.parametrize(
    ("kernel", "degree"),
    [
        *[(k, d) for k in ("mass", "stiffness", "load") for d in (1, 2, 3)],
        ("det", 1),
        ("nested", 1),
    ],
)
def test_mlir_and_c_kernels_agree(lagrange_element, code_dir, tmp_path, kernel, degree):
    """The MLIR kernel computes the same element tensor as the C kernel."""
    form, shape = _form(lagrange_element, kernel, degree)
    module, _ = uflx_codegeneration.generate(form, language="MLIR")
    code, signature = uflx_codegeneration.generate(form)

    result = _tabulate_mlir(module, str(tmp_path), shape)
    expected = _tabulate_c(code, signature, f"test_mlir_c_{kernel}_{degree}", code_dir, shape)
    # The C compiler may contract a * b + c into a fused multiply-add where llc does not.
    np.testing.assert_allclose(result, expected, rtol=1e-12, atol=1e-14 * np.max(np.abs(expected)))
