import shutil

import pytest
from basix import CellType
from xdsl.dialects.builtin import ModuleOp
from xdsl.transforms.mlir_opt import DEFAULT_MLIR_OPT_EXECUTABLE

from ffcy.pipeline import compile_form, run_mlir_opt
from forms.mass import mass_form

pytestmark = pytest.mark.skipif(
    shutil.which(DEFAULT_MLIR_OPT_EXECUTABLE) is None,
    reason="needs mlir-opt on PATH or in XDSL_MLIR_OPT",
)


def test_pipeline_runs():
    module = compile_form(mass_form(2, CellType.hexahedron), 5)
    text = run_mlir_opt(module, "canonicalize,cse")
    assert "linalg.generic" in text


def test_failure_reports_pipeline():
    with pytest.raises(RuntimeError, match="no-such-pass"):
        run_mlir_opt(ModuleOp([]), "no-such-pass")
