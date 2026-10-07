"""Code generation."""

import quadraturerules
from uflx.algorithms import pull_back_to_reference, simplify
from uflx.geometry import (
    expand_geometry,
)
from uflx.graphs import (
    GraphNode,
)
from uflx.integrals import AbstractMeasure, dx
from uflx.maps import apply_push_forwards

from uflx_codegeneration import symbols
from uflx_codegeneration.algorithms import (
    expand_inner_products,
    insert_geometry_functions,
    tabulate_finite_elements,
)
from uflx_codegeneration.c import CGenerator, tables_to_c
from uflx_codegeneration.mlir import MLIRGenerator
from uflx_codegeneration.quadrature import (
    QuadratureRule,
    integrals_to_quadrature,
    quadrature_rule,
    tabulate_quadrature,
)
from uflx_codegeneration.utils import indented


def generate(
    form: GraphNode,
    language: str = "C",
    cse: bool = True,
) -> tuple[str, str]:
    """Generate code.

    Args:
        form: The form or other object to be assembled
        language: The language to generate: "C", or "MLIR" (func, scf, arith, math and memref
            dialects)
        cse: Evaluate each common subexpression of a statement once, into a temporary. This
            only affects C: in MLIR's SSA form every value is computed once regardless.

    Returns:
        The code, and the signature of the kernel (a C declaration, or an MLIR function type)
    """
    if language not in ("C", "MLIR"):
        raise NotImplementedError(f"Generation of {language} is not supported")

    # TODO: get this from somewhere
    rules: dict[AbstractMeasure, QuadratureRule] = {}
    # For now, use a degree 10 rule:
    points, weights = quadraturerules.single_integral_quadrature(
        quadraturerules.QuadratureRule.XiaoGimbutas,
        quadraturerules.Domain.Triangle,
        10,
    )
    rules[dx] = quadrature_rule([p[1:] for p in points], 0.5 * weights)

    # Apply algorithms from UFLx
    form = pull_back_to_reference(form)
    form = apply_push_forwards(form)
    form = simplify(form)

    # Apply codegeneration algorithms
    form = integrals_to_quadrature(form, rules)
    geometry_functions, form = insert_geometry_functions(form)
    form = expand_geometry(form)
    form = expand_inner_products(form)

    # Tabulate quadrature rules and finite element functions
    q_tables, form = tabulate_quadrature(form)
    fe_tables, form = tabulate_finite_elements(form)
    tables = {**q_tables, **fe_tables}

    # Tabulate the finite elements used by each geometry function
    lowered_geometry_functions = {}
    for fname, (dtype, inputs, function) in geometry_functions.items():
        ftables, function = tabulate_finite_elements(function)
        lowered_geometry_functions[fname] = (dtype, inputs, ftables, function)

    if language == "MLIR":
        all_tables = {**tables}
        for _, _, ftables, _ in lowered_geometry_functions.values():
            all_tables.update(ftables)
        mlir = MLIRGenerator(all_tables).module(lowered_geometry_functions, tables, form)
        mlir_signature = (
            "(memref<?xf64>, memref<?xf64>, memref<?xf64>, memref<?xf64>, memref<?xi32>, "
            "memref<?xi8>, !llvm.ptr) -> ()"
        )
        return mlir, mlir_signature

    generator = CGenerator(cse)
    code = ""
    for fname, (dtype, inputs, ftables, function) in lowered_geometry_functions.items():
        code += f"{dtype} {fname}("
        code += ", ".join(f"{i._dtype} {i._variable}" for i in inputs)
        code += ") {\n"
        code += indented(tables_to_c(ftables), 2)
        code += "\n\n"
        declarations, value = generator.statement(function)
        code += "".join(f"  {d}\n" for d in declarations)
        code += f"  return {value};\n"
        code += "}\n\n"
    code += (
        "void tabulate_tensor_f64(\n"
        f"    double* restrict {symbols.local_tensor},\n"
        f"    const double* restrict {symbols.coefficients},\n"
        f"    const double* restrict {symbols.constants},\n"
        f"    const double* restrict {symbols.coordinate_dofs},\n"
        f"    const int* restrict {symbols.entity_local_index},\n"
        f"    const uint8_t* restrict {symbols.quadrature_permutation},\n"
        f"    void* {symbols.custom_data}\n"
        ") {\n"
    )

    code += indented(tables_to_c(tables), 2)
    code += "\n\n"
    code += indented(generator.code(form), 2)
    code += "\n}\n"

    signature = (
        "void tabulate_tensor_f64(double* restrict, const double* restrict, "
        "const double* restrict, const double* restrict, const int* restrict, "
        "const uint8_t* restrict, void*);"
    )

    return code, signature
