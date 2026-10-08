# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

# Options, plugins and pipeline steps shared by the benchmarks of every target.
#
# Configure from a shell where the Spack environment and the venv are active, so
# that the compiler finds UFLx, xDSL and Basix.
set(MLIR_ROOT "/home/jhale/llvm-mlir-22.1.2/install" CACHE PATH "MLIR installation")
find_program(FFCY_PYTHON python3 HINTS "$ENV{VIRTUAL_ENV}/bin" NO_DEFAULT_PATH)
find_program(FFCY_PYTHON python3 REQUIRED)
set(FFCY_PYTHONPATH "$ENV{PYTHONPATH}" CACHE STRING "PYTHONPATH for FFCy")
set(FFCY_FORMS mass laplace helmholtz CACHE STRING "Forms to build")
set(FFCY_CELL hexahedron CACHE STRING "quadrilateral or hexahedron")
set(FFCY_DEGREES 1 2 3 4 5 6 7 CACHE STRING "Polynomial degrees to build")
# inline computes the geometry in the action, precomputed reads it as qdata.
set(FFCY_GEOMETRIES precomputed CACHE STRING "Geometry evaluations to build")
# Degrees of the geometry describing the cells, which keep straight sides. Variants
# with geometry of degree k above 1 have _gk after the geometry evaluation.
set(FFCY_GEOMETRY_DEGREES 1 CACHE STRING "Geometry degrees to build")
option(FFCY_DUMP_IR "Write the IR after every pass of a kernel's pipeline" OFF)
# The plugins must match the compiler that built MLIR.
set(FFCY_PLUGIN_CXX_COMPILER /usr/bin/c++ CACHE FILEPATH "Compiler for mlir-opt plugins")
set(CMAKE_CXX_STANDARD 20)
# The drivers do work on the host that is timed, such as the CPU scatter.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(FFCY_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
set(FFCY_MLIR "${FFCY_ROOT}/mlir")
set(FFCY_BLOCKS "${FFCY_MLIR}/shared/blocks.mlir")
file(GLOB_RECURSE FFCY_SOURCES CONFIGURE_DEPENDS
     "${FFCY_ROOT}/ffcy/*.py" "${FFCY_ROOT}/forms/*.py" "${FFCY_ROOT}/bench/*.py")

if(FFCY_CELL STREQUAL "hexahedron")
  set(FFCY_DIM 3)
else()
  set(FFCY_DIM 2)
endif()

include(ExternalProject)

# Builds the plugin in `source` as `target`, setting `library` to its path.
function(add_plugin target source library_name library)
  ExternalProject_Add(
    ${target}
    SOURCE_DIR ${source}
    CMAKE_ARGS -DMLIR_DIR=${MLIR_ROOT}/lib/cmake/mlir
               -DCMAKE_CXX_COMPILER=${FFCY_PLUGIN_CXX_COMPILER}
               -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
    BUILD_BYPRODUCTS <BINARY_DIR>/lib${library_name}.so
    # Rebuild when the plugin source changes.
    BUILD_ALWAYS ON
    INSTALL_COMMAND "")
  ExternalProject_Get_Property(${target} BINARY_DIR)
  set(${library} ${BINARY_DIR}/lib${library_name}.so PARENT_SCOPE)
endfunction()

add_plugin(ffcy_shared_plugin ${FFCY_MLIR}/shared FFCySharedPlugin
           FFCY_SHARED_PLUGIN)

# A variant is a geometry evaluation, then _gk for geometry of degree k above 1, then
# _assembled or _fused for actions on assembled vectors.
function(variant_flags variant flags)
  set(result)
  if(variant MATCHES "^precomputed")
    list(APPEND result --precompute-geometry)
  endif()
  if(variant MATCHES "_g([0-9]+)")
    list(APPEND result --geometry-degree ${CMAKE_MATCH_1})
  endif()
  if(variant MATCHES "_assembled$")
    list(APPEND result --assemble)
  elseif(variant MATCHES "_fused$")
    list(APPEND result --assemble --atomic-scatter)
  endif()
  set(${flags} ${result} PARENT_SCOPE)
endfunction()

# The arguments and expected result of the action for a few cells, shared by every
# block size of the same form, degree and variant.
function(add_case form degree variant)
  set(case "${CMAKE_CURRENT_BINARY_DIR}/cases/${form}_${FFCY_CELL}_q${degree}_${variant}")
  variant_flags(${variant} emit_flags)
  add_custom_command(
    OUTPUT ${case}.case ${case}.bin
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/cases
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${FFCY_PYTHONPATH} ${FFCY_PYTHON} -m bench.cases
            --form ${form} --cell ${FFCY_CELL} --degree ${degree} ${emit_flags} ${case}
    WORKING_DIRECTORY ${FFCY_ROOT}
    DEPENDS ${FFCY_SOURCES}
    VERBATIM)
  set_property(GLOBAL APPEND PROPERTY FFCY_CASES ${case}.case)
endfunction()

# Emits the action as base.mlir.
function(emit_action base form degree variant)
  variant_flags(${variant} emit_flags)
  add_custom_command(
    OUTPUT ${base}.mlir
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${FFCY_PYTHONPATH} ${FFCY_PYTHON} -m bench.emit
            --form ${form} --cell ${FFCY_CELL} --degree ${degree} ${emit_flags} ${ARGN}
            -o ${base}.mlir
    WORKING_DIRECTORY ${FFCY_ROOT}
    DEPENDS ${FFCY_SOURCES}
    VERBATIM)
endfunction()

# Compiles base.mlir to base.llvm.mlir with `pipeline` in one mlir-opt run, loading
# the shared plugin and `plugin`, built by `plugin_target`. With FFCY_DUMP_IR the IR
# after every pass goes to base.ir.
function(run_pipeline base pipeline plugin_target plugin)
  set(dump)
  if(FFCY_DUMP_IR)
    set(dump --mlir-print-ir-after-all --mlir-print-ir-tree-dir=${base}.ir)
  endif()
  add_custom_command(
    OUTPUT ${base}.llvm.mlir
    COMMAND ${MLIR_ROOT}/bin/mlir-opt ${base}.mlir --load-pass-plugin=${FFCY_SHARED_PLUGIN}
            --load-pass-plugin=${plugin} "--pass-pipeline=builtin.module(${pipeline})" ${dump}
            -o ${base}.llvm.mlir
    DEPENDS ${base}.mlir ${FFCY_BLOCKS} ${ARGN} ffcy_shared_plugin ${FFCY_SHARED_PLUGIN}
            ${plugin_target} ${plugin}
    VERBATIM)
endfunction()

# The variant name of a geometry evaluation with geometry of `degree`.
function(geometry_variant geometry degree variant)
  if(degree EQUAL 1)
    set(${variant} ${geometry} PARENT_SCOPE)
  else()
    set(${variant} ${geometry}_g${degree} PARENT_SCOPE)
  endif()
endfunction()

function(add_cases_target)
  get_property(cases GLOBAL PROPERTY FFCY_CASES)
  add_custom_target(cases ALL DEPENDS ${cases})
endfunction()
