# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

cmake_minimum_required(VERSION 3.24)
project(ffcy_cpu_benchmarks LANGUAGES CXX)

set(FFCY_CELLS_PER_BLOCK 1 2 4 8 16 32 64 CACHE STRING "Cells per block to build")
# evector acts on E-vectors, assembled on assembled vectors through a dofmap, with
# the driver summing the values on cells into the result. Its variants end in
# _assembled.
set(FFCY_VECTORS evector CACHE STRING "Vectors to act on")
set(FFCY_CPU znver2 CACHE STRING "CPU to generate code for")
# Without a triple LLVM does not know the C library, and expands memset into a loop.
set(FFCY_TRIPLE x86_64-unknown-linux-gnu CACHE STRING "Target triple")
# A block's intermediates go on the stack, so they must fit in it.
set(FFCY_STACK_BYTES 1048576 CACHE STRING "Largest intermediate on the stack")
include(../common.cmake)

add_plugin(ffcy_cpu_plugin ${FFCY_MLIR}/cpu FFCyCpuPlugin FFCY_CPU_PLUGIN)
set(FFCY_SCHEDULE ${FFCY_MLIR}/cpu/schedule.mlir)

function(add_benchmark form degree variant cells)
  set(name ${form}_${FFCY_CELL}_q${degree}_${variant}_c${cells})
  set(base "${CMAKE_CURRENT_BINARY_DIR}/${name}")
  emit_action(${base} ${form} ${degree} ${variant} --lanes ${cells})
  run_pipeline(
    ${base}
    "ffcy-cpu-pipeline{blocks=${FFCY_BLOCKS} schedule=${FFCY_SCHEDULE} stack-bytes=${FFCY_STACK_BYTES}}"
    ffcy_cpu_plugin ${FFCY_CPU_PLUGIN} ${FFCY_SCHEDULE})
  add_custom_command(
    OUTPUT ${base}.o
    COMMAND ${MLIR_ROOT}/bin/mlir-translate --mlir-to-llvmir ${base}.llvm.mlir -o ${base}.ll
    COMMAND ${MLIR_ROOT}/bin/opt -O3 -mtriple=${FFCY_TRIPLE} -mcpu=${FFCY_CPU} ${base}.ll -o ${base}.opt.bc
    COMMAND ${MLIR_ROOT}/bin/llc -O3 -mtriple=${FFCY_TRIPLE} -mcpu=${FFCY_CPU} -filetype=obj -relocation-model=pic
            ${base}.opt.bc -o ${base}.o
    DEPENDS ${base}.llvm.mlir
    VERBATIM)
  set_source_files_properties(${base}.o PROPERTIES EXTERNAL_OBJECT TRUE GENERATED TRUE)

  # The driver loads each kernel as a library.
  add_library(${name} SHARED ${base}.o)
  set_target_properties(${name} PROPERTIES LINKER_LANGUAGE CXX)
endfunction()

find_package(OpenMP REQUIRED)
add_executable(ffcy_driver driver.cpp)
target_link_libraries(ffcy_driver PRIVATE OpenMP::OpenMP_CXX ${CMAKE_DL_LIBS})

foreach(form IN LISTS FFCY_FORMS)
  foreach(degree IN LISTS FFCY_DEGREES)
    foreach(geometry IN LISTS FFCY_GEOMETRIES)
      foreach(geometry_degree IN LISTS FFCY_GEOMETRY_DEGREES)
        geometry_variant(${geometry} ${geometry_degree} base)
        foreach(vectors IN LISTS FFCY_VECTORS)
          set(variant ${base})
          if(vectors STREQUAL "assembled")
            set(variant ${base}_assembled)
          endif()
          add_case(${form} ${degree} ${variant})
          foreach(cells IN LISTS FFCY_CELLS_PER_BLOCK)
            add_benchmark(${form} ${degree} ${variant} ${cells})
          endforeach()
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

add_cases_target()

# Checks of the CPU plugin's passes, with FileCheck from the MLIR installation.
enable_testing()
set(FFCY_PLUGIN_TESTS ${FFCY_MLIR}/cpu/test)
add_test(NAME tile_stage_pencils
         COMMAND sh -c "${MLIR_ROOT}/bin/mlir-opt $0 --load-pass-plugin=$1 --pass-pipeline='builtin.module(ffcy-tile-stage-pencils)' | ${MLIR_ROOT}/bin/FileCheck $0"
                 ${FFCY_PLUGIN_TESTS}/tile_stage_pencils.mlir ${FFCY_CPU_PLUGIN})
