# Copyright (C) 2026 Jack S. Hale
#
# This file is part of FFCy (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    MIT

cmake_minimum_required(VERSION 3.24)
# The driver scatters with CUDA kernels, and nvcc supports this host compiler.
set(CMAKE_CUDA_COMPILER /usr/bin/nvcc CACHE FILEPATH "CUDA compiler")
set(CMAKE_CUDA_HOST_COMPILER /usr/bin/g++-12 CACHE FILEPATH "Host compiler for nvcc")
set(CMAKE_CUDA_ARCHITECTURES 70 CACHE STRING "CUDA architectures")
project(ffcy_benchmarks LANGUAGES CXX CUDA)

# evector acts on E-vectors. assembled and fused act on assembled vectors through a
# dofmap, scattering in the driver or with atomic adds in the action.
set(FFCY_VECTORS evector CACHE STRING "Vectors to act on")
set(FFCY_CELLS_PER_BLOCK 1 2 4 8 16 32 64 CACHE STRING "Cells per block to build")
set(CMAKE_CUDA_STANDARD 20)
include(../common.cmake)

find_package(CUDAToolkit REQUIRED)
add_plugin(ffcy_gpu_plugin ${FFCY_MLIR}/gpu FFCyGpuPlugin FFCY_GPU_PLUGIN)
set(FFCY_GPU_CHIP sm_70 CACHE STRING "GPU to compile for")
set(FFCY_SCHEDULE ${FFCY_MLIR}/gpu/schedule.mlir)

function(add_benchmark form degree variant cells)
  set(name ${form}_${FFCY_CELL}_q${degree}_${variant}_c${cells})
  set(base "${CMAKE_CURRENT_BINARY_DIR}/${name}")
  emit_action(${base} ${form} ${degree} ${variant} --cells-per-block ${cells})
  run_pipeline(
    ${base}
    "ffcy-gpu-pipeline{blocks=${FFCY_BLOCKS} schedule=${FFCY_SCHEDULE} chip=${FFCY_GPU_CHIP}}"
    ffcy_gpu_plugin ${FFCY_GPU_PLUGIN} ${FFCY_SCHEDULE})
  add_custom_command(
    OUTPUT ${base}.o
    COMMAND ${MLIR_ROOT}/bin/mlir-translate --mlir-to-llvmir ${base}.llvm.mlir -o ${base}.ll
    COMMAND ${MLIR_ROOT}/bin/llc -O3 -filetype=obj -relocation-model=pic ${base}.ll -o ${base}.o
    DEPENDS ${base}.llvm.mlir
    VERBATIM)
  set_source_files_properties(${base}.o PROPERTIES EXTERNAL_OBJECT TRUE GENERATED TRUE)

  # The driver loads each kernel as a library.
  add_library(${name} SHARED ${base}.o)
  set_target_properties(${name} PROPERTIES LINKER_LANGUAGE CXX BUILD_RPATH ${MLIR_ROOT}/lib)
  target_link_libraries(${name} PRIVATE ${MLIR_ROOT}/lib/libmlir_cuda_runtime.so)
endfunction()

add_executable(ffcy_driver driver.cu)
target_link_libraries(ffcy_driver PRIVATE CUDA::cudart ${CMAKE_DL_LIBS})

foreach(form IN LISTS FFCY_FORMS)
  foreach(degree IN LISTS FFCY_DEGREES)
    # Each cell takes n^(d - 1) threads, and a block holds at most 1024.
    math(EXPR max_cells "1024 / (${degree} + 1)")
    if(FFCY_DIM EQUAL 3)
      math(EXPR max_cells "${max_cells} / (${degree} + 1)")
    endif()
    foreach(geometry IN LISTS FFCY_GEOMETRIES)
      foreach(geometry_degree IN LISTS FFCY_GEOMETRY_DEGREES)
        geometry_variant(${geometry} ${geometry_degree} base)
        foreach(vectors IN LISTS FFCY_VECTORS)
          set(variant ${base})
          if(NOT vectors STREQUAL "evector")
            set(variant ${base}_${vectors})
          endif()
          add_case(${form} ${degree} ${variant})
          foreach(cells IN LISTS FFCY_CELLS_PER_BLOCK)
            if(cells LESS_EQUAL max_cells)
              add_benchmark(${form} ${degree} ${variant} ${cells})
            endif()
          endforeach()
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

add_cases_target()
