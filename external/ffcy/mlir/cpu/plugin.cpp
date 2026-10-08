// Copyright (C) 2026 Jack S. Hale
//
// This file is part of FFCy (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    MIT

// Registers the passes for CPU targets.

#include "passes.h"

#include "mlir/Tools/Plugins/PassPlugin.h"

extern "C" LLVM_ATTRIBUTE_WEAK mlir::PassPluginLibraryInfo mlirGetPassPluginInfo()
{
  return {MLIR_PLUGIN_API_VERSION, "FFCyCpuPasses", LLVM_VERSION_STRING,
          []()
          {
            register_slice_constant_indices();
            register_tile_stage_pencils();
            register_cpu_pipeline();
          }};
}
