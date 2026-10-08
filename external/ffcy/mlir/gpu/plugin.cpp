// Copyright (C) 2026 Jack S. Hale
//
// This file is part of FFCy (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    MIT

// Registers the passes for GPU targets.

#include "passes.h"

#include "mlir/Tools/Plugins/PassPlugin.h"

extern "C" LLVM_ATTRIBUTE_WEAK mlir::PassPluginLibraryInfo mlirGetPassPluginInfo()
{
  return {MLIR_PLUGIN_API_VERSION, "FFCyGpuPasses", LLVM_VERSION_STRING,
          []()
          {
            register_workgroup_allocations();
            register_atomic_scatter();
            register_pencil_schedule();
            register_gpu_pipeline();
          }};
}
