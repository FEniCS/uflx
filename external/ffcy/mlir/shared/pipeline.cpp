// Copyright (C) 2026 Jack S. Hale
//
// This file is part of FFCy (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    MIT

// Registers the pipeline steps shared by every target.

#include "passes.h"

#include "mlir/Pass/PassRegistry.h"

using namespace mlir;

void register_shared_pipelines()
{
  // Size queries on intermediates would block fusion, so they are resolved first.
  PassPipelineRegistration<>(
      "ffcy-fuse", "Fuse elementwise ops once size queries are resolved",
      [](OpPassManager& pm)
      {
        if (failed(parsePassPipeline(
                "func.func(resolve-shaped-type-result-dims,canonicalize,cse,ffcy-fuse-elementwise)",
                pm)))
          llvm::report_fatal_error("invalid ffcy-fuse pipeline");
      });
}
