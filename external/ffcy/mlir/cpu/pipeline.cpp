// Copyright (C) 2026 Jack S. Hale
//
// This file is part of FFCy (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    MIT

// Registers ffcy-cpu-pipeline, which compiles an emitted action on interleaved cells
// to the LLVM dialect.
//
// The pipeline fuses elementwise ops, tiles the cell blocks with @__blocks, fuses
// each stage into one loop over pencils, vectorises across the cells of a block,
// once inputs read at constant indices are read through slices instead, then
// bufferizes with each block's intermediates allocated once on the stack and reused
// by every block, and lowers. Contractions become outer products, which are FMAs
// across lanes.

#include "passes.h"

#include "mlir/Pass/PassOptions.h"
#include "mlir/Pass/PassRegistry.h"

#include <string>

using namespace mlir;

namespace
{
struct Options : PassPipelineOptions<Options>
{
  Option<std::string> blocks{*this, "blocks",
                             llvm::cl::desc("Transform library with @__blocks")};
  Option<std::string> schedule{
      *this, "schedule", llvm::cl::desc("Transform library with @__vectorize and @__bufferize")};
  Option<int64_t> stack_bytes{*this, "stack-bytes",
                              llvm::cl::desc("Largest intermediate on the stack"),
                              llvm::cl::init(1048576)};
};

void build(OpPassManager& pm, const Options& options)
{
  const std::string pipeline
      = "ffcy-fuse,transform-preload-library{transform-library-paths=" + options.blocks + ","
        + options.schedule
        + "},transform-interpreter{entry-point=__blocks},ffcy-sort-stages,"
          "ffcy-tile-stage-pencils,ffcy-slice-constant-indices,"
          "transform-interpreter{entry-point=__vectorize},"
          "transform-interpreter{entry-point=__bufferize},"
          // Stage loops are scf.forall until they are bufferized.
          "scf-forall-to-for,canonicalize,cse,"
          "buffer-results-to-out-params{hoist-dynamic-allocs=true modify-public-functions=true},"
          "resolve-ranked-shaped-type-result-dims,canonicalize,"
          "func.func(buffer-loop-hoisting,promote-buffers-to-stack{max-alloc-size-in-bytes="
        + std::to_string(options.stack_bytes)
        + "}),convert-vector-to-scf{full-unroll=true},convert-linalg-to-loops,"
          "fold-memref-alias-ops,expand-strided-metadata,lower-affine,convert-scf-to-cf,"
          "convert-vector-to-llvm{vector-contract-lowering=outerproduct},convert-to-llvm,"
          "reconcile-unrealized-casts";
  if (failed(parsePassPipeline(pipeline, pm)))
    llvm::report_fatal_error("invalid ffcy-cpu-pipeline");
}
} // namespace

void register_cpu_pipeline()
{
  PassPipelineRegistration<Options>(
      "ffcy-cpu-pipeline", "Compile an emitted action to the LLVM dialect for CPUs", build);
}
