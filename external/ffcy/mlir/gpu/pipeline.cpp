// Registers ffcy-gpu-pipeline, which compiles an emitted action to the LLVM dialect
// with a GPU binary, and ffcy-check-mapped, which it uses.
//
// The pipeline fuses elementwise ops, tiles the cell blocks with @__blocks, tiles
// each stage onto threads, bufferizes with @__bufferize, maps the block loops to a
// dynamic grid and the thread loops to a block, moves shared buffers to workgroup
// memory and private ones to the stack, outlines the kernel and lowers it.

#include "passes.h"

#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"
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
  Option<std::string> schedule{*this, "schedule",
                               llvm::cl::desc("Transform library with @__bufferize")};
  Option<std::string> chip{*this, "chip", llvm::cl::desc("GPU to compile for"),
                           llvm::cl::init("sm_70")};
};

void build(OpPassManager& pm, const Options& options)
{
  const std::string pipeline
      = "ffcy-fuse,transform-preload-library{transform-library-paths=" + options.blocks + ","
        + options.schedule
        + "},transform-interpreter{entry-point=__blocks},ffcy-sort-stages,ffcy-tile-stages,"
          "transform-interpreter{entry-point=__bufferize},"
          // The caller allocates the result, and loop bounds depend only on the arguments,
          // so that they are visible before the launch to be mapped.
          "canonicalize,cse,"
          "buffer-results-to-out-params{hoist-dynamic-allocs=true modify-public-functions=true},"
          "resolve-ranked-shaped-type-result-dims,canonicalize,"
          "func.func(gpu-map-parallel-loops,convert-parallel-loops-to-gpu),canonicalize,"
          // An action tagged for an atomic scatter adds its values into y while threads are
          // still loops. Allocations left after moving shared ones to workgroup memory are
          // each thread's own, so they go on its stack.
          "ffcy-atomic-scatter,ffcy-map-threads,ffcy-workgroup-allocations,"
          "func.func(promote-buffers-to-stack{max-alloc-size-in-bytes=1024}),cse,canonicalize,"
          "gpu-kernel-outlining,canonicalize,ffcy-check-mapped,"
          "convert-linalg-to-loops,fold-memref-alias-ops,expand-strided-metadata,lower-affine,"
          "convert-scf-to-cf,gpu-lower-to-nvvm-pipeline{cubin-chip="
        + options.chip + " cubin-format=fatbin opt-level=3}";
  if (failed(parsePassPipeline(pipeline, pm)))
    llvm::report_fatal_error("invalid ffcy-gpu-pipeline");
}

// A loop that could not be mapped would silently run on the host.
struct CheckMapped : PassWrapper<CheckMapped, OperationPass<ModuleOp>>
{
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(CheckMapped)

  StringRef getArgument() const final { return "ffcy-check-mapped"; }

  StringRef getDescription() const final
  {
    return "Fail if a parallel loop is left to run on the host";
  }

  void runOnOperation() override
  {
    getOperation().walk(
        [&](scf::ParallelOp loop)
        {
          loop.emitError("parallel loop left on the host");
          signalPassFailure();
        });
  }
};
} // namespace

void register_gpu_pipeline()
{
  PassRegistration<CheckMapped>();
  PassPipelineRegistration<Options>(
      "ffcy-gpu-pipeline", "Compile an emitted action to the LLVM dialect for GPUs", build);
}
