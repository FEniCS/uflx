// Registers the passes shared by every target.

#include "passes.h"

#include "mlir/Tools/Plugins/PassPlugin.h"

extern "C" LLVM_ATTRIBUTE_WEAK mlir::PassPluginLibraryInfo mlirGetPassPluginInfo()
{
  return {MLIR_PLUGIN_API_VERSION, "FFCySharedPasses", LLVM_VERSION_STRING,
          []()
          {
            register_fuse_elementwise();
            register_sort_stages();
            register_shared_pipelines();
          }};
}
