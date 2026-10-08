# An mlir-opt pass plugin. mlir-opt provides the MLIR symbols when it loads it.
find_package(MLIR REQUIRED CONFIG)

function(add_ffcy_plugin name)
  add_library(${name} MODULE ${ARGN})
  target_include_directories(${name} PRIVATE ${LLVM_INCLUDE_DIRS} ${MLIR_INCLUDE_DIRS})
  target_compile_definitions(${name} PRIVATE ${LLVM_DEFINITIONS})
  target_compile_features(${name} PRIVATE cxx_std_17)
  # LLVM is built without RTTI.
  target_compile_options(${name} PRIVATE -fno-rtti)
endfunction()
