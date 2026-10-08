// The GPU parts of the staged pencil schedule that do not depend on the form.
//
// @__bufferize bufferizes and prepares the block loops from @__blocks in
// mlir/shared/blocks.mlir for mapping. Between them ffcy-sort-stages makes each stage
// contiguous and ffcy-tile-stages tiles each stage onto threads, and after them
// ffcy-map-threads maps the thread loops.
//
// Stage outputs live in shared memory, and the rest of a stage stays in each
// thread's private memory.
module attributes {transform.with_named_sequence} {
  transform.named_sequence @__bufferize(%root: !transform.any_op {transform.consumed}) {
    // Intermediates a thread computes and uses itself stay in its registers.
    %thread_loops = transform.structured.match attributes{ffcy.threads} in %root
      : (!transform.any_op) -> !transform.any_op
    transform.foreach %thread_loops : !transform.any_op {
    ^bb0(%loop: !transform.any_op):
      %empties = transform.structured.match ops{["tensor.empty"]} in %loop
        : (!transform.any_op) -> !transform.any_op
      // Only ops writing into an empty tensor can be given a buffer.
      %privates = transform.get_consumers_of_result %empties[0]
        : (!transform.any_op) -> !transform.any_op
      %buffers, %new = transform.structured.bufferize_to_allocation %privates
        {memory_space = #gpu.address_space<private>, alloc_op = "memref.alloca",
         bufferize_destination_only}
        : !transform.any_op
    }

    // GPU code cannot call the runtime's memref copy, so copies are linalg.copy.
    %bufferized = transform.bufferization.one_shot_bufferize %root {
      bufferize_function_boundaries = true,
      function_boundary_type_conversion = 1 : i32,
      memcpy_op = "linalg.copy"
    } : (!transform.any_op) -> !transform.any_op

    // The block loop becomes scf.parallel, which can be mapped with a dynamic grid.
    %block_loops = transform.structured.match attributes{ffcy.block} in %bufferized
      : (!transform.any_op) -> !transform.any_op
    transform.foreach %block_loops : !transform.any_op {
    ^bb0(%loop: !transform.any_op):
      %parallel = transform.loop.forall_to_parallel %loop : (!transform.any_op) -> !transform.any_op
    }
    transform.yield
  }
}
