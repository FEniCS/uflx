// The CPU schedule after @__blocks in mlir/shared/blocks.mlir, on cells interleaved by
// ffcy-interleave-cells. ffcy-sort-stages, ffcy-tile-stage-pencils and
// ffcy-slice-constant-indices run before @__vectorize. After @__bufferize,
// scf-forall-to-for makes the block and stage loops run in turn on one thread, so a
// caller can run blocks on other threads or ranks by passing a range of them.
module attributes {transform.with_named_sequence} {
  // ffcy-tile-stage-pencils has tiled every op of a block to one pencil of its lanes,
  // so each tile becomes vector ops whose innermost dimension holds a cell in each
  // lane.
  transform.named_sequence @__vectorize(%root: !transform.any_op {transform.readonly}) {
    %functions = transform.structured.match ops{["func.func"]} in %root
      : (!transform.any_op) -> !transform.any_op
    transform.foreach %functions : !transform.any_op {
    ^bb0(%function: !transform.any_op):
      // Gathers read the assembled vector at indices from the dofmap.
      %vectorized = transform.structured.vectorize_children_and_apply_patterns %function
        {vectorize_nd_extract}
        : (!transform.any_op) -> !transform.any_op
      // A tile is one pencil, so its vectors have unit dimensions for the other axes.
      transform.apply_patterns to %vectorized {
        transform.apply_patterns.vector.drop_unit_dims_with_shape_cast
        transform.apply_patterns.vector.cast_away_vector_leading_one_dim
        transform.apply_patterns.canonicalization
      } : !transform.any_op
    }
    transform.yield
  }

  transform.named_sequence @__bufferize(%root: !transform.any_op {transform.consumed}) {
    // Kernels are not linked with the runtime's memref copy, so copies are linalg.copy.
    %bufferized = transform.bufferization.one_shot_bufferize %root {
      bufferize_function_boundaries = true,
      function_boundary_type_conversion = 1 : i32,
      memcpy_op = "linalg.copy"
    } : (!transform.any_op) -> !transform.any_op
    transform.yield
  }
}
