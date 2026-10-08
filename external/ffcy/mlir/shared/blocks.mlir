// Copyright (C) 2026 Jack S. Hale
//
// This file is part of FFCy (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    MIT

// Tiles the cells into blocks and fuses each block's whole computation, for every
// target. Each block handles the cells of one block from ffcy-split-cells.
module attributes {transform.with_named_sequence} {
  transform.named_sequence @__blocks(%root: !transform.any_op {transform.readonly}) {
    %functions = transform.structured.match ops{["func.func"]} in %root
      : (!transform.any_op) -> !transform.any_op
    transform.foreach %functions : !transform.any_op {
    ^bb0(%function: !transform.any_op):
      transform.annotate %function "llvm.emit_c_interface" : !transform.any_op
      %return = transform.structured.match ops{["func.return"]} in %function
        : (!transform.any_op) -> !transform.any_op
      %last = transform.get_producer_of_operand %return[0] : (!transform.any_op) -> !transform.any_op
      %block_tile, %blocks = transform.structured.fuse %last tile_sizes [1] {use_forall, apply_cleanup}
        : (!transform.any_op) -> (!transform.any_op, !transform.any_op)
      transform.annotate %blocks "ffcy.block" : !transform.any_op
      // Fusion clones a producer for each of its uses, so the clones are merged.
      transform.apply_cse to %function : !transform.any_op

      // Slices of full-size empty outputs become per-block empty tensors.
      // Canonicalisation hoists constants, so the tables are cloned in afterwards.
      transform.apply_patterns to %function {
        transform.apply_patterns.canonicalization
        transform.apply_patterns.tensor.fold_tensor_empty
      } : !transform.any_op
      %tables = transform.structured.match attributes{ffcy.table} in %function
        : (!transform.any_op) -> !transform.any_op
      %fused_tables, %blocks_with_tables = transform.structured.fuse_into_containing_op %tables into %blocks
        : (!transform.any_op, !transform.any_op) -> (!transform.any_op, !transform.any_op)
    }
    transform.yield
  }
}
