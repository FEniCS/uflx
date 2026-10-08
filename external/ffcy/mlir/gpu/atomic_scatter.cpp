// Copyright (C) 2026 Jack S. Hale
//
// This file is part of FFCy (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    MIT

// Sums an action's values on cells into an assembled vector with atomic adds.
//
// ffcy-assemble tags an action whose values are to be scattered with
// ffcy.scatter_dofmap, the index of its dofmap argument. After bufferization the
// values are the last argument, written through views in each thread's loop of
// the last stage, possibly many times as a contraction accumulates. Each thread
// instead writes its part into a private buffer, then adds it into the assembled
// vector at the dofs the same view of the dofmap gives. The argument becomes that
// vector, which the caller zeroes.

#include "passes.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/Linalg/IR/Linalg.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/IR/Builders.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace
{
constexpr StringLiteral scatter_dofmap = "ffcy.scatter_dofmap";

// The same chain of subviews as `view` has from `root`, applied to `source`.
Value rebuild_view(OpBuilder& builder, Value view, Value root, Value source)
{
  if (view == root)
    return source;
  auto subview = cast<memref::SubViewOp>(view.getDefiningOp());
  Value inner = rebuild_view(builder, subview.getSource(), root, source);
  auto type = memref::SubViewOp::inferRankReducedResultType(
      subview.getType().getShape(), cast<MemRefType>(inner.getType()),
      subview.getMixedOffsets(), subview.getMixedSizes(), subview.getMixedStrides());
  return memref::SubViewOp::create(builder, subview.getLoc(), cast<MemRefType>(type), inner,
                                   subview.getMixedOffsets(), subview.getMixedSizes(),
                                   subview.getMixedStrides());
}

LogicalResult scatter_atomically(func::FuncOp function)
{
  auto dofmap_index = function->getAttrOfType<IntegerAttr>(scatter_dofmap);
  Block& entry = function.getBody().front();
  BlockArgument values = entry.getArguments().back();
  BlockArgument dofmap = entry.getArgument(dofmap_index.getInt());
  auto values_type = cast<MemRefType>(values.getType());
  auto dofmap_type = cast<MemRefType>(dofmap.getType());
  MLIRContext* context = function.getContext();
  OpBuilder builder(context);
  Location loc = function.getLoc();

  // The dofmap has the values' cell axes, then the flat index of a dof on its cell.
  const int64_t cell_axes = dofmap_type.getRank() - 1;
  SmallVector<ReassociationIndices> reassociation;
  for (int64_t a = 0; a < cell_axes; ++a)
    reassociation.push_back({a});
  reassociation.emplace_back();
  for (int64_t a = cell_axes; a < values_type.getRank(); ++a)
    reassociation.back().push_back(a);
  builder.setInsertionPointToStart(&entry);
  SmallVector<OpFoldResult> output_shape;
  for (int64_t a = 0; a < values_type.getRank(); ++a)
  {
    if (values_type.isDynamicDim(a))
      output_shape.push_back(memref::DimOp::create(builder, loc, dofmap, a).getResult());
    else
      output_shape.push_back(builder.getIndexAttr(values_type.getDimSize(a)));
  }
  Value dofs = memref::ExpandShapeOp::create(
      builder, loc, values_type.clone(dofmap_type.getElementType()), dofmap, reassociation,
      output_shape);
  BlockArgument vector = entry.addArgument(
      MemRefType::get({ShapedType::kDynamic}, values_type.getElementType()), loc);

  // Views of the values that ops other than subviews write through.
  SmallVector<Value> written, worklist{values};
  while (!worklist.empty())
  {
    Value view = worklist.pop_back_val();
    for (Operation* user : view.getUsers())
    {
      if (auto subview = dyn_cast<memref::SubViewOp>(user))
        worklist.push_back(subview);
      else if (!llvm::is_contained(written, view))
        written.push_back(view);
    }
  }

  auto private_space = gpu::AddressSpaceAttr::get(context, gpu::AddressSpace::Private);
  for (Value view : written)
  {
    auto view_type = cast<MemRefType>(view.getType());
    Block* block = view.getParentBlock();
    auto loop = dyn_cast<scf::ForallOp>(block->getParentOp());
    if (!view_type.hasStaticShape() or !loop or !loop->hasAttr("ffcy.threads"))
      return view.getDefiningOp()->emitError(
          "scattered values must be written by a thread loop through static views");

    builder.setInsertionPointAfterValue(view);
    Value buffer = memref::AllocaOp::create(
        builder, loc,
        MemRefType::get(view_type.getShape(), view_type.getElementType(), nullptr,
                        private_space));
    view.replaceAllUsesWith(buffer);

    builder.setInsertionPoint(block->getTerminator());
    Value dof_view = rebuild_view(builder, view, values, dofs);
    // A linalg op rather than loops keeps its bounds constant, since loops are only
    // made inside the outlined kernel.
    const int64_t rank = view_type.getRank();
    linalg::GenericOp::create(
        builder, loc, TypeRange{}, ValueRange{dof_view}, ValueRange{buffer},
        SmallVector<AffineMap>(2, builder.getMultiDimIdentityMap(rank)),
        SmallVector<utils::IteratorType>(rank, utils::IteratorType::parallel),
        [&](OpBuilder& b, Location l, ValueRange args)
        {
          Value dof = arith::IndexCastOp::create(b, l, b.getIndexType(), args[0]);
          memref::AtomicRMWOp::create(b, l, arith::AtomicRMWKind::addf, args[1], vector,
                                      ValueRange{dof});
          linalg::YieldOp::create(b, l, args[1]);
        });
  }

  // The views of the values are now unused.
  for (bool erased = true; erased;)
  {
    erased = false;
    function.walk(
        [&](memref::SubViewOp subview)
        {
          if (subview->use_empty())
          {
            subview->erase();
            erased = true;
          }
        });
  }
  if (!values.use_empty())
    return function.emitError("scattered values have uses other than views");
  entry.eraseArgument(values.getArgNumber());
  function.setType(
      FunctionType::get(context, entry.getArgumentTypes(), function.getResultTypes()));
  function->removeAttr(scatter_dofmap);
  return success();
}

struct AtomicScatter : PassWrapper<AtomicScatter, OperationPass<ModuleOp>>
{
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(AtomicScatter)

  StringRef getArgument() const final { return "ffcy-atomic-scatter"; }

  StringRef getDescription() const final
  {
    return "Sum values on cells into an assembled vector with atomic adds";
  }

  void getDependentDialects(DialectRegistry& registry) const override
  {
    registry.insert<arith::ArithDialect, gpu::GPUDialect, linalg::LinalgDialect,
                    memref::MemRefDialect>();
  }

  void runOnOperation() override
  {
    getOperation().walk(
        [&](func::FuncOp function)
        {
          if (function->hasAttr(scatter_dofmap) and failed(scatter_atomically(function)))
            signalPassFailure();
        });
  }
};
} // namespace

void register_atomic_scatter() { PassRegistration<AtomicScatter>(); }
