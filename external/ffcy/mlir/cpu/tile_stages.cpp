// Copyright (C) 2026 Jack S. Hale
//
// This file is part of FFCy (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    MIT

// Tiles each stage from ffcy-annotate-stages into one loop over pencils, on cells
// interleaved by ffcy-interleave-cells.
//
// The sinks of a stage are tiled to one pencil of their lanes and their loops fused,
// then the rest of the stage is fused into the loop, so a pencil's values pass from
// op to op without going through the block's buffers. Fusion drops a producer once
// its uses in the loop are fused, so producers are listed after their users, in
// reverse program order, for all their uses to be in the loop by then. The loops are
// unmapped scf.forall, which run in turn once bufferized.

#include "passes.h"

#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/Linalg/IR/Linalg.h"
#include "mlir/Dialect/Linalg/TransformOps/LinalgTransformOps.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/SCF/TransformOps/SCFTransformOps.h"
#include "mlir/Dialect/Tensor/TransformOps/TensorTransformOps.h"
#include "mlir/Dialect/Transform/IR/TransformDialect.h"
#include "mlir/Dialect/Transform/IR/TransformOps.h"
#include "mlir/Dialect/Transform/IR/TransformTypes.h"
#include "mlir/Dialect/Transform/Interfaces/TransformInterfaces.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"

#include <map>

using namespace mlir;

namespace
{
// A stage's sinks along each pencil axis, in order of first appearance, and the
// number of its other ops.
struct Stage
{
  SmallVector<std::pair<int64_t, int64_t>> sinks;
  int64_t members = 0;
};

std::map<int64_t, Stage> stages_of(Operation* blocks, int64_t& dim)
{
  std::map<int64_t, Stage> stages;
  blocks->walk(
      [&](Operation* op)
      {
        if (auto sink = op->getAttrOfType<IntegerAttr>("ffcy.sink"))
        {
          // Tensors are [blocks, ..., lanes] after ffcy-interleave-cells.
          dim = cast<ShapedType>(op->getResult(0).getType()).getRank() - 2;
          auto& sinks = stages[sink.getInt()].sinks;
          const int64_t axis = op->getAttrOfType<IntegerAttr>("ffcy.pencil").getInt();
          auto found = llvm::find_if(sinks, [&](auto& s) { return s.first == axis; });
          if (found == sinks.end())
            sinks.push_back({axis, 1});
          else
            ++found->second;
        }
        else if (auto stage = op->getAttrOfType<IntegerAttr>("ffcy.stage"))
          ++stages[stage.getInt()].members;
      });
  return stages;
}

struct Builder
{
  OpBuilder& b;
  Location loc;
  Type any = transform::AnyOpType::get(b.getContext());

  Value match(Value target, ArrayRef<StringRef> ops, ArrayRef<NamedAttribute> attrs)
  {
    return transform::MatchOp::create(b, loc, any, target,
                                      ops.empty() ? nullptr : b.getStrArrayAttr(ops), nullptr,
                                      attrs.empty() ? nullptr : b.getDictionaryAttr(attrs),
                                      nullptr, nullptr);
  }

  // `count` handles to the ops of `handle`, one each.
  SmallVector<Value> split(Value handle, int64_t count)
  {
    if (count == 1)
      return {handle};
    auto op = transform::SplitHandleOp::create(b, loc, handle, count);
    return SmallVector<Value>(op.getResults());
  }
};

void tile_stages(Builder& t, Value function, const std::map<int64_t, Stage>& stages,
                 int64_t dim)
{
  OpBuilder& b = t.b;
  auto attr = [&](StringRef name, int64_t value)
  { return b.getNamedAttr(name, b.getI64IntegerAttr(value)); };
  Value blocks = t.match(function, {}, {b.getNamedAttr("ffcy.block", b.getUnitAttr())});

  for (const auto& [stage, info] : stages)
  {
    if (info.sinks.empty())
      continue;
    SmallVector<Value> loops;
    for (auto [axis, count] : info.sinks)
    {
      // The block axis is already tiled, and the reductions and lanes come last.
      SmallVector<int64_t> sizes{0};
      for (int64_t a = 1; a < dim + 1; ++a)
        sizes.push_back(a == axis ? 0 : 1);
      Value sinks = t.match(blocks, {"linalg.generic"},
                            {attr("ffcy.sink", stage), attr("ffcy.pencil", axis)});
      auto tiled = transform::TileUsingForallOp::create(b, t.loc, sinks, ArrayRef(sizes),
                                                        transform::TileSizesSpec());
      llvm::append_range(loops, t.split(tiled.getForallOp(), count));
    }
    // A loop is fused into a later one, which its operands dominate.
    Value loop = loops.front();
    for (Value other : ArrayRef(loops).drop_front())
      loop = transform::LoopFuseSiblingOp::create(b, t.loc, t.any, loop, other);
    transform::AnnotateOp::create(b, t.loc, loop, "ffcy.pencils", Value());

    // A producer is fused once per use in the loop, so the copies are merged before
    // the next producer is fused, or they would multiply through shared producers.
    SmallVector<Value> members
        = t.split(t.match(blocks, {}, {attr("ffcy.stage", stage)}), info.members);
    for (Value producer : llvm::reverse(members))
    {
      loop = transform::FuseIntoContainingOp::create(b, t.loc, producer, loop)
                 .getNewContainingOp();
      transform::ApplyCommonSubexpressionEliminationOp::create(b, t.loc, loop);
    }
  }

  // CSE only covers the block loop, since over the whole function it would merge
  // the tables back into their originals outside the block loop.
  transform::ApplyCommonSubexpressionEliminationOp::create(b, t.loc, blocks);
  transform::ApplyPatternsOp::create(b, t.loc, function,
                                     [](OpBuilder& b, Location loc)
                                     { transform::ApplyFoldTensorEmptyPatternsOp::create(b, loc); });
  // Applying patterns hoists constants, so every constant the blocks use is cloned
  // in again. Constants used only outside the blocks cannot be fused, which is
  // harmless.
  transform::SequenceOp::create(
      b, t.loc, TypeRange{}, transform::FailurePropagationMode::Suppress, blocks,
      [&](OpBuilder& b, Location loc, BlockArgument block)
      {
        Builder inner{b, loc};
        Value constants = inner.match(function, {"arith.constant"}, {});
        transform::FuseIntoContainingOp::create(b, loc, constants, block);
        transform::YieldOp::create(b, loc);
      });
  // The block-level fills are dead once fused into the stage loops.
  transform::ApplyDeadCodeEliminationOp::create(b, t.loc, function);
}

LogicalResult apply(Operation* payload, function_ref<void(Builder&, Value)> body)
{
  MLIRContext* context = payload->getContext();
  Location loc = payload->getLoc();
  OwningOpRef<ModuleOp> holder = ModuleOp::create(loc);
  OpBuilder b = OpBuilder::atBlockEnd(holder->getBody());
  auto sequence = transform::SequenceOp::create(
      b, loc, TypeRange{}, transform::FailurePropagationMode::Propagate,
      transform::AnyOpType::get(context),
      [&](OpBuilder& b, Location loc, BlockArgument root)
      {
        Builder t{b, loc};
        body(t, root);
        transform::YieldOp::create(b, loc);
      });
  return transform::applyTransforms(payload, sequence);
}

struct TileStages : PassWrapper<TileStages, OperationPass<ModuleOp>>
{
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(TileStages)

  StringRef getArgument() const final { return "ffcy-tile-stage-pencils"; }

  StringRef getDescription() const final
  {
    return "Tile each stage of interleaved cells into one loop over pencils";
  }

  void getDependentDialects(DialectRegistry& registry) const override
  {
    registry.insert<transform::TransformDialect, scf::SCFDialect>();
  }

  void runOnOperation() override
  {
    for (auto function : getOperation().getOps<func::FuncOp>())
    {
      SmallVector<Operation*> blocks;
      function.walk(
          [&](scf::ForallOp loop)
          {
            if (loop->hasAttr("ffcy.block"))
              blocks.push_back(loop);
          });
      if (blocks.size() != 1)
      {
        function.emitError("expected one block loop from the schedule");
        return signalPassFailure();
      }
      int64_t dim = 0;
      const auto stages = stages_of(blocks.front(), dim);
      if (failed(apply(function, [&](Builder& t, Value root)
                       { tile_stages(t, root, stages, dim); })))
        return signalPassFailure();
    }
  }
};
} // namespace

void register_tile_stage_pencils() { PassRegistration<TileStages>(); }
