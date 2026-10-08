// Copyright (C) 2026 Jack S. Hale
//
// This file is part of FFCy (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    MIT

// The parts of the staged pencil schedule that depend on the form, built as
// transform dialect sequences from the IR's annotations and applied directly.
//
// ffcy-tile-stages tiles each stage from ffcy-annotate-stages onto threads. Each
// thread loops over one pencil, along the axis recorded in an op's ffcy.pencil.
// ffcy-map-threads maps the thread loops of each launch onto a block whose
// dimensions are the largest thread loops along each mapping dimension.

#include "passes.h"

#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/GPU/TransformOps/GPUTransformOps.h"
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
          // Tensors are [blocks, cells, ...] after ffcy-split-cells.
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

// Tiles each stage into one thread loop.
//
// The sinks of a stage are tiled and their loops fused, then the rest of the stage
// and its fills are fused into the loop, so a thread passes its pencils along
// without going through shared memory. Fusion drops a producer once its uses in
// the loop are fused, so producers are listed after their users, in reverse program
// order, for all their uses to be in the loop by then.
void tile_stages(Builder& t, Value function, const std::map<int64_t, Stage>& stages,
                 int64_t dim)
{
  OpBuilder& b = t.b;
  MLIRContext* context = b.getContext();
  auto attr = [&](StringRef name, int64_t value)
  { return b.getNamedAttr(name, b.getI64IntegerAttr(value)); };
  Value blocks = t.match(function, {}, {b.getNamedAttr("ffcy.block", b.getUnitAttr())});

  // The cell index takes the slowest thread dimension and the last axis the fastest.
  SmallVector<Attribute> mapping;
  for (int64_t a = dim; a-- > 0;)
    mapping.push_back(gpu::GPUThreadMappingAttr::get(context, gpu::MappingId(a)));
  ArrayAttr mapping_attr = b.getArrayAttr(mapping);

  for (const auto& [stage, info] : stages)
  {
    if (info.sinks.empty())
      continue;
    SmallVector<Value> loops;
    for (auto [axis, count] : info.sinks)
    {
      SmallVector<int64_t> sizes{0, 1};
      for (int64_t a = 2; a < dim + 2; ++a)
        sizes.push_back(a == axis ? 0 : 1);
      Value sinks = t.match(blocks, {"linalg.generic"},
                            {attr("ffcy.sink", stage), attr("ffcy.pencil", axis)});
      auto tiled = transform::TileUsingForallOp::create(
          b, t.loc, sinks, sizes, transform::TileSizesSpec(), mapping_attr);
      llvm::append_range(loops, t.split(tiled.getForallOp(), count));
    }
    // A loop is fused into a later one, which its operands dominate.
    Value loop = loops.front();
    for (Value other : ArrayRef(loops).drop_front())
      loop = transform::LoopFuseSiblingOp::create(b, t.loc, t.any, loop, other);
    transform::AnnotateOp::create(b, t.loc, loop, "ffcy.threads", Value());

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

  // Producers used more than once in a stage are fused once per use. Fused ops write
  // into slices of empty tensors, which become empty tensors of their own. CSE over
  // the whole function would merge the tables back into their originals outside the
  // block loop, so it only covers the block loop.
  transform::ApplyCommonSubexpressionEliminationOp::create(b, t.loc, blocks);
  transform::ApplyPatternsOp::create(b, t.loc, function,
                                     [](OpBuilder& b, Location loc)
                                     { transform::ApplyFoldTensorEmptyPatternsOp::create(b, loc); });
  // Applying patterns hoists constants, and folds slices of tables into new constants,
  // so every constant the blocks use is cloned in again. Constants used only outside
  // the blocks cannot be fused, which is harmless.
  transform::SequenceOp::create(
      b, t.loc, TypeRange{}, transform::FailurePropagationMode::Suppress, blocks,
      [&](OpBuilder& b, Location loc, BlockArgument block)
      {
        Builder inner{b, loc};
        Value constants = inner.match(function, {"arith.constant"}, {});
        transform::FuseIntoContainingOp::create(b, loc, constants, block);
        transform::YieldOp::create(b, loc);
      });
  // The block-level fills are dead once fused into the thread loops. DCE rather than
  // canonicalisation removes them, since canonicalisation would hoist the tables.
  transform::ApplyDeadCodeEliminationOp::create(b, t.loc, function);
}

// Builds a sequence with `body` and applies it to `payload`.
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

void transform_dialects(DialectRegistry& registry)
{
  registry.insert<transform::TransformDialect, gpu::GPUDialect, scf::SCFDialect>();
}

struct TileStages : PassWrapper<TileStages, OperationPass<ModuleOp>>
{
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(TileStages)

  StringRef getArgument() const final { return "ffcy-tile-stages"; }

  StringRef getDescription() const final
  {
    return "Tile each stage of a pencil schedule onto threads";
  }

  void getDependentDialects(DialectRegistry& registry) const override
  {
    transform_dialects(registry);
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

struct MapThreads : PassWrapper<MapThreads, OperationPass<ModuleOp>>
{
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(MapThreads)

  StringRef getArgument() const final { return "ffcy-map-threads"; }

  StringRef getDescription() const final
  {
    return "Map the thread loops of each launch onto the threads of a block";
  }

  void getDependentDialects(DialectRegistry& registry) const override
  {
    transform_dialects(registry);
  }

  void runOnOperation() override
  {
    SmallVector<gpu::LaunchOp> launches;
    getOperation().walk([&](gpu::LaunchOp launch) { launches.push_back(launch); });
    for (gpu::LaunchOp launch : launches)
    {
      SmallVector<int64_t> block_dims(3, 1);
      launch.walk(
          [&](scf::ForallOp loop)
          {
            if (!loop->hasAttr("ffcy.threads") or !loop.getMapping())
              return;
            for (auto [id, bound] : llvm::zip(*loop.getMapping(), loop.getStaticUpperBound()))
            {
              auto index = cast<gpu::GPUThreadMappingAttr>(id).getMappingId();
              block_dims[index] = std::max(block_dims[index], bound);
            }
          });
      if (failed(apply(launch,
                       [&](Builder& t, Value root)
                       {
                         transform::MapNestedForallToThreads::create(t.b, t.loc, t.any, root,
                                                                     block_dims);
                       })))
        return signalPassFailure();
    }
  }
};
} // namespace

void register_pencil_schedule()
{
  PassRegistration<TileStages>();
  PassRegistration<MapThreads>();
}
