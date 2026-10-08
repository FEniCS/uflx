// Sorts the ops of each block by the stages that ffcy-annotate-stages assigned, so
// that the ops of a stage are contiguous and their loops can be fused as siblings.
//
// An op tagged ffcy.stage or ffcy.sink goes in its stage. Any other op goes with the
// latest stage it depends on, after that stage's sinks if it reads one, since only
// later stages read sinks. Ties keep program order, so the order stays topological.

#include "passes.h"

#include "mlir/IR/BuiltinAttributes.h"
#include "mlir/Pass/Pass.h"
#include "mlir/Transforms/RegionUtils.h"

#include <algorithm>
#include <tuple>
#include <vector>

using namespace mlir;

namespace
{
struct Key
{
  std::int64_t stage = 0;
  int after_sinks = 0;
  bool operator<(const Key& other) const
  {
    return std::tie(stage, after_sinks) < std::tie(other.stage, other.after_sinks);
  }
};

void sort_stages(Block& block)
{
  std::vector<Operation*> ops;
  for (Operation& op : block.without_terminator())
    ops.push_back(&op);

  llvm::DenseMap<Operation*, Key> keys;
  for (Operation* op : ops)
  {
    if (auto sink = op->getAttrOfType<IntegerAttr>("ffcy.sink"))
      keys[op] = {sink.getInt(), 0};
    else if (auto stage = op->getAttrOfType<IntegerAttr>("ffcy.stage"))
      keys[op] = {stage.getInt(), 0};
    else
    {
      SetVector<Value> used(op->getOperands().begin(), op->getOperands().end());
      getUsedValuesDefinedAbove(op->getRegions(), used);
      Key key;
      for (Value value : used)
      {
        Operation* def = value.getDefiningOp();
        if (!def or def->getBlock() != &block)
          continue;
        Key dependency = keys.at(def);
        if (def->hasAttr("ffcy.sink"))
          dependency.after_sinks = 1;
        key = std::max(key, dependency);
      }
      keys[op] = key;
    }
  }

  std::stable_sort(ops.begin(), ops.end(),
                   [&](Operation* a, Operation* b) { return keys.at(a) < keys.at(b); });
  for (Operation* op : ops)
    op->moveBefore(block.getTerminator());
}

struct SortStages : PassWrapper<SortStages, OperationPass<>>
{
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(SortStages)

  StringRef getArgument() const final { return "ffcy-sort-stages"; }

  StringRef getDescription() const final
  {
    return "Sort ops by the stages of ffcy-annotate-stages, keeping a topological order";
  }

  void runOnOperation() override
  {
    SetVector<Block*> blocks;
    getOperation()->walk(
        [&](Operation* op)
        {
          if (op->hasAttr("ffcy.sink"))
            blocks.insert(op->getBlock());
        });
    for (Block* block : blocks)
      sort_stages(*block);
  }
};
} // namespace

void register_sort_stages() { PassRegistration<SortStages>(); }
