// An mlir-opt pass plugin that moves static allocations inside GPU launches into
// workgroup memory, reusing buffers whose lifetimes do not overlap.

#include "passes.h"

#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Interfaces/ViewLikeInterface.h"
#include "mlir/Pass/Pass.h"

#include <algorithm>
#include <limits>

using namespace mlir;

namespace
{
// The top-level ops of the launch body that use `value` or a view of it.
void collect_uses(Value value, Block& body, SmallVectorImpl<Operation*>& uses)
{
  for (Operation* user : value.getUsers())
  {
    if (isa<ViewLikeOpInterface>(user))
      for (Value result : user->getResults())
        collect_uses(result, body, uses);
    uses.push_back(body.findAncestorOpInBlock(*user));
  }
}

struct Interval
{
  memref::AllocOp alloc;
  int first;
  int last;
};

// Upstream lowering only handles shared memory declared as a workgroup attribution
// of the launch, so each static allocation in a launch becomes one. The schedule
// puts a barrier after every phase, so a buffer can be reused once the last phase
// touching its previous contents has finished. Phases are counted by barriers,
// since a phase may span many ops.
struct WorkgroupAllocations
    : PassWrapper<WorkgroupAllocations, OperationPass<ModuleOp>>
{
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(WorkgroupAllocations)

  StringRef getArgument() const final { return "ffcy-workgroup-allocations"; }

  StringRef getDescription() const final
  {
    return "Move static allocations in GPU launches to reused workgroup memory";
  }

  void getDependentDialects(DialectRegistry& registry) const override
  {
    registry.insert<gpu::GPUDialect, memref::MemRefDialect>();
  }

  void runOnOperation() override
  {
    getOperation().walk([&](gpu::LaunchOp launch) { allocate(launch); });
  }

  void allocate(gpu::LaunchOp launch)
  {
    Block& body = launch.getBody().front();
    // The phase of each top-level op: the number of barriers before it.
    DenseMap<Operation*, int> phase;
    int barriers = 0;
    for (Operation& op : body)
    {
      phase[&op] = barriers;
      if (isa<gpu::BarrierOp>(op))
        ++barriers;
    }

    SmallVector<Interval> intervals;
    launch.getBody().walk(
        [&](memref::AllocOp alloc)
        {
          MemRefType type = alloc.getType();
          if (type.getMemorySpace() or !type.hasStaticShape()
              or !type.getLayout().isIdentity())
            return;
          SmallVector<Operation*> uses;
          collect_uses(alloc.getResult(), body, uses);
          if (uses.empty())
            return;
          Interval interval{alloc, std::numeric_limits<int>::max(), 0};
          for (Operation* use : uses)
          {
            interval.first = std::min(interval.first, phase[use]);
            interval.last = std::max(interval.last, phase[use]);
          }
          intervals.push_back(interval);
        });
    llvm::sort(intervals, [](const Interval& a, const Interval& b)
               { return a.first < b.first; });

    auto workgroup
        = gpu::AddressSpaceAttr::get(&getContext(), gpu::AddressSpace::Workgroup);
    struct Buffer
    {
      Value value;
      MemRefType type;
      int last;
    };
    SmallVector<Buffer> buffers;
    for (Interval& interval : intervals)
    {
      MemRefType type = interval.alloc.getType();
      auto free = llvm::find_if(buffers, [&](const Buffer& b)
                                { return b.type == type and b.last < interval.first; });
      if (free == buffers.end())
      {
        auto shared = MemRefType::get(type.getShape(), type.getElementType(),
                                      type.getLayout(), workgroup);
        buffers.push_back(
            {launch.addWorkgroupAttribution(shared, interval.alloc.getLoc()), type, 0});
        free = std::prev(buffers.end());
      }
      free->last = interval.last;
      OpBuilder builder(interval.alloc);
      auto cast = memref::MemorySpaceCastOp::create(builder, interval.alloc.getLoc(), type,
                                                    free->value);
      interval.alloc.replaceAllUsesWith(cast.getResult());
      interval.alloc.erase();
    }
  }
};
} // namespace

void register_workgroup_allocations() { PassRegistration<WorkgroupAllocations>(); }
