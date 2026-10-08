// Upstream elementwise fusion, limited to pointwise consumers and producers with
// one use. Fusing into a contraction would recompute the producer for each
// reduction iteration, and a producer with other uses would be computed twice.
// The upstream pass also folds constants and deduplicates operands, which fuses
// inline geometry into one large op that runs slower.

#include "passes.h"

#include "mlir/Dialect/Linalg/IR/Linalg.h"
#include "mlir/Dialect/Linalg/Transforms/Transforms.h"
#include "mlir/IR/PatternMatch.h"
#include "mlir/Pass/Pass.h"
#include "mlir/Transforms/GreedyPatternRewriteDriver.h"

using namespace mlir;

namespace
{
struct FuseElementwiseOps : OpRewritePattern<linalg::GenericOp>
{
  using OpRewritePattern<linalg::GenericOp>::OpRewritePattern;

  LogicalResult matchAndRewrite(linalg::GenericOp consumer,
                                PatternRewriter& rewriter) const override
  {
    if (consumer.getNumReductionLoops() > 0)
      return failure();
    for (OpOperand& operand : consumer->getOpOperands())
    {
      Operation* producer = operand.get().getDefiningOp();
      if (!linalg::areElementwiseOpsFusable(&operand) or !producer->hasOneUse())
        continue;
      FailureOr<linalg::ElementwiseOpFusionResult> fused
          = linalg::fuseElementwiseOps(rewriter, &operand);
      if (failed(fused))
        return failure();
      for (auto [original, replacement] : fused->replacements)
        rewriter.replaceUsesWithIf(original, replacement, [&](OpOperand& use)
                                   { return use.get().getDefiningOp() != producer; });
      rewriter.eraseOp(consumer);
      return success();
    }
    return failure();
  }
};

struct FuseElementwise : PassWrapper<FuseElementwise, OperationPass<>>
{
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(FuseElementwise)

  StringRef getArgument() const final { return "ffcy-fuse-elementwise"; }

  StringRef getDescription() const final
  {
    return "Fuse elementwise linalg ops into pointwise consumers";
  }

  void getDependentDialects(DialectRegistry& registry) const override
  {
    registry.insert<linalg::LinalgDialect>();
  }

  void runOnOperation() override
  {
    RewritePatternSet patterns(&getContext());
    patterns.add<FuseElementwiseOps>(&getContext());
    if (failed(applyPatternsGreedily(getOperation(), std::move(patterns))))
      signalPassFailure();
  }
};
} // namespace

void register_fuse_elementwise() { PassRegistration<FuseElementwise>(); }
