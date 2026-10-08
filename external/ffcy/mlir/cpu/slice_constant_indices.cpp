// Reads operands of linalg.generic at constant indices through slices instead.
//
// An indexing map with a constant result, such as one field of the geometry data,
// is not a projected permutation, so the op cannot be vectorised. Each such input
// becomes a rank-reduced slice at those indices, and its map drops them.

#include "passes.h"

#include "mlir/Dialect/Linalg/IR/Linalg.h"
#include "mlir/Dialect/Tensor/IR/Tensor.h"
#include "mlir/IR/PatternMatch.h"
#include "mlir/Pass/Pass.h"
#include "mlir/Transforms/GreedyPatternRewriteDriver.h"

using namespace mlir;

namespace
{
struct SliceConstantIndicesPattern : OpRewritePattern<linalg::GenericOp>
{
  using OpRewritePattern<linalg::GenericOp>::OpRewritePattern;

  LogicalResult matchAndRewrite(linalg::GenericOp op, PatternRewriter& rewriter) const override
  {
    SmallVector<AffineMap> maps = op.getIndexingMapsArray();
    bool changed = false;
    for (OpOperand* operand : op.getDpsInputOperands())
    {
      Value value = operand->get();
      auto type = dyn_cast<RankedTensorType>(value.getType());
      AffineMap map = maps[operand->getOperandNumber()];
      if (!type or llvm::none_of(map.getResults(), llvm::IsaPred<AffineConstantExpr>))
        continue;

      rewriter.setInsertionPoint(op);
      SmallVector<OpFoldResult> offsets, sizes;
      SmallVector<OpFoldResult> strides(type.getRank(), rewriter.getIndexAttr(1));
      SmallVector<AffineExpr> kept;
      SmallVector<int64_t> shape;
      for (auto [axis, expr] : llvm::enumerate(map.getResults()))
      {
        if (auto constant = dyn_cast<AffineConstantExpr>(expr))
        {
          offsets.push_back(rewriter.getIndexAttr(constant.getValue()));
          sizes.push_back(rewriter.getIndexAttr(1));
          continue;
        }
        offsets.push_back(rewriter.getIndexAttr(0));
        sizes.push_back(tensor::getMixedSize(rewriter, op.getLoc(), value, axis));
        kept.push_back(expr);
        shape.push_back(type.getDimSize(axis));
      }
      Value slice = tensor::ExtractSliceOp::create(
          rewriter, op.getLoc(), RankedTensorType::get(shape, type.getElementType()), value,
          offsets, sizes, strides);
      rewriter.modifyOpInPlace(op, [&] { operand->set(slice); });
      maps[operand->getOperandNumber()]
          = AffineMap::get(map.getNumDims(), map.getNumSymbols(), kept, op.getContext());
      changed = true;
    }
    if (!changed)
      return failure();
    rewriter.modifyOpInPlace(op, [&]
                             { op.setIndexingMapsAttr(rewriter.getAffineMapArrayAttr(maps)); });
    return success();
  }
};

struct SliceConstantIndices : PassWrapper<SliceConstantIndices, OperationPass<>>
{
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(SliceConstantIndices)

  StringRef getArgument() const final { return "ffcy-slice-constant-indices"; }

  StringRef getDescription() const final
  {
    return "Read linalg.generic inputs at constant indices through slices";
  }

  void getDependentDialects(DialectRegistry& registry) const override
  {
    registry.insert<linalg::LinalgDialect, tensor::TensorDialect>();
  }

  void runOnOperation() override
  {
    RewritePatternSet patterns(&getContext());
    patterns.add<SliceConstantIndicesPattern>(&getContext());
    if (failed(applyPatternsGreedily(getOperation(), std::move(patterns))))
      signalPassFailure();
  }
};
} // namespace

void register_slice_constant_indices() { PassRegistration<SliceConstantIndices>(); }
