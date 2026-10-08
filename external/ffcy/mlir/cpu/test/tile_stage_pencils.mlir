// Tiles each stage of a block of interleaved cells into one loop over pencils.
// The input is Laplace on quadrilaterals at Q1 with inline geometry, after
// @__blocks and ffcy-sort-stages. Its second stage computes the inverse Jacobian,
// whose four divisions share producers, so each must be fused once.

// CHECK-LABEL: func.func @action
// CHECK:       scf.forall
// CHECK-NOT:   linalg.generic
// CHECK:       scf.forall
// CHECK:       } {ffcy.pencils}
// CHECK-COUNT-4: arith.divf
// CHECK-NOT:   arith.divf
// CHECK:       } {ffcy.pencils}
// CHECK:       } {ffcy.pencils}
// CHECK:       } {ffcy.pencils}
// CHECK-NOT:   linalg.generic
// CHECK:       } {ffcy.block}

#map = affine_map<(d0, d1, d2, d3, d4) -> (d1, d3)>
#map1 = affine_map<(d0, d1, d2, d3, d4) -> (d0, d3, d2, d4)>
#map2 = affine_map<(d0, d1, d2, d3, d4) -> (d0, d1, d2, d4)>
#map3 = affine_map<(d0, d1, d2, d3) -> (d0, 0, 0, 1, d3)>
#map4 = affine_map<(d0, d1, d2, d3) -> (d0, 0, 1, 1, d3)>
#map5 = affine_map<(d0, d1, d2, d3) -> (d0, 1, 0, 1, d3)>
#map6 = affine_map<(d0, d1, d2, d3) -> (d0, 1, 1, 1, d3)>
#map7 = affine_map<(d0, d1, d2, d3) -> (d1, 0)>
#map8 = affine_map<(d0, d1, d2, d3) -> (d1, 1)>
#map9 = affine_map<(d0, d1, d2, d3) -> (d2, 0)>
#map10 = affine_map<(d0, d1, d2, d3) -> (d2, 1)>
#map11 = affine_map<(d0, d1, d2, d3) -> (d0, d1, d2, d3)>
#map12 = affine_map<(d0, d1, d2, d3) -> (d0, 0, 0, 0, d3)>
#map13 = affine_map<(d0, d1, d2, d3) -> (d0, 0, 1, 0, d3)>
#map14 = affine_map<(d0, d1, d2, d3) -> (d0, 1, 0, 0, d3)>
#map15 = affine_map<(d0, d1, d2, d3) -> (d0, 1, 1, 0, d3)>
#map16 = affine_map<(d0, d1, d2, d3, d4) -> (d2, d3)>
#map17 = affine_map<(d0, d1, d2, d3, d4) -> (d0, d1, d3, d4)>
#map18 = affine_map<(d0, d1, d2, d3) -> (d1, d2)>
module {
  func.func @action(%arg0: tensor<?x4x2x4xf64>, %arg1: tensor<?x4x4xf64>) -> tensor<?x2x2x4xf64> attributes {llvm.emit_c_interface} {
    %cst = arith.constant {ffcy.table} dense<[[-0.99999999999999988, -0.99999999999999988], [1.000000e+00, 1.000000e+00]]> : tensor<2x2xf64>
    %cst_0 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518719], [0.21132486540518713, 0.78867513459481287]]> : tensor<2x2xf64>
    %cst_1 = arith.constant {ffcy.table} dense<2.500000e-01> : tensor<2x2xf64>
    %cst_2 = arith.constant -0.99999999999999988 : f64
    %cst_3 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
    %cst_4 = arith.constant 0.000000e+00 : f64
    %cst_5 = arith.constant {ffcy.table} dense<[[-0.99999999999999988, 1.000000e+00], [-0.99999999999999988, 1.000000e+00]]> : tensor<2x2xf64>
    %c0 = arith.constant 0 : index
    %dim = tensor.dim %arg0, %c0 : tensor<?x4x2x4xf64>
    %0 = tensor.empty(%dim) {ffcy.buffer = 42 : i64} : tensor<?x2x2x4xf64>
    %1 = scf.forall (%arg2) in (%dim) shared_outs(%arg3 = %0) -> (tensor<?x2x2x4xf64>) {
      %extracted_slice = tensor.extract_slice %arg0[%arg2, 0, 0, 0] [1, 4, 2, 4] [1, 1, 1, 1] : tensor<?x4x2x4xf64> to tensor<1x4x2x4xf64>
      %expanded = tensor.expand_shape %extracted_slice [[0], [1, 2], [3], [4]] output_shape [1, 2, 2, 2, 4] : tensor<1x4x2x4xf64> into tensor<1x2x2x2x4xf64>
      %2 = tensor.empty() : tensor<1x2x2x4xf64>
      %cst_6 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %cst_7 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %cst_8 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %cst_9 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %3 = tensor.empty() : tensor<1x2x2x4xf64>
      %cst_10 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %cst_11 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %cst_12 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %cst_13 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %4 = tensor.empty() : tensor<1x2x2x4xf64>
      %cst_14 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %cst_15 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %cst_16 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %cst_17 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %5 = tensor.empty() : tensor<1x2x2x4xf64>
      %cst_18 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %cst_19 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %cst_20 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %cst_21 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %6 = tensor.empty() : tensor<1x2x2x4xf64>
      %7 = tensor.empty() : tensor<1x2x2x4xf64>
      %8 = tensor.empty() : tensor<1x2x2x4xf64>
      %9 = tensor.empty() : tensor<1x2x2x4xf64>
      %extracted_slice_22 = tensor.extract_slice %arg1[%arg2, 0, 0] [1, 4, 4] [1, 1, 1] : tensor<?x4x4xf64> to tensor<1x4x4xf64>
      %expanded_23 = tensor.expand_shape %extracted_slice_22 [[0], [1, 2], [3]] output_shape [1, 2, 2, 4] : tensor<1x4x4xf64> into tensor<1x2x2x4xf64>
      %10 = tensor.empty() : tensor<1x2x2x4xf64>
      %11 = linalg.fill {ffcy.stage = 0 : i64} ins(%cst_4 : f64) outs(%10 : tensor<1x2x2x4xf64>) -> tensor<1x2x2x4xf64>
      %cst_24 = arith.constant {ffcy.table} dense<[[-0.99999999999999988, 1.000000e+00], [-0.99999999999999988, 1.000000e+00]]> : tensor<2x2xf64>
      %12 = linalg.generic {indexing_maps = [#map, #map1, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction", "parallel"]} ins(%cst_24, %expanded_23 : tensor<2x2xf64>, tensor<1x2x2x4xf64>) outs(%11 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 1 : i64, ffcy.sink = 0 : i64} {
      ^bb0(%in: f64, %in_35: f64, %out: f64):
        %53 = arith.mulf %in, %in_35 fastmath<contract> : f64
        %54 = arith.addf %out, %53 fastmath<contract> : f64
        linalg.yield %54 : f64
      } -> tensor<1x2x2x4xf64>
      %13 = tensor.empty() : tensor<1x2x2x4xf64>
      %cst_25 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %14 = tensor.empty() : tensor<1x2x2x4xf64>
      %15 = tensor.empty() : tensor<1x2x2x4xf64>
      %16 = tensor.empty() : tensor<1x2x2x4xf64>
      %17 = tensor.empty() : tensor<1x2x2x4xf64>
      %18 = linalg.fill {ffcy.stage = 0 : i64} ins(%cst_4 : f64) outs(%17 : tensor<1x2x2x4xf64>) -> tensor<1x2x2x4xf64>
      %cst_26 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518713], [0.21132486540518719, 0.78867513459481287]]> : tensor<2x2xf64>
      %19 = linalg.generic {indexing_maps = [#map, #map1, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction", "parallel"]} ins(%cst_26, %expanded_23 : tensor<2x2xf64>, tensor<1x2x2x4xf64>) outs(%18 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 1 : i64, ffcy.sink = 0 : i64} {
      ^bb0(%in: f64, %in_35: f64, %out: f64):
        %53 = arith.mulf %in, %in_35 fastmath<contract> : f64
        %54 = arith.addf %out, %53 fastmath<contract> : f64
        linalg.yield %54 : f64
      } -> tensor<1x2x2x4xf64>
      %20 = tensor.empty() : tensor<1x2x2x4xf64>
      %cst_27 = arith.constant {ffcy.table} dense<[[-0.99999999999999988, 1.000000e+00], [-0.99999999999999988, 1.000000e+00]]> : tensor<2x2xf64>
      %21 = tensor.empty() : tensor<1x2x2x4xf64>
      %cst_28 = arith.constant {ffcy.table} dense<2.500000e-01> : tensor<2x2xf64>
      %22 = tensor.empty() : tensor<1x2x2x4xf64>
      %cst_29 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518719], [0.21132486540518713, 0.78867513459481287]]> : tensor<2x2xf64>
      %23 = tensor.empty() : tensor<1x2x2x4xf64>
      %cst_30 = arith.constant {ffcy.table} dense<[[-0.99999999999999988, -0.99999999999999988], [1.000000e+00, 1.000000e+00]]> : tensor<2x2xf64>
      %24 = tensor.empty() : tensor<1x2x2x4xf64>
      %cst_31 = arith.constant {ffcy.table} dense<2.500000e-01> : tensor<2x2xf64>
      %25 = tensor.empty() : tensor<1x2x2x4xf64>
      %cst_32 = arith.constant {ffcy.table} dense<[[-0.99999999999999988, -0.99999999999999988], [1.000000e+00, 1.000000e+00]]> : tensor<2x2xf64>
      %26 = tensor.empty() : tensor<1x2x2x4xf64>
      %cst_33 = arith.constant {ffcy.table} dense<[[0.78867513459481287, 0.21132486540518719], [0.21132486540518713, 0.78867513459481287]]> : tensor<2x2xf64>
      %extracted_slice_34 = tensor.extract_slice %arg3[%arg2, 0, 0, 0] [1, 2, 2, 4] [1, 1, 1, 1] : tensor<?x2x2x4xf64> to tensor<1x2x2x4xf64>
      %27 = linalg.generic {indexing_maps = [#map3, #map4, #map5, #map6, #map7, #map8, #map9, #map10, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%expanded, %expanded, %expanded, %expanded, %cst_6, %cst_7, %cst_8, %cst_9 : tensor<1x2x2x2x4xf64>, tensor<1x2x2x2x4xf64>, tensor<1x2x2x2x4xf64>, tensor<1x2x2x2x4xf64>, tensor<2x2xf64>, tensor<2x2xf64>, tensor<2x2xf64>, tensor<2x2xf64>) outs(%2 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.stage = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %in_36: f64, %in_37: f64, %in_38: f64, %in_39: f64, %in_40: f64, %in_41: f64, %out: f64):
        %53 = arith.mulf %in, %cst_2 fastmath<contract> : f64
        %54 = arith.addf %53, %in_35 fastmath<contract> : f64
        %55 = arith.mulf %in_38, %54 fastmath<contract> : f64
        %56 = arith.mulf %in_36, %cst_2 fastmath<contract> : f64
        %57 = arith.addf %56, %in_37 fastmath<contract> : f64
        %58 = arith.mulf %in_39, %57 fastmath<contract> : f64
        %59 = arith.addf %55, %58 fastmath<contract> : f64
        linalg.yield %59 : f64
      } -> tensor<1x2x2x4xf64>
      %28 = linalg.generic {indexing_maps = [#map12, #map13, #map14, #map15, #map7, #map8, #map9, #map10, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%expanded, %expanded, %expanded, %expanded, %cst_10, %cst_11, %cst_12, %cst_13 : tensor<1x2x2x2x4xf64>, tensor<1x2x2x2x4xf64>, tensor<1x2x2x2x4xf64>, tensor<1x2x2x2x4xf64>, tensor<2x2xf64>, tensor<2x2xf64>, tensor<2x2xf64>, tensor<2x2xf64>) outs(%3 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.stage = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %in_36: f64, %in_37: f64, %in_38: f64, %in_39: f64, %in_40: f64, %in_41: f64, %out: f64):
        %53 = arith.mulf %in, %cst_2 fastmath<contract> : f64
        %54 = arith.addf %53, %in_36 fastmath<contract> : f64
        %55 = arith.mulf %in_40, %54 fastmath<contract> : f64
        %56 = arith.mulf %in_35, %cst_2 fastmath<contract> : f64
        %57 = arith.addf %56, %in_37 fastmath<contract> : f64
        %58 = arith.mulf %in_41, %57 fastmath<contract> : f64
        %59 = arith.addf %55, %58 fastmath<contract> : f64
        linalg.yield %59 : f64
      } -> tensor<1x2x2x4xf64>
      %29 = linalg.generic {indexing_maps = [#map3, #map4, #map5, #map6, #map7, #map8, #map9, #map10, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%expanded, %expanded, %expanded, %expanded, %cst_14, %cst_15, %cst_16, %cst_17 : tensor<1x2x2x2x4xf64>, tensor<1x2x2x2x4xf64>, tensor<1x2x2x2x4xf64>, tensor<1x2x2x2x4xf64>, tensor<2x2xf64>, tensor<2x2xf64>, tensor<2x2xf64>, tensor<2x2xf64>) outs(%4 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.stage = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %in_36: f64, %in_37: f64, %in_38: f64, %in_39: f64, %in_40: f64, %in_41: f64, %out: f64):
        %53 = arith.mulf %in, %cst_2 fastmath<contract> : f64
        %54 = arith.addf %53, %in_36 fastmath<contract> : f64
        %55 = arith.mulf %in_40, %54 fastmath<contract> : f64
        %56 = arith.mulf %in_35, %cst_2 fastmath<contract> : f64
        %57 = arith.addf %56, %in_37 fastmath<contract> : f64
        %58 = arith.mulf %in_41, %57 fastmath<contract> : f64
        %59 = arith.addf %55, %58 fastmath<contract> : f64
        linalg.yield %59 : f64
      } -> tensor<1x2x2x4xf64>
      %30 = linalg.generic {indexing_maps = [#map12, #map13, #map14, #map15, #map7, #map8, #map9, #map10, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%expanded, %expanded, %expanded, %expanded, %cst_18, %cst_19, %cst_20, %cst_21 : tensor<1x2x2x2x4xf64>, tensor<1x2x2x2x4xf64>, tensor<1x2x2x2x4xf64>, tensor<1x2x2x2x4xf64>, tensor<2x2xf64>, tensor<2x2xf64>, tensor<2x2xf64>, tensor<2x2xf64>) outs(%5 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.stage = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %in_36: f64, %in_37: f64, %in_38: f64, %in_39: f64, %in_40: f64, %in_41: f64, %out: f64):
        %53 = arith.mulf %in, %cst_2 fastmath<contract> : f64
        %54 = arith.addf %53, %in_35 fastmath<contract> : f64
        %55 = arith.mulf %in_38, %54 fastmath<contract> : f64
        %56 = arith.mulf %in_36, %cst_2 fastmath<contract> : f64
        %57 = arith.addf %56, %in_37 fastmath<contract> : f64
        %58 = arith.mulf %in_39, %57 fastmath<contract> : f64
        %59 = arith.addf %55, %58 fastmath<contract> : f64
        linalg.yield %59 : f64
      } -> tensor<1x2x2x4xf64>
      %31 = linalg.generic {indexing_maps = [#map11, #map11, #map11, #map11, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%28, %27, %29, %30 : tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>) outs(%6 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.stage = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %in_36: f64, %in_37: f64, %out: f64):
        %53 = arith.mulf %in_36, %in_37 fastmath<contract> : f64
        %54 = arith.mulf %in, %in_35 fastmath<contract> : f64
        %55 = arith.subf %54, %53 fastmath<contract> : f64
        linalg.yield %55 : f64
      } -> tensor<1x2x2x4xf64>
      %32 = linalg.generic {indexing_maps = [#map11, #map11, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%27, %31 : tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>) outs(%7 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.stage = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %out: f64):
        %53 = arith.divf %in, %in_35 fastmath<contract> : f64
        linalg.yield %53 : f64
      } -> tensor<1x2x2x4xf64>
      %33 = linalg.generic {indexing_maps = [#map11, #map11, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%30, %31 : tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>) outs(%8 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.stage = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %out: f64):
        %53 = arith.negf %in : f64
        %54 = arith.divf %53, %in_35 fastmath<contract> : f64
        linalg.yield %54 : f64
      } -> tensor<1x2x2x4xf64>
      %34 = linalg.generic {indexing_maps = [#map11, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%31 : tensor<1x2x2x4xf64>) outs(%9 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.stage = 1 : i64} {
      ^bb0(%in: f64, %out: f64):
        %53 = math.absf %in : f64
        linalg.yield %53 : f64
      } -> tensor<1x2x2x4xf64>
      %35 = linalg.fill {ffcy.stage = 1 : i64} ins(%cst_4 : f64) outs(%13 : tensor<1x2x2x4xf64>) -> tensor<1x2x2x4xf64>
      %36 = linalg.generic {indexing_maps = [#map16, #map17, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction", "parallel"]} ins(%cst_25, %12 : tensor<2x2xf64>, tensor<1x2x2x4xf64>) outs(%35 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.stage = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %out: f64):
        %53 = arith.mulf %in, %in_35 fastmath<contract> : f64
        %54 = arith.addf %out, %53 fastmath<contract> : f64
        linalg.yield %54 : f64
      } -> tensor<1x2x2x4xf64>
      %37 = linalg.generic {indexing_maps = [#map11, #map11, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%29, %31 : tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>) outs(%14 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.stage = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %out: f64):
        %53 = arith.negf %in : f64
        %54 = arith.divf %53, %in_35 fastmath<contract> : f64
        linalg.yield %54 : f64
      } -> tensor<1x2x2x4xf64>
      %38 = linalg.generic {indexing_maps = [#map11, #map11, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%28, %31 : tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>) outs(%15 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.stage = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %out: f64):
        %53 = arith.divf %in, %in_35 fastmath<contract> : f64
        linalg.yield %53 : f64
      } -> tensor<1x2x2x4xf64>
      %39 = linalg.generic {indexing_maps = [#map11, #map11, #map11, #map11, #map11, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%32, %37, %33, %38, %34 : tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>) outs(%16 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.stage = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %in_36: f64, %in_37: f64, %in_38: f64, %out: f64):
        %53 = arith.mulf %in_36, %in_37 fastmath<contract> : f64
        %54 = arith.mulf %in, %in_35 fastmath<contract> : f64
        %55 = arith.addf %54, %53 fastmath<contract> : f64
        %56 = arith.mulf %55, %in_38 fastmath<contract> : f64
        linalg.yield %56 : f64
      } -> tensor<1x2x2x4xf64>
      %40 = linalg.fill {ffcy.stage = 1 : i64} ins(%cst_4 : f64) outs(%20 : tensor<1x2x2x4xf64>) -> tensor<1x2x2x4xf64>
      %41 = linalg.generic {indexing_maps = [#map16, #map17, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction", "parallel"]} ins(%cst_27, %19 : tensor<2x2xf64>, tensor<1x2x2x4xf64>) outs(%40 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.stage = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %out: f64):
        %53 = arith.mulf %in, %in_35 fastmath<contract> : f64
        %54 = arith.addf %out, %53 fastmath<contract> : f64
        linalg.yield %54 : f64
      } -> tensor<1x2x2x4xf64>
      %42 = linalg.generic {indexing_maps = [#map11, #map11, #map11, #map11, #map11, #map11, #map11, #map11, #map18, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%32, %32, %33, %33, %34, %36, %39, %41, %cst_28 : tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<2x2xf64>) outs(%21 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.sink = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %in_36: f64, %in_37: f64, %in_38: f64, %in_39: f64, %in_40: f64, %in_41: f64, %in_42: f64, %out: f64):
        %53 = arith.mulf %in_36, %in_37 fastmath<contract> : f64
        %54 = arith.mulf %in, %in_35 fastmath<contract> : f64
        %55 = arith.addf %54, %53 fastmath<contract> : f64
        %56 = arith.mulf %in_40, %in_41 fastmath<contract> : f64
        %57 = arith.mulf %55, %in_38 fastmath<contract> : f64
        %58 = arith.mulf %57, %in_39 fastmath<contract> : f64
        %59 = arith.addf %58, %56 fastmath<contract> : f64
        %60 = arith.mulf %59, %in_42 fastmath<contract> : f64
        linalg.yield %60 : f64
      } -> tensor<1x2x2x4xf64>
      %43 = linalg.generic {indexing_maps = [#map11, #map11, #map11, #map11, #map11, #map11, #map11, #map11, #map18, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%39, %36, %37, %37, %38, %38, %34, %41, %cst_31 : tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>, tensor<2x2xf64>) outs(%24 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.sink = 1 : i64} {
      ^bb0(%in: f64, %in_35: f64, %in_36: f64, %in_37: f64, %in_38: f64, %in_39: f64, %in_40: f64, %in_41: f64, %in_42: f64, %out: f64):
        %53 = arith.mulf %in_38, %in_39 fastmath<contract> : f64
        %54 = arith.mulf %in_36, %in_37 fastmath<contract> : f64
        %55 = arith.addf %54, %53 fastmath<contract> : f64
        %56 = arith.mulf %55, %in_40 fastmath<contract> : f64
        %57 = arith.mulf %56, %in_41 fastmath<contract> : f64
        %58 = arith.mulf %in, %in_35 fastmath<contract> : f64
        %59 = arith.addf %58, %57 fastmath<contract> : f64
        %60 = arith.mulf %59, %in_42 fastmath<contract> : f64
        linalg.yield %60 : f64
      } -> tensor<1x2x2x4xf64>
      %44 = linalg.fill {ffcy.stage = 2 : i64} ins(%cst_4 : f64) outs(%22 : tensor<1x2x2x4xf64>) -> tensor<1x2x2x4xf64>
      %45 = linalg.generic {indexing_maps = [#map16, #map17, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction", "parallel"]} ins(%cst_29, %42 : tensor<2x2xf64>, tensor<1x2x2x4xf64>) outs(%44 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.sink = 2 : i64} {
      ^bb0(%in: f64, %in_35: f64, %out: f64):
        %53 = arith.mulf %in, %in_35 fastmath<contract> : f64
        %54 = arith.addf %out, %53 fastmath<contract> : f64
        linalg.yield %54 : f64
      } -> tensor<1x2x2x4xf64>
      %46 = linalg.fill {ffcy.stage = 2 : i64} ins(%cst_4 : f64) outs(%25 : tensor<1x2x2x4xf64>) -> tensor<1x2x2x4xf64>
      %47 = linalg.generic {indexing_maps = [#map16, #map17, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction", "parallel"]} ins(%cst_32, %43 : tensor<2x2xf64>, tensor<1x2x2x4xf64>) outs(%46 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 2 : i64, ffcy.sink = 2 : i64} {
      ^bb0(%in: f64, %in_35: f64, %out: f64):
        %53 = arith.mulf %in, %in_35 fastmath<contract> : f64
        %54 = arith.addf %out, %53 fastmath<contract> : f64
        linalg.yield %54 : f64
      } -> tensor<1x2x2x4xf64>
      %48 = linalg.fill {ffcy.stage = 3 : i64} ins(%cst_4 : f64) outs(%23 : tensor<1x2x2x4xf64>) -> tensor<1x2x2x4xf64>
      %49 = linalg.generic {indexing_maps = [#map, #map1, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction", "parallel"]} ins(%cst_30, %45 : tensor<2x2xf64>, tensor<1x2x2x4xf64>) outs(%48 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 1 : i64, ffcy.stage = 3 : i64} {
      ^bb0(%in: f64, %in_35: f64, %out: f64):
        %53 = arith.mulf %in, %in_35 fastmath<contract> : f64
        %54 = arith.addf %out, %53 fastmath<contract> : f64
        linalg.yield %54 : f64
      } -> tensor<1x2x2x4xf64>
      %50 = linalg.fill {ffcy.stage = 3 : i64} ins(%cst_4 : f64) outs(%26 : tensor<1x2x2x4xf64>) -> tensor<1x2x2x4xf64>
      %51 = linalg.generic {indexing_maps = [#map, #map1, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction", "parallel"]} ins(%cst_33, %47 : tensor<2x2xf64>, tensor<1x2x2x4xf64>) outs(%50 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 1 : i64, ffcy.stage = 3 : i64} {
      ^bb0(%in: f64, %in_35: f64, %out: f64):
        %53 = arith.mulf %in, %in_35 fastmath<contract> : f64
        %54 = arith.addf %out, %53 fastmath<contract> : f64
        linalg.yield %54 : f64
      } -> tensor<1x2x2x4xf64>
      %52 = linalg.generic {indexing_maps = [#map11, #map11, #map11], iterator_types = ["parallel", "parallel", "parallel", "parallel"]} ins(%49, %51 : tensor<1x2x2x4xf64>, tensor<1x2x2x4xf64>) outs(%extracted_slice_34 : tensor<1x2x2x4xf64>) attrs =  {ffcy.pencil = 1 : i64, ffcy.sink = 3 : i64} {
      ^bb0(%in: f64, %in_35: f64, %out: f64):
        %53 = arith.addf %in, %in_35 fastmath<contract> : f64
        linalg.yield %53 : f64
      } -> tensor<1x2x2x4xf64>
      scf.forall.in_parallel {
        tensor.parallel_insert_slice %52 into %arg3[%arg2, 0, 0, 0] [1, 2, 2, 4] [1, 1, 1, 1] : tensor<1x2x2x4xf64> into tensor<?x2x2x4xf64>
      }
    } {ffcy.block}
    return %1 : tensor<?x2x2x4xf64>
  }
}

