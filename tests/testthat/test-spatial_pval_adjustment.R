# context("Testing spatial_pval_adjustment")
library(miloDE)
library(miloR)
library(SingleCellExperiment)

# load data
data("sce_mouseEmbryo", package = "miloDE")
set.seed(32)
sce = assign_neighbourhoods(sce_mouseEmbryo , k = 25, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
nhoods_sce = nhoods(sce)
pvals_1 = rep(1,0.0001,ncol(nhoods_sce))
pvals_2 = rep(1,0.001,ncol(nhoods_sce))
pvals_3 = rep(1,0.01,ncol(nhoods_sce))
pvals_4 = rep(1,0.1,ncol(nhoods_sce))
pvals_5 = pvals_4
pvals_5[c(4,1,8)] = NaN


test_that("Right length", {
  out = spatial_pval_adjustment(nhoods_sce , pvals_2)
  expect_equal(ncol(nhoods_sce) , length(out))
})

test_that("NaNs at right places", {

  out = spatial_pval_adjustment(nhoods_sce , pvals_5)
  expect_identical(out[1] , NaN)
  expect_identical(out[4] , NaN)
  expect_identical(out[8] , NaN)
})

