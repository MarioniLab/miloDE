# context("Testing spatial_pval_adjustment")
library(miloDE)
library(miloR)
library(SingleCellExperiment)
library(Matrix)
library(BiocNeighbors)

# create minimal valuable data
n_row = 500
n_col = 100
n_latent = 5
set.seed(32)

sce <- SingleCellExperiment(assays = list(counts = floor(matrix(rnorm(n_row * n_col), ncol = n_col)) + 4))
rownames(sce) <- as.factor(1:n_row)
colnames(sce) <- as.character(c(1:n_col))
sce$cell <- colnames(sce)
reducedDim(sce, "reduced_dim") <- matrix(rnorm(n_col * n_latent), ncol = n_latent)



test_that("... successfully passes BNPARAM to buildGraph", {

  # 1. Run the function with default parameters (uses KmknnParam)
  set.seed(32)
  out_default <- assign_neighbourhoods(sce,
                                       reducedDim_name = "reduced_dim",
                                       verbose = FALSE,
                                       k = 10, k_init = 10) # Using smaller k for speed

  # 2. Run the function passing other options to BNPARAM via ...
  set.seed(32)
  out_annoy <- assign_neighbourhoods(sce,
                                     reducedDim_name = "reduced_dim",
                                     verbose = FALSE,
                                     k = 10, k_init = 10,
                                     BNPARAM = AnnoyParam())
  set.seed(32)
  out_hnsw <- assign_neighbourhoods(sce,
                               reducedDim_name = "reduced_dim",
                               verbose = FALSE,
                               k = 10, k_init = 10,
                               BNPARAM = HnswParam())
  set.seed(32)
  out_default_2 <- assign_neighbourhoods(sce,
                               reducedDim_name = "reduced_dim",
                               verbose = FALSE,
                               k = 10, k_init = 10,
                               BNPARAM = KmknnParam())


  # 3. Check that the outputs are different for different graph constructions
  expect_false(
    identical(miloR::graph(out_default), miloR::graph(out_annoy))
  )
  expect_false(
    identical(miloR::nhoods(out_default), miloR::nhoods(out_hnsw))
  )
  expect_true(
    identical(miloR::nhoods(out_default), miloR::nhoods(out_default_2))
  )

})

test_that("... passes invalid arguments resulting in an error", {

  # 'not_an_arg' is not a real argument for buildGraph or its dependencies
  expect_error(
    assign_neighbourhoods(sce,
                          reducedDim_name = "reduced_dim",
                          verbose = FALSE,
                          not_an_arg = "foo"),
    "unused argument" # Error message from R for unused arguments
  )

  # Test for a bad value to a real arg
  expect_error(
    assign_neighbourhoods(sce,
                          reducedDim_name = "reduced_dim",
                          verbose = FALSE,
                          BNPARAM = "not_a_param_object")
    # The specific error will come from BiocNeighbors,
    # so we just check that *an* error occurs.
  )
})
