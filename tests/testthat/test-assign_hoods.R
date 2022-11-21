# context("Testing assign_hoods")
library(miloDE)

# load data
data("sce_mouseEmbryo", package = "miloDE")
sce_mouseEmbryo = add_embedding(sce = sce_mouseEmbryo , genes = NULL, n_hvgs = 500,
                                      assay_type = "logcounts" ,
                                      reduction_type = "MNN" ,
                                      reducedDim_name = "pca_mnn",
                                      sample_id = "sample",
                                      cells_ref = colnames(sce_mouseEmbryo)[sce_mouseEmbryo$type == "wt"],
                                      cells_query = colnames(sce_mouseEmbryo)[sce_mouseEmbryo$type == "chimera"],
                                      cell_id = NULL ,
                                      d = 10)
test_that("Wrong input gives errors", {
  # k should be positive
  expect_error(assign_hoods(sce_mouseEmbryo , k = 0, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca_mnn", k_init = 50, d = 30),
               "Check k - should be positive integer",
               fixed=TRUE
  )
  # reduced dim should be in reducedDimNames
  expect_error(assign_hoods(sce = sce_mouseEmbryo , k = 10, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca_azimuth", k_init = 50, d = 30),
               "reducedDim_name should be in reducedDimNames(sce). If you do not have embedding precalculated, run first 'add_embedding'.",
               fixed=TRUE
  )
  # order can be only 1 or 2
  expect_error(assign_hoods(sce_mouseEmbryo , k = 10, prop = 0.2, order = 3, filtering = T, reducedDim_name = "pca_azimuth", k_init = 50, d = 30),
               "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)",
               fixed=TRUE
  )
  expect_error(assign_hoods(sce_mouseEmbryo , k = 10, prop = 0.2, order = "a", filtering = T, reducedDim_name = "pca_azimuth", k_init = 50, d = 30),
               "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)",
               fixed=TRUE
  )

})
