# context("Testing assign_neighbourhoods")
library(miloDE)

# load data
data("sce_mouseEmbryo", package = "miloDE")

test_that("Wrong input gives errors", {
  # k should be positive
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 0, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca_mnn", k_init = 50, d = 30),
               "Check k - should be positive integer",
               fixed=TRUE
  )
  # reduced dim should be in reducedDimNames
  expect_error(assign_neighbourhoods(sce = sce_mouseEmbryo , k = 10, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca_azimuth", k_init = 50, d = 30),
               "reducedDim_name should be in reducedDimNames(sce).",
               fixed=TRUE
  )
  # order can be only 1 or 2
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 10, prop = 0.2, order = 3, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 10, prop = 0.2, order = "a", filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)",
               fixed=TRUE
  )

})
