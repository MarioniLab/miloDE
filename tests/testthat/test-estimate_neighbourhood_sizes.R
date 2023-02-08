# context("Testing assign_neighbourhoods")
library(miloDE)

# load data
data("sce_mouseEmbryo", package = "miloDE")

stat = estimate_neighbourhood_sizes(x = sce_mouseEmbryo , k = c(5,10,20), prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)

# error msgs
test_that("Wrong input gives errors", {

  # x should be of the right format
  expect_error(estimate_neighbourhood_sizes(x = 1 , k_grid = c(10,30,50), prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "x should be a SingleCellExperiment or Milo object.",
               fixed=TRUE
  )
  sce_test = sce_mouseEmbryo
  colnames(sce_test) = rep(1,1,ncol(sce_test))
  expect_error(estimate_neighbourhood_sizes(sce_test , k_grid = c(10,30,50), prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "If colnames(x) exist, they should be unique.",
               fixed=TRUE
  )


  # k_grid should be of the right format
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c("10" , "30" , "50"), prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check k_grid - should be numeric vector",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(-1,0,10), prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Values of k_grid should be positive integers. Please enter valid k_grid.",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,20.5,30), prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Values of k_grid should be positive integers. Please enter valid k_grid.",
               fixed=TRUE
  )
  expect_warning(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = 20, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
                 "You only selected one value for k. If it is intended, we recommend to run directly 'assign_neighbourhoods'."
  )
  expect_warning(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(5,20,1500), prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
                 "The highest selected value is > 1000. It is gonna cost computationally, and we generally do not recommend such high k. Consider reducing."
  )


  # prop should be positive number between 0 and 1
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = "0.2", order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check prop - should be positive number between 0 and 1",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check prop - should be positive number between 0 and 1",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = -1, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check prop - should be positive number between 0 and 1",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check prop - should be positive number between 0 and 1",
               fixed=TRUE
  )


  # order should be 1 or 2
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.1, order = 0, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30), prop = 0.1, order = 3, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30), prop = 0.1, order = 1.5, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.1, order = "x", filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)",
               fixed=TRUE
  )

  # check filtering should be T or F
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.1, order = 2, filtering = 2, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check filtering - should be either TRUE or FALSE",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.1, order = 2, filtering = "aa", reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check filtering - should be either TRUE or FALSE",
               fixed=TRUE
  )

  # reducedDim should be character
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.1, order = 2, filtering = TRUE, reducedDim_name = 1, k_init = 50, d = 30),
               "Check reducedDim_name - should be character vector",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.1, order = 2, reducedDim_name = TRUE, k_init = 50, d = 30),
               "Check reducedDim_name - should be character vector",
               fixed=TRUE
  )


  # k_init - positive integer
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = "50", d = 30),
               "Check k_init - should be positive integer",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = 0, d = 30),
               "Check k_init - should be positive integer",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = -10, d = 30),
               "Check k_init - should be positive integer",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k = c(10,30,50), prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = 10.5, d = 30),
               "Check k_init - should be positive integer",
               fixed=TRUE
  )


  # d - positive integer
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = 50, d = "30"),
               "Check d - should be positive integer",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = 50, d = 0),
               "Check d - should be positive integer",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = 50, d = -1),
               "Check d - should be positive integer",
               fixed=TRUE
  )
  expect_error(estimate_neighbourhood_sizes(sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = 50, d = 1.5),
               "Check d - should be positive integer",
               fixed=TRUE
  )


  # reduced dim should be in reducedDimNames
  expect_error(estimate_neighbourhood_sizes(x = sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.2, k_init = 50, d = 30)
  )

  expect_error(estimate_neighbourhood_sizes(x = sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca_corrected", k_init = 50, d = 30),
               "reducedDim_name should be in reducedDimNames(x).",
               fixed=TRUE
  )


  # cluster_id -- NULL or in the columns of colData
  expect_error(estimate_neighbourhood_sizes(x = sce_mouseEmbryo , k_grid = c(10,30,50), prop = 0.2, order = 2, filtering = T,
                                            reducedDim_name = "pca.corrected", k_init = 50, d = 30 , cluster_id = "ct"),
               "If cluster_id not NULL, it should be in colnames(colData(x))",
               fixed=TRUE
  )

})



# return of the correct output
test_that("Return is the correct class", {
  # right class
  expect_s3_class(stat, "data.frame")

  # right colnames
  expected_colnames = c("k" , "min" , "q25" , "med" , "q75" , "max")
  expect_identical(colnames(stat) , expected_colnames)

  # right rownames
  expect_equal(nrow(stat) , 3)
})



# values increasing in both directions
test_that("Values increase as expected", {

  expect_gt(stat$min[2] , stat$min[1])
  expect_gt(stat$min[3] , stat$min[2])
  expect_gt(stat$q75[2] , stat$q25[2])
  expect_gt(stat$med[3] , stat$q25[3])
  expect_gt(stat$max[1] , stat$med[1])
  expect_gt(stat$q25[2] , stat$min[3])

})


test_that("Finishes for diff cluster_id", {
  out_null = estimate_neighbourhood_sizes(x = sce_mouseEmbryo , k_grid = c(10,30), prop = 0.2, order = 2, filtering = T,
                                         reducedDim_name = "pca.corrected", k_init = 50, d = 30 , cluster_id = NULL)
  out_ct = estimate_neighbourhood_sizes(x = sce_mouseEmbryo , k_grid = c(10,30), prop = 0.2, order = 2, filtering = T,
                                          reducedDim_name = "pca.corrected", k_init = 50, d = 30 , cluster_id = "celltype.mapped")

  expect_s3_class(out_null, "data.frame")
  expect_s3_class(out_ct, "data.frame")

})


test_that("CT grid should be reasonable", {
  expect_error(estimate_neighbourhood_sizes(x = sce_mouseEmbryo , k_grid = c(10,30,1000), prop = 0.2, order = 2, filtering = T,
                                            reducedDim_name = "pca.corrected", k_init = 50, d = 30 , cluster_id = "ct"),
               "All specified clusters have # cells < 2*max(k). We recommed to provide lower clustering resolution, decreasing max(k) or set cluster_id = NULL.",
               fixed=TRUE
  )
})



