# context("Testing assign_neighbourhoods")
library(miloDE)

# load data
data("sce_mouseEmbryo", package = "miloDE")


# error msgs
test_that("Wrong input gives errors", {

  # sce should be of the right format
  expect_error(assign_neighbourhoods(x = 1 , k = 10, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "x should be a SingleCellExperiment or Milo object.",
               fixed=TRUE
  )
  sce_test = sce_mouseEmbryo
  colnames(sce_test) = rep(1,1,ncol(sce_test))
  expect_error(assign_neighbourhoods(sce_test , k = 10, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "If colnames(x) exist, they should be unique.",
               fixed=TRUE
  )


  # k should be positive integer
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 0, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check k - should be positive integer",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = "1", prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check k - should be positive integer",
               fixed=TRUE
  )

  # prop should be positive number between 0 and 1
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = "0.2", order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check prop - should be positive number between 0 and 1",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check prop - should be positive number between 0 and 1",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = -1, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check prop - should be positive number between 0 and 1",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check prop - should be positive number between 0 and 1",
               fixed=TRUE
  )


  # order should be 1 or 2
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 0, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 3, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 1.5, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = "x", filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)",
               fixed=TRUE
  )

  # check filtering should be T or F
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 2, filtering = 2, reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check filtering - should be either TRUE or FALSE",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 2, filtering = "aa", reducedDim_name = "pca.corrected", k_init = 50, d = 30),
               "Check filtering - should be either TRUE or FALSE",
               fixed=TRUE
  )

 # reducedDim should be character
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 2, filtering = TRUE, reducedDim_name = 1, k_init = 50, d = 30),
               "Check reducedDim_name - should be character vector",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 2, reducedDim_name = TRUE, k_init = 50, d = 30),
               "Check reducedDim_name - should be character vector",
               fixed=TRUE
  )


  # k_init - positive integer
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = "50", d = 30),
               "Check k_init - should be positive integer",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = 0, d = 30),
               "Check k_init - should be positive integer",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = -10, d = 30),
               "Check k_init - should be positive integer",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = 10.5, d = 30),
               "Check k_init - should be positive integer",
               fixed=TRUE
  )


  # d - positive integer
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = 50, d = "30"),
               "Check d - should be positive integer",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = 50, d = 0),
               "Check d - should be positive integer",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = 50, d = -1),
               "Check d - should be positive integer",
               fixed=TRUE
  )
  expect_error(assign_neighbourhoods(sce_mouseEmbryo , k = 2, prop = 0.1, order = 2, reducedDim_name = "pca.corrected", k_init = 50, d = 1.5),
               "Check d - should be positive integer",
               fixed=TRUE
  )


  # reduced dim should be in reducedDimNames
  expect_error(assign_neighbourhoods(x = sce_mouseEmbryo , k = 10, prop = 0.2, k_init = 50, d = 30)
  )

  expect_error(assign_neighbourhoods(x = sce_mouseEmbryo , k = 10, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca_corrected", k_init = 50, d = 30),
               "reducedDim_name should be in reducedDimNames(x).",
               fixed=TRUE
  )

})



# return of the correct output
test_that("Return is the correct class", {
  require(miloR)
  # right class
  out = assign_neighbourhoods(x = sce_mouseEmbryo , k = 10, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  expect_s4_class(out, "Milo")

  # nhoods have data for all the cells
  nhoods_out = nhoods(out)
  expect_equal(nrow(nhoods_out) , ncol(out))
  expect_identical(sort(rownames(nhoods_out)) , sort(colnames(out)))

})


# n-hoods/hood-size ~ k, order, filtering
test_that("Scaling of neighbourhood sizes and numbers with k, order, filtering", {

  # smaller k - more neighbourhoods
  sce_1 = assign_neighbourhoods(x = sce_mouseEmbryo , k = 10, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  sce_2 = assign_neighbourhoods(x = sce_mouseEmbryo , k = 20, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  nhoods_1 = nhoods(sce_1)
  nhoods_2 = nhoods(sce_2)
  expect_gt(ncol(nhoods_1) , ncol(nhoods_2))

  # smaller k - more neighbourhoods
  sce_1 = assign_neighbourhoods(x = sce_mouseEmbryo , k = 25, prop = 0.2, order = 1, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  sce_2 = assign_neighbourhoods(x = sce_mouseEmbryo , k = 75, prop = 0.2, order = 1, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  nhoods_1 = nhoods(sce_1)
  nhoods_2 = nhoods(sce_2)
  expect_gt(ncol(nhoods_1) , ncol(nhoods_2))


  # order 1 - more neighbourhoods
  sce_1 = assign_neighbourhoods(x = sce_mouseEmbryo , k = 25, prop = 0.1, order = 1, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  sce_2 = assign_neighbourhoods(x = sce_mouseEmbryo , k = 25, prop = 0.1, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  nhoods_1 = nhoods(sce_1)
  nhoods_2 = nhoods(sce_2)
  expect_gt(ncol(nhoods_1) , ncol(nhoods_2))


  # order 1 - more neighbourhoods
  sce_1 = assign_neighbourhoods(x = sce_mouseEmbryo , k = 5, prop = 0.2, order = 1, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  sce_2 = assign_neighbourhoods(x = sce_mouseEmbryo , k = 5, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  nhoods_1 = nhoods(sce_1)
  nhoods_2 = nhoods(sce_2)
  expect_gt(ncol(nhoods_1) , ncol(nhoods_2))


  # filtering - less neighbourhoods
  sce_1 = assign_neighbourhoods(x = sce_mouseEmbryo , k = 5, prop = 0.2, order = 1, filtering = F, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  sce_2 = assign_neighbourhoods(x = sce_mouseEmbryo , k = 5, prop = 0.2, order = 1, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  nhoods_1 = nhoods(sce_1)
  nhoods_2 = nhoods(sce_2)
  expect_gt(ncol(nhoods_1) , ncol(nhoods_2))


  # filtering - less neighbourhoods
  sce_1 = assign_neighbourhoods(x = sce_mouseEmbryo , k = 20, prop = 0.1, order = 2, filtering = F, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  sce_2 = assign_neighbourhoods(x = sce_mouseEmbryo , k = 20, prop = 0.1, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  nhoods_1 = nhoods(sce_1)
  nhoods_2 = nhoods(sce_2)
  expect_gt(ncol(nhoods_1) , ncol(nhoods_2))


})





