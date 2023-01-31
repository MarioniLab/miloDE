# context("Testing filter_neighbourhoods")
library(miloDE)
library(miloR)

# load data
data("sce_mouseEmbryo", package = "miloDE")


# error msgs
test_that("Wrong input gives errors", {
  # sce should be of the right format
  expect_error(filter_neighbourhoods(x = 1),
               "x should be a Milo object. Please run `assign_neighbourhoods` first.",
               fixed=TRUE
  )
  expect_error(filter_neighbourhoods(sce_mouseEmbryo),
               "x should be a Milo object. Please run `assign_neighbourhoods` first.",
               fixed=TRUE
  )
  sce_test = Milo(sce_mouseEmbryo)
  expect_error(filter_neighbourhoods(sce_test),
               "x should contain non-trivial graph. Please run `assign_neighbourhoods` first.",
               fixed=TRUE
  )
})

# return correct class
test_that("Return is the correct class", {
  # right class
  out = assign_neighbourhoods(x = sce_mouseEmbryo , k = 10, prop = 0.2, order = 2, filtering = FALSE, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
  out = filter_neighbourhoods(out)
  expect_s4_class(out, "Milo")
})



# filtering filtered data gives same result
test_that("Filtering filtered data is futile", {
  set.seed(15)
  sce_milo = assign_neighbourhoods(sce_mouseEmbryo , k = 25, filtering = TRUE, reducedDim_name = "pca.corrected")
  sce_milo_filtered = filter_neighbourhoods(sce_milo)

  # sce should be of the right format
  expect_identical(nhoods(sce_milo) , nhoods(sce_milo_filtered))
  expect_identical(as.numeric(nhoodIndex(sce_milo)) , as.numeric(nhoodIndex(sce_milo_filtered)))

})


# initial filtering = post hoc filtering
test_that("Filtering filtered data is futile", {
  set.seed(32)
  sce_milo_1 = assign_neighbourhoods(sce_mouseEmbryo , k = 25, filtering = TRUE, reducedDim_name = "pca.corrected")

  set.seed(32)
  sce_milo_2 = assign_neighbourhoods(sce_mouseEmbryo , k = 25, filtering = FALSE, reducedDim_name = "pca.corrected")
  sce_milo_2 = filter_neighbourhoods(sce_milo_2)

  # sce should be of the right format
  expect_identical(nhoods(sce_milo_1) , nhoods(sce_milo_2))
  expect_identical( as.numeric(nhoodIndex(sce_milo_1)) , as.numeric(nhoodIndex(sce_milo_2)))

})

