# context("Testing rank_neighbourhoods_by_DE_magnitude")
library(miloDE)

# load data
data("sce_mouseEmbryo", package = "miloDE")
set.seed(32)
sce_milo = assign_neighbourhoods(sce_mouseEmbryo , k = 25, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
meta_sce = as.data.frame(colData(sce_milo))
de_stat = de_test_neighbourhoods(sce_milo , sample_id = "sample", design = ~tomato, covariates = c("tomato"), output_type = "SCE")


# wrong inputs return error
# error msgs
test_that("Wrong input gives errors", {
  require(SingleCellExperiment)
  # de_stat should be right format
  expect_error(rank_neighbourhoods_by_DE_magnitude(1),
               "de_stat should be either data.frame or SingleCellExperiment object.\n
         To get valid de_stat object, please run 'de_test_neighbourhoods.R'",
               fixed=TRUE
  )
  expect_error(rank_neighbourhoods_by_DE_magnitude(sce_milo),
               "de_stat should be either data.frame or SingleCellExperiment object.\n
         To get valid de_stat object, please run 'de_test_neighbourhoods.R'",
               fixed=TRUE
  )


  out = data.frame(Nhood = c(1:10) , Nhood_center = c(1:10))
  out$logFC = 0
  out$pval = 1
  expect_error(rank_neighbourhoods_by_DE_magnitude(out),
               "colnames(de_stat) missing some of the assay_names or coldata_names.",
               fixed=TRUE
  )

  out$pval_corrected_across_nhoods = 1
  out$pval_corrected_across_genes = 1
  out$Nhood = as.character(out$Nhood)
  expect_error(rank_neighbourhoods_by_DE_magnitude(out),
               "Nhood field should be numeric. To get valid de_stat object, please run 'de_test_neighbourhoods.R'",
               fixed=TRUE
  )


  out = SingleCellExperiment(list(logFC = assay(de_stat , "logFC") ,
                                  pval = assay(de_stat , "pval") ,
                                  pval_corrected_across_genes = assay(de_stat , "pval_corrected_across_genes") ,
                                  pval_corrected_across_nhoods = assay(de_stat , "pval_corrected_across_nhoods")
                                  ))
  rownames(out) = rownames(de_stat)
  out$Nhood = de_stat$Nhood
  expect_error(rank_neighbourhoods_by_DE_magnitude(out),
               "de_stat missing some of the coldata_names",
               fixed=TRUE
  )
  out$Nhood = as.character(out$Nhood)
  out$Nhood_center = out$Nhood
  expect_error(rank_neighbourhoods_by_DE_magnitude(out),
               "colData field 'Nhood' should be numeric. To get valid de_stat object, please run 'de_test_neighbourhoods.R'",
               fixed=TRUE
  )


  # pvalue - positive number between 0 and 1
  expect_error(rank_neighbourhoods_by_DE_magnitude(de_stat , pval.thresh = "a"),
               "pval.thresh should be numeric",
               fixed=TRUE
  )
  expect_error(rank_neighbourhoods_by_DE_magnitude(de_stat , pval.thresh = FALSE),
               "pval.thresh should be numeric",
               fixed=TRUE
  )
  expect_error(rank_neighbourhoods_by_DE_magnitude(de_stat , pval.thresh = c(1:10)),
               "pval.thresh should be a single number",
               fixed=TRUE
  )
  expect_error(rank_neighbourhoods_by_DE_magnitude(de_stat , pval.thresh = -1),
               "pval.thresh should be between 0 and 1",
               fixed=TRUE
  )
  expect_error(rank_neighbourhoods_by_DE_magnitude(de_stat , pval.thresh = 0),
               "pval.thresh should be between 0 and 1",
               fixed=TRUE
  )
  expect_error(rank_neighbourhoods_by_DE_magnitude(de_stat , pval.thresh = 1),
               "pval.thresh should be between 0 and 1",
               fixed=TRUE
  )


  # z-thresh - neagtive number
  expect_error(rank_neighbourhoods_by_DE_magnitude(de_stat , z.thresh = "a"),
               "z.thresh should be numeric",
               fixed=TRUE
  )
  expect_error(rank_neighbourhoods_by_DE_magnitude(de_stat , z.thresh = FALSE),
               "z.thresh should be numeric",
               fixed=TRUE
  )
  expect_error(rank_neighbourhoods_by_DE_magnitude(de_stat , z.thresh = c(1:10)),
               "z.thresh should be a single number",
               fixed=TRUE
  )
  expect_error(rank_neighbourhoods_by_DE_magnitude(de_stat , z.thresh = 0.01),
               "z.thresh should be not higher than 0",
               fixed=TRUE
  )

})


test_that("right format", {
  out = rank_neighbourhoods_by_DE_magnitude(de_stat)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out) , ncol(nhoods(sce_milo)))

  nhood_stat = data.frame(Nhood = c(1:ncol(nhoods(sce_milo))) , Nhood_center = colnames(nhoods(sce_milo)))
  nhood_stat$Nhood_center = as.character(nhood_stat$Nhood_center)

  out = out[, c("Nhood" , "Nhood_center")]
  out$Nhood_center = as.character(out$Nhood_center)
  out = out[order(out$Nhood) , ]
  expect_identical(out$Nhood , nhood_stat$Nhood)
  expect_identical(out$Nhood_center , nhood_stat$Nhood_center)
})



test_that("NaNs give 0s", {
  # de_stat should be right format
  test = convert_de_stat(de_stat)
  test$pval = NaN
  test$pval_corrected_across_nhoods = NaN
  test$pval_corrected_across_genes = NaN

  out = rank_neighbourhoods_by_DE_magnitude(test)
  expect_s3_class(out, "data.frame")
  expect_equal(1 , mean(out$n_DE_genes == 0))
  expect_equal(1 , mean(out$n_specific_DE_genes == 0))
})



test_that("lower pvals - less genes", {
  out_1 = rank_neighbourhoods_by_DE_magnitude(de_stat , pval.thresh = 0.01)
  out_2 = rank_neighbourhoods_by_DE_magnitude(de_stat , pval.thresh = 0.1)
  out_3 = rank_neighbourhoods_by_DE_magnitude(de_stat , pval.thresh = 0.5)
  expect_gte(mean(out_3$n_DE_genes) , mean(out_2$n_DE_genes))
  expect_gte(mean(out_2$n_DE_genes) , mean(out_1$n_DE_genes))
})

test_that("lower z - less genes", {
  out_1 = rank_neighbourhoods_by_DE_magnitude(de_stat , z.thresh = -1)
  out_2 = rank_neighbourhoods_by_DE_magnitude(de_stat , z.thresh = -2)
  out_3 = rank_neighbourhoods_by_DE_magnitude(de_stat , z.thresh = -3)
  expect_gte(mean(out_1$n_specific_DE_genes) , mean(out_2$n_specific_DE_genes))
  expect_gte(mean(out_2$n_specific_DE_genes) , mean(out_3$n_specific_DE_genes))
})


