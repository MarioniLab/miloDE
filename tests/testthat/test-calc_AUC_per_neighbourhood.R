# context("Testing calc_AUC_per_neighbourhood")
library(miloDE)

# load data
data("sce_mouseEmbryo", package = "miloDE")
set.seed(32)
sce_mouseEmbryo_milo = assign_neighbourhoods(sce_mouseEmbryo, reducedDim_name = "pca.corrected")
require(miloR)
stat_auc = calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , condition_id = "tomato" , min_n_cells_per_sample = 1)
stat_auc_w_annotated_hoods =  annotateNhoods(sce_mouseEmbryo_milo , da.res = stat_auc , coldata_col = "celltype.mapped")

# generate toy data
require(SingleCellExperiment)
n_row = 50
n_col = 100
n_latent = 5
sce = SingleCellExperiment(list(counts = floor(matrix(rnorm(n_row*n_col), ncol=n_col)) + 4  ,
                           logcounts = floor(matrix(rnorm(n_row*n_col), ncol=n_col)) + 4) )
rownames(sce) = as.factor(1:n_row)
colnames(sce) = c(1:n_col)
sce$cell = colnames(sce)
sce$sample = c(1:ncol(sce))
reducedDim(sce , "reduced_dim") = matrix(rnorm(n_col*n_latent), ncol=n_latent)
sce = assign_neighbourhoods(sce, reducedDim_name = "reduced_dim")

sce_1 = sce
sce_1$condition = ifelse(logcounts(sce_1[1,]) >= 5 , 1 , 0)
sce_2 = sce
sce_2$condition = sample(c(0,1) , ncol(sce_2) , replace = T)

sce_3 = sce
sce_3$condition = sapply(c(1:ncol(sce_3)) , function(i){
  if (logcounts(sce_3[1,i]) <= 2){
    return(0)
  } else if (logcounts(sce_3[1,i]) <= 4){
    return(1)
  } else {
    return(2)
  }
})
sce_4 = sce
sce_4$condition = sample(c(0,1,2) , ncol(sce_4) , replace = T)




# error msgs
test_that("Wrong input gives errors", {
  # x should be of the right format
  expect_error(calc_AUC_per_neighbourhood(x = 1 , condition_id = "test"),
               "x should be a SingleCellExperiment or Milo object.",
               fixed=TRUE
  )
  sce_test = sce_mouseEmbryo_milo
  colnames(sce_test) = rep(1,1,ncol(sce_test))
  expect_error(calc_AUC_per_neighbourhood(sce_test , condition_id = "test2"),
               "If colnames(x) exist, they should be unique.",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo , condition_id = "test3"),
               "x should be a Milo object. Please run `assign_neighbourhoods` first.",
               fixed=TRUE
  )


  # genes should be a subset of sce's rows
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , genes = c("1" , "11") , condition_id = "tomato"),
               "Some gene names are missing from x",
               fixed=TRUE
  )

  # sample_id should be of the right format
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , condition_id = "tomato" , sample_id = 1),
               "Check sample_id - should be character vector",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , condition_id = "tomato" , sample_id = TRUE),
               "Check sample_id - should be character vector",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , condition_id = "tomato" , sample_id = "test"),
               "'sample_id' should be in colData(x)",
               fixed=TRUE
  )


  # condition_id should be of the right format
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , sample_id = "sample" , condition_id = "sample" ),
               "'sample_id' and 'condition_id' can not be the same",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , sample_id = "sample" , condition_id = 1 ),
               "Check condition_id - should be character vector",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , sample_id = "sample" , condition_id = "test" ),
               "'condition_id' should be in colData(x)",
                fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , sample_id = "sample" , condition_id = "stage" ),
               "x should have at least two levels for tested conditions.",
               fixed=TRUE
  )


  # conditions are of the right format
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , sample_id = "sample" ,
                                          condition_id = "tomato" , conditions = c(1,"2") ),
               "All specified conditions should be present.",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(sce_4 , sample_id = "sample" ,
                                          condition_id = "condition" , conditions = NULL , min_n_cells_per_sample = 1),
               "If conditions == NULL, there should be exactly two levels for tested conditions.",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(sce_3 , sample_id = "sample" ,
                                          condition_id = "condition" , conditions = NULL , min_n_cells_per_sample = 1),
               "If conditions == NULL, there should be exactly two levels for tested conditions.",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(sce_3 , sample_id = "sample" ,
                                          condition_id = "condition" , conditions = c(0,1,2) , min_n_cells_per_sample = 1),
               "Conditions should be a vector of 2 elements.",
               fixed=TRUE
  )


  # min_n_cells_per_sample - positive integer
  expect_error(calc_AUC_per_neighbourhood(x = sce_mouseEmbryo_milo , sample_id = "sample", condition_id = "tomato", min_n_cells_per_sample = "aa"),
               "Check min_n_cells_per_sample - should be positive integer",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(x = sce_mouseEmbryo_milo , sample_id = "sample", condition_id = "tomato", min_n_cells_per_sample = 1.5),
               "Check min_n_cells_per_sample - should be positive integer",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , sample_id = "sample", condition_id = "tomato", min_n_cells_per_sample = 0),
               "Check min_n_cells_per_sample - should be positive integer",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , sample_id = "sample", condition_id = "tomato", min_n_cells_per_sample = -100),
               "Check min_n_cells_per_sample - should be positive integer",
               fixed=TRUE
  )

  # n_threads - positive integer
  expect_error(calc_AUC_per_neighbourhood(x = sce_mouseEmbryo_milo , sample_id = "sample", condition_id = "tomato", n_threads = "aa"),
               "Check n_threads - should be positive integer",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(x = sce_mouseEmbryo_milo , sample_id = "sample", condition_id = "tomato", n_threads = 1.5),
               "Check n_threads - should be positive integer",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , sample_id = "sample", condition_id = "tomato", n_threads = 0),
               "Check n_threads - should be positive integer",
               fixed=TRUE
  )
  expect_error(calc_AUC_per_neighbourhood(sce_mouseEmbryo_milo , sample_id = "sample", condition_id = "tomato", n_threads = -100),
               "Check n_threads - should be positive integer",
               fixed=TRUE
  )
})




# right input returns the expected format (in default output)
test_that("Return is the correct class", {
  # right class
  expect_s3_class(stat_auc , "data.frame")

  out = calc_AUC_per_neighbourhood(sce_1 , condition_id = "condition" , min_n_cells_per_sample = 1)
  expect_s3_class(out , "data.frame")
  out = calc_AUC_per_neighbourhood(sce_2 , condition_id = "condition" , conditions = c(0,1))
  expect_s3_class(out , "data.frame")
  out = calc_AUC_per_neighbourhood(sce_3 , condition_id = "condition" , conditions = c(0,1))
  expect_s3_class(out , "data.frame")
  out = calc_AUC_per_neighbourhood(sce_4 , condition_id = "condition" , conditions = c(2,1) , min_n_cells_per_sample = 1)
  expect_s3_class(out , "data.frame")


  # right dimensions
  expect_equal(ncol(nhoods(sce_mouseEmbryo_milo)) , nrow(stat_auc))
  expect_equal(4, ncol(stat_auc))
  cols = c("Nhood" , "Nhood_center" , "auc" , "auc_calculated")
  expect_identical(cols , colnames(stat_auc))

})



# check that AUCs are approx of what we expect
## real mouse data
test_that("AUCs make sense for mouse", {

  # AUCs for BPs are either NaNs or > 0.5
  aucs = stat_auc_w_annotated_hoods$auc[stat_auc_w_annotated_hoods$celltype.mapped == "Blood progenitors 2"]
  expect_equal(0, length(which(aucs < 0.5)))

  # AUCs for Endo are all > 0.5 (and not nan)
  aucs = stat_auc_w_annotated_hoods$auc[stat_auc_w_annotated_hoods$celltype.mapped == "Endothelium"]
  expect_equal(0, sum(is.na(aucs)))
  expect_equal(length(aucs), sum(aucs > 0.5))
})


## simulations - 2 conditions
test_that("AUCs make sense for sim1", {
  stat_1 = calc_AUC_per_neighbourhood(sce_1 , condition_id = "condition" , min_n_cells_per_sample = 1)
  stat_1_gene_1 = calc_AUC_per_neighbourhood(sce_1 , genes = c("1" , "2") , condition_id = "condition" , min_n_cells_per_sample = 1)
  stat_1_min_cells_3 = calc_AUC_per_neighbourhood(sce_1 , condition_id = "condition" , min_n_cells_per_sample = 3)

  stat_2 = calc_AUC_per_neighbourhood(sce_2 , condition_id = "condition" , min_n_cells_per_sample = 1)
  stat_2_gene_1 = calc_AUC_per_neighbourhood(sce_2 , genes = c("1" , "2") , condition_id = "condition" , min_n_cells_per_sample = 1)
  stat_2_min_cells_3 = calc_AUC_per_neighbourhood(sce_2 , condition_id = "condition" , min_n_cells_per_sample = 3)

  expect_gt(stat_1$auc , stat_2$auc)
  expect_gt(stat_1$auc , 0.5)
  expect_gt(stat_1_gene_1$auc , stat_1$auc)
  expect_true(is.nan(stat_1_min_cells_3$auc))
  expect_true(is.nan(stat_2_min_cells_3$auc))
  expect_true(!is.nan(stat_2$auc))
  expect_true(!is.nan(stat_2_gene_1$auc))
})


## simulations - 3 conditions
test_that("AUCs make sense for sim2", {

  stat_3_12 = calc_AUC_per_neighbourhood(sce_3 , genes = c("1" , "2") , condition_id = "condition" , min_n_cells_per_sample = 1 , conditions = c(1,2))
  stat_3_01 = calc_AUC_per_neighbourhood(sce_3 , genes = c("1" , "2") , condition_id = "condition" , min_n_cells_per_sample = 1 , conditions = c(1,0))
  stat_3_20 = calc_AUC_per_neighbourhood(sce_3 , genes = c("1" , "2") , condition_id = "condition" , min_n_cells_per_sample = 1 , conditions = c(0,2))

  stat_4_12 = calc_AUC_per_neighbourhood(sce_4 , genes = c("1" , "2") , condition_id = "condition" , min_n_cells_per_sample = 1 , conditions = c(1,2))
  stat_4_01 = calc_AUC_per_neighbourhood(sce_4 , genes = c("1" , "2") , condition_id = "condition" , min_n_cells_per_sample = 1 , conditions = c(1,0))
  stat_4_20 = calc_AUC_per_neighbourhood(sce_4 , genes = c("1" , "2") , condition_id = "condition" , min_n_cells_per_sample = 1 , conditions = c(0,2))

  expect_gte(stat_3_20$auc , stat_3_12$auc)
  expect_gte(stat_3_20$auc , stat_3_01$auc)

  expect_gt(stat_3_12$auc , stat_4_12$auc)
  expect_gt(stat_3_01$auc , stat_4_20$auc)
  expect_gt(stat_3_20$auc , stat_4_01$auc)

})






