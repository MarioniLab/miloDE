# context("Testing convert_de_stat")
library(miloDE)
library(miloR)
library(SingleCellExperiment)

# load data
data("sce_mouseEmbryo", package = "miloDE")
set.seed(32)
sce = assign_neighbourhoods(sce_mouseEmbryo , k = 25, prop = 0.2, order = 2, filtering = T, reducedDim_name = "pca.corrected", k_init = 50, d = 30)
meta_sce = as.data.frame(colData(sce))
de_stat_sce = de_test_neighbourhoods(sce , sample_id = "sample", design = ~tomato, covariates = c("tomato"), output_type = "SCE")
de_stat_df = de_test_neighbourhoods(sce , sample_id = "sample", design = ~tomato, covariates = c("tomato"),  output_type = "data.frame")


# wrong inputs return error
# error msgs
test_that("Wrong input gives errors", {

  # de_stat should be of the right format
  expect_error(convert_de_stat(de_stat = 1),
               "de_stat should be either data.frame or SingleCellExperiment object.\n
         To get valid de_stat object, please run 'de_test_neighbourhoods.R'",
               fixed=TRUE
  )
  expect_error(convert_de_stat(de_stat = de_stat_df , assay_names = c("logFC") , coldata_names = c("logFC")),
               "assay_names and coldata_names can not overlap",
               fixed=TRUE
  )
  expect_error(convert_de_stat(de_stat = de_stat_df , assay_names = c("test")),
               "colnames(de_stat) missing some of the assay_names or coldata_names.",
               fixed=TRUE
  )
  expect_error(convert_de_stat(de_stat = de_stat_df , coldata_names = c("test")),
               "colnames(de_stat) missing some of the assay_names or coldata_names.",
               fixed=TRUE
  )
  expect_error(convert_de_stat(de_stat = de_stat_sce , coldata_names = c("test")),
               "de_stat missing some of the coldata_names",
               fixed=TRUE
  )
  expect_error(convert_de_stat(de_stat = de_stat_sce , assay_names = c("test")),
               "de_stat missing some of the required assays.",
               fixed=TRUE
  )
})



# right input returns the expected format (in default output)
test_that("Return is the correct class", {
  # right class
  expect_s3_class(convert_de_stat(de_stat_sce) , "data.frame")
  expect_s4_class(convert_de_stat(de_stat_df) , "SingleCellExperiment")

  # right assay names, sce to df
  de_stat_test = convert_de_stat(de_stat_sce)
  out = mean(c("logFC" , "pval" , "pval_corrected_across_genes" , "pval_corrected_across_nhoods" ,
               "Nhood" , "Nhood_center" , "test_performed") %in% colnames(de_stat_test))
  expect_equal(1 , out)

  n_genes = length(unique(de_stat_test$gene))
  expect_equal(n_genes , nrow(de_stat_sce))
  n_hoods = length(unique(de_stat_test$Nhood))
  expect_equal(n_hoods , ncol(de_stat_sce))


  # right assay names, df to sce
  de_stat_test = convert_de_stat(de_stat_df)
  n_genes = length(unique(de_stat_df$gene))
  expect_equal(n_genes , nrow(de_stat_test))
  n_hoods = length(unique(de_stat_df$Nhood))
  expect_equal(n_hoods , ncol(de_stat_test))

})


## expect identical output if conversion is done post-hoc
test_that("Conversion ad hoc is equal to conversion post hoc", {

  expect_identical(de_stat_df , convert_de_stat(de_stat_sce))

  de_stat_sce_2 = convert_de_stat(de_stat_df)
  de_stat_sce = de_stat_sce[order(rownames(de_stat_sce)) , ]
  de_stat_sce_2 = de_stat_sce_2[order(rownames(de_stat_sce_2)) , ]

  expect_identical(assay(de_stat_sce, "logFC") , assay(de_stat_sce_2 , "logFC"))
  expect_identical(assay(de_stat_sce, "pval") , assay(de_stat_sce_2 , "pval"))
  expect_identical(assay(de_stat_sce, "pval_corrected_across_genes") , assay(de_stat_sce_2, "pval_corrected_across_genes"))
  expect_identical(assay(de_stat_sce, "pval_corrected_across_nhoods") , assay(de_stat_sce_2 , "pval_corrected_across_nhoods"))
  expect_identical(as.data.frame(colData(de_stat_sce)) , as.data.frame(colData(de_stat_sce_2)))

})


## expect conversion of right variables (coldata and assays) to the right format
test_that("Conversion of additional metadata", {

  de_stat_df_test = de_stat_df
  de_stat_df_test$celltype = 0
  de_stat_df_test$celltype[de_stat_df_test$Nhood %in% c(1,2,5)] = 1

  de_stat_sce_test = convert_de_stat(de_stat_df_test , coldata_names = "celltype")
  meta_test = as.data.frame(colData(de_stat_sce_test))

  out = rep(0,1,ncol(de_stat_sce_test))
  out[c(1,2,5)] = 1
  expect_identical(meta_test$celltype , out)

})


## expect error if coldata is not legit
test_that("Conversion of additional metadata legit only if metadata is legit", {

  de_stat_df_test = de_stat_df
  de_stat_df_test$celltype = 0
  unq_genes = unique(de_stat_df_test$gene)
  de_stat_df_test$celltype[de_stat_df_test$gene == unq_genes[1]] = 1

  expect_error(convert_de_stat(de_stat_df_test , coldata_names = "celltype"))
  expect_error(convert_de_stat(de_stat_df_test , assay_names = "celltype") , NA)

  de_stat_df_test = de_stat_sce
  assay(de_stat_df_test , "logFC_2") = assay(de_stat_sce , "logFC")
  out_1 = convert_de_stat(de_stat_df_test)
  out_2 = convert_de_stat(de_stat_df_test , assay_names = "logFC_2")

  expect_equal(FALSE , "logFC_2" %in% colnames(out_1))
  expect_equal(TRUE , "logFC_2" %in% colnames(out_2))
})



