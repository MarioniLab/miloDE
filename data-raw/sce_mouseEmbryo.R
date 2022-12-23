## code to prepare `sce_mouseEmbryo` dataset goes here

usethis::use_data(sce_mouseEmbryo, overwrite = TRUE)


# load libraries
library(SingleCellExperiment)
library(MouseGastrulationData)
library(geneBasisR)

# subset of 3 cell types
cts = c("Endothelium" , "PGC", "Neural crest")

# load chimera Tal1
sce = Tal1ChimeraData()
# select CTs
sce = sce[, sce$celltype.mapped %in% cts]
# delete row for tomato
sce = sce[!rownames(sce) == "tomato-td" , ]

# add logcounts
sce = scuttle::logNormCounts(sce)


# add covariates
sce$sex = sapply(1:ncol(sce) , function(i) ifelse(sce$sample[i] %in% c("3" , "4") , "F" , "M"))
set.seed(32)
sce$toy_cov_1 = sample(1:5 , ncol(sce),1)

# select only 1000 genes
sce = retain_informative_genes(sce , n = 500)
sce_mouseEmbryo = SingleCellExperiment(list(counts = counts(sce)),
                                       colData = colData(sce),
                                       rowData = rowData(sce)
                                       )
reducedDim(sce_mouseEmbryo , "pca.corrected") = reducedDim(sce , "pca.corrected")
usethis::use_data(sce_mouseEmbryo , overwrite = TRUE , compress = "xz")



