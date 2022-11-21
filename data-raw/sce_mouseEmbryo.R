## code to prepare `sce_mouseEmbryo` dataset goes here

usethis::use_data(sce_mouseEmbryo, overwrite = TRUE)


# load libraries
library(SingleCellExperiment)
library(MouseGastrulationData)
library(geneBasisR)

# subset of 3 cell types
cts = c("Spinal cord" , "Cardiomyocytes" , "NMP")

# load ref WT
samples = c(17,29,36,37)
sce_ref <- lapply(samples , function(sample){
  out <- EmbryoAtlasData(samples = sample)
  out <- out[, out$celltype %in% cts]
  return(out)
})
sce_ref <- do.call(cbind , sce_ref)
# add tal1 and type ids
sce_ref$type = "wt"
sce_ref$tal1 = "plus"


# load chimera Tal1
sce_chimera = Tal1ChimeraData()
# select CTs
sce_chimera = sce_chimera[, sce_chimera$celltype.mapped %in% cts]
# delete row for tomato
sce_chimera = sce_chimera[!rownames(sce_chimera) == "tomato-td" , ]

# add chimera_ prefix to avoid conclusion
colnames(sce_chimera) = paste0("chimera_" , colnames(sce_chimera) )
sce_chimera$cell = paste0("chimera_" ,sce_chimera$cell)
sce_chimera$sample = paste0("chimera_" , sce_chimera$sample )

# add type and tal1
sce_chimera$type = "chimera"
sce_chimera$tal1 = ifelse(sce_chimera$tomato , "minus" , "plus")

# sample same number of cells that is to sce
sce_chimera = sce_chimera[, sample(colnames(sce_chimera) , ncol(sce))]

# concatenate
## order rownames, check that dss are concatenable
sce_ref = sce_ref[order(rownames(sce_ref)) , ]
sce_chimera = sce_chimera[order(rownames(sce_chimera)) , ]
print(mean(rownames(sce_ref) == rownames(sce_chimera)))

## concatenate
sce = SingleCellExperiment(list(counts = cbind(counts(sce_ref) , counts(sce_chimera))))
## add metadata
sce$cell = c(sce_ref$cell , sce_chimera$cell)
sce$celltype = c(sce_ref$celltype , sce_chimera$celltype.mapped)
sce$sample = c(sce_ref$sample , sce_chimera$sample)
sce$type = c(sce_ref$type , sce_chimera$type)
sce$tal1 = c(sce_ref$tal1 , sce_chimera$tal1)

# add rowdata
rowData(sce) = rowData(sce_ref)

# add logcounts
sce = scuttle::logNormCounts(sce)


# add covariates
sce$stage = "E8.5"
sce$sex = sapply(1:ncol(sce) , function(i) ifelse(sce$sample[i] %in% c(36,37,"chimera_3" , "chimera_4") , "F" , "M"))
set.seed(32)
sce$toy_cov_1 = sample(1:5 , ncol(sce),1)

# select only 1000 genes
sce = retain_informative_genes(sce , n = 1000)
sce_mouseEmbryo = sce
usethis::use_data(sce_mouseEmbryo , overwrite = TRUE , compress = "xz")



