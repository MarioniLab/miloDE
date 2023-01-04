# milo_DE
Framework for sensitive DE testing (using neighbourhoods).

miloDE builds on existing framework for DA testing called [Milo](https://pubmed.ncbi.nlm.nih.gov/34594043/). 
It exploits the notion of overlapping neighborhoods (nhoods) of homogeneous cells, constructed from graph-representation of scRNA-seq data, and performs testing within each neighborhood. Multiple testing correction is performed either across nhoods or across genes. 

In addition to DE testing, we provide functional to rank neighbourhoods by degree of the DE as well as plotting functions to visualise results. In vignette, we showcase how you can carry out clustering analysis to group genes in co-regulated transcriptional modules.



## Installation

```
# Install development version:
library(devtools)
devtools::install_github("MarioniLab/miloDE") 
library(miloDE)

## If you plan to use parallelisation (desired for big datasets), please install BiocParallel:
BiocManager::install("BiocParallel")


```




## Pipeline

1. *Input*. Input of `miloDE` is scRNA-seq data, provided as `SingleCellExperiment` object. 
Additionally, we require that:

* Latent dimension (used for graph construction) is pre-computed and stored in `reducedDim(sce)`.
* `colData(sce)` has to contain metadata corresponding to individual replicates (passed to `sample_id`) and tested condition (e.g. healthy or disease, passed to `condition_id`).

You can explore toy data here:

```
data("sce_mouseEmbryo", package = "miloDE")
print(sce_mouseEmbryo)
# `pca.corrected` in reducedDim(sce) - batch corrected PCs that we will use for graph construction

head(colData(sce_mouseEmbryo))
# `tomato` corresponds to condition id  
# `sample` corresponds to individual replicates. There are 2 samples per each condition:
table(sce_mouseEmbryo$sample , sce_mouseEmbryo$sample$tomato)


```

2. *Neighbourhood assignment*: First step is to assign neighbourhoods using graph representation of scRNA-seq data'

```

sce_mouseEmbryo = assign_neighbourhoods(sce_mouseEmbryo, k = 25, order = 2, filtering = TRUE, reducedDim_name = "pca.corrected")

```

2. *DE testing*: Once neighbourhoods are assigned, we can carry out DE testing. Output is returned in either `data.frame` or `SingleCellExperiment format`. For each tested gene-nhood, we return `logFC`, `pvalue`, `pvalue corrected across genes` and `pvalue corrected across nhoods`. We also return two boolean flags (for each nhood): 

* `sufficient_n_samples` indicates whether testing can be carried out given sample composition within a nhood (i.e. do we have cells from enough samples for both conditions); 
* `design_matrix_suitable` indicates whether testing can be carried out with included covariates. 

In other words, `sufficient_n_samples` reflects whether testing can be carried out without covariates, and `design_matrix_suitable` reflects whether testing can be carried out with covariates.



```

de_stat = de_test_neighbourhoods = function(sce_mouseEmbryo , sample_id = "sample", condition_id = "tomato",
                                  gene_selection = "none", output_type = "data.frame" )


```


Please check the vignette to grasp on additional functions aiding interpretation and analysis of miloDE output.


## Vignette

@@ Add link to vignette here

