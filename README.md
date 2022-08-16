# milo_DE
Framework for sensitive DE testing (using neighbourhoods).

miloDE builds on existing framework for DA testing called [Milo](https://pubmed.ncbi.nlm.nih.gov/34594043/). 
It exploits the notion of overlapping neighborhoods (henceforth refered to as 'hoods') of homogeneous cells, constructed from graph-representation of scRNA-seq datasets, and performs testing within each neighborhood.



### Installation

```
## Install development version
library(devtools)
devtools::install_github("MarioniLab/miloDE") 
```


### Main steps

1. We need to define which latent space will be used for graph construction. It is either provided apriori (and should be a field in reducedDim(SCE)) or we can calculate one using `add-embedding`. For construction of latent space we can use either `Azimuth` or `MNN`(needs to be specified in `reduction_type`). 
Other essential parameters are: `reducedDim.name`, `cells_ref` - specifying colnames that correspond to reference cells, `cells_query` - specifying colnames that correspond to query cells.


```
library(miloDE)
sce = add_embedding(sce , 
                    reduction_type = "MNN", 
                    reducedDim.name = "pca.corrected" , 
                    cells_ref = colnames(sce)[sce$type == "reference"],
                    cells_query = colnames(sce)[sce$type == "query"])

```

2. We need to assign hoods to our SCE object using `assign_neighbourhoods`. We allow a standard milo assignment and a more optimal (i.e. minimise redundancy) assignment by specifying parameter `filtering = T`. Note that the optimal assignment can be also carried out post hoc using `filter_neighbourhoods`.

```

sce = assign_neighbourhoods(sce , reducedDim.name = "pca.corrected")

```

3. Now we can carry out DE testing for each hood using `de_test_all_hoods`. That returns table for each hood x gene with calculated statsitcis: logFC, Pvalue, corrected Pvalue across genes (FDR) and corrected Pvalue across hoods (SpatialFDR).


```

sce = de_test_all_hoods(sce , reducedDim.name = "pca.corrected" , condition.id = "type")

```

