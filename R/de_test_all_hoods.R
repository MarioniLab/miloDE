



#' de_test_all_hoods
#'
#' Tests all hoods for DE + performs spatialFDR correction after
#' @param sce SingleCellExperiment object
#' @param genes Subset of rownames(sce) for which we will perform DE testing. Default = rownames(sce)
#' @param nhoods_sce Can be extracted from sce as nhoods(sce) prior to running the function. The reason the argument is passed by itself is to avoid calculating it every time.
#' @param hood.id Character specifying for which hood we should perform testing. Should be in colnames(nhoods_sce)
#' @param sample.id Character specifying which variable should be used as a sample/replica id. Should be in colData(sce)
#' @param condition.id Character specifying which variable should be used as a condition id. Should be in colData(sce)
#' @param covariates Vector specifying if additional covariates should be passed into experimental design. Default = NULL (no covariates)
#' @param min_n_cells_per_sample positive integer specifying the minimun number of cells per replica to be included in testing. Default = 2
#' @param gene_selection In {"all" , "per_hood"}. If == "per_hood", it will further discard genes not relevant for this particular hood.
#' If == "all", it will perform overall selection and then test each hood across same genes.
#' @param min.count Positive integer, specifying min.count for gene selection. Default = 5
#'
#' @return
#' @export
#' @importFrom SingleCellExperiment colData SingleCellExperiment counts assay
#' @import edgeR
#' @importFrom tibble rownames_to_column
#' @importFrom scuttle summarizeAssayByGroup
#'
#' @examples
de_test_all_hoods = function(sce , genes = rownames(sce) , sample.id , condition.id , covariates = NULL , min_n_cells_per_sample = 2, gene_selection = "all" , min.count){

  sce = sce[genes , ]
  # assign condition and sample ids
  coldata <- as.data.frame(colData(sce))
  sce$condition.id <- as.factor( coldata[, condition.id] )
  sce$sample.id <- as.factor( coldata[, sample.id] )

  nhoods_sce = nhoods(sce)

  if (gene_selection == "all"){
    summed = summarizeAssayByGroup(counts(sce), colData(sce)[,c("condition.id", "sample.id")])
    summed = SingleCellExperiment(list(counts=assay(summed, "sum")), colData=colData(summed))
    y <- DGEList(counts(summed), samples=colData(summed), lib.size = colSums(counts(summed)))
    keep <- filterByExpr(y, group=summed$type , min.count = min.count , min.total.count = round(min.count * 1.5))
    genes = names(keep)[keep]
    sce = sce[genes , ]
  }

  # get de stat for all hoods
  de_stat = lapply(colnames(nhoods_sce) , function(hood.id){
    out = de_test_single_hood(sce , genes = genes , nhoods_sce = nhoods_sce, hood.id = hood.id,
                        sample.id = sample.id, condition.id = condition.id, covariates = covariates,
                        min_n_cells_per_sample = min_n_cells_per_sample ,
                        gene_selection = gene_selection , min.count = min.count)
    return(out)
  })
  de_stat = do.call(rbind , de_stat)

  # add Nhood
  nhood_meta = data.frame(Nhood = 1:ncol(nhoods_sce) , Nhood_id = colnames(nhoods_sce))
  de_stat = merge(de_stat , nhood_meta , all.x = T)

  # add spatial correction
  ## get weights from nhoods_sce
  intersect_mat <- crossprod(nhoods_sce)
  t.connect <- unname(rowSums(intersect_mat))
  weights<- 1/unlist(t.connect)

  unq.genes = unique(de_stat_selected_hoods$gene)

  de_stat_w_correction = lapply(unq.genes , function(gene){
    current.stat = de_stat[de_stat$gene == gene , ]
    current.stat = current.stat[order(current.stat$Nhood) , ]
    current.stat$SpatialFDR = spatial_pval_adjustment(weights = weights , nhoods_sce = nhoods_sce , pvalues = current.stat$PValue)
    current.stat$gene = gene
    return(current.stat)
  })
  de_stat_w_correction = do.call(rbind , de_stat_w_correction)
  return(de_stat_w_correction)
}

