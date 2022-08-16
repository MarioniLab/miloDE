


#' filter_not_perturbed_hoods
#'
#' Pre-filtering hoods where we expect to have no DE (using classifier based method Augur that returns AUCs. AUC = 0.5)
#' @param sce sce_milo A Milo object (SCE + assigned hoods)
#' @param min_cells A positive integer specifying min number of cells in the hood for which we calculate AUC. Default = 3 (avoid changing this parameter).
#' @param auc.thresh A positive numeric specifying an AUC threshold below each we assign hoods to be no DE. Default = 0.5 (avoid changing this parameter)
#'
#' @return
#' @export
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom miloR buildNhoodGraph nhoodIndex nhoods
#'
#' @examples
filter_not_perturbed_hoods <- function(sce , min_cells = 3 , auc.thresh = 0.5){
  # assign condition and sample ids
  coldata <- as.data.frame(colData(sce))
  sce$condition.id <- as.factor( coldata[, condition.id] )
  sce$sample.id <- as.factor( coldata[, sample.id] )

  nhoods_sce = nhoods(sce)

  auc_stat = lapply(colnames(nhoods_sce) , function(hood.id){
    out = .get_auc_single_hood(sce , nhoods_sce , hood.id , min_cells = min_cells)
    return(out)
  })
  auc_stat = do.call(rbind , auc_stat)
  auc_stat$Nhood_id = colnames(nhoods_sce)

  hoods_keep = auc_stat$Nhood_id[auc_stat$auc > auc.thresh]

  if (length(hoods_keep) >= 1){
    if (length(hoods_keep) > 1){
      nhoods_sce = nhoods_sce[, colnames(nhoods_sce) %in% hoods_keep]
    }
    else if (length(hoods_keep) == 1){
      nhoods_sce = as.matrix(nhoods_sce[, colnames(nhoods_sce) %in% hoods_keep])
      colnames(nhoods_sce) = hoods_keep
    }
    nhoods(sce_milo) = nhoods_sce
    sce_milo = buildNhoodGraph(sce_milo)
    nhoodIndex(sce_milo) = nhoodIndex(sce_milo)[match(colnames(nhoods(sce_milo)) , nhoodIndex(sce_milo))]
    return(sce_milo)
  }
  else {
    message("All hoods have AUC less than selected threshold. No DE to be performed.")
  }

}


#' @importFrom SingleCellExperiment logcounts
#' @importFrom Augur calculate_auc
.get_auc_single_hood = function(sce , nhoods_sce , hood.id , min_cells = 3){

  # select cells
  current.cells = which(nhoods_sce[,hood.id] == 1)
  current.sce = sce[,current.cells]
  current.sce = .filter_samples_with_low_n_cells_in_hood(current.sce , min_n_cells_per_sample = min_n_cells_per_sample)

  current.sce$celltype.dummy = "dummy"
  meta = as.data.frame(colData(current.sce))
  if (ncol(current.sce) > 0){
    auc = calculate_auc(logcounts(current.sce), meta, cell_type_col = "celltype.dummy",
                        label_col = "condition.id" , n_subsamples = 0 ,
                        subsample_size = min_cells , min_cells = min_cells ,
                        feature_perc = 1)
    out = as.data.frame(auc$AUC)
    return(out)
  }
}
