



#' de_test_all_hoods
#'
#' Tests all hoods for DE + performs spatialFDR correction after
#' @param sce Milo object
#' @param genes Subset of rownames(sce) for which we will perform DE testing. Default = rownames(sce)
#' @param sample.id Character specifying which variable should be used as a sample/replica id. Should be in colData(sce)
#' @param condition.id Character specifying which variable should be used as a condition id. Should be in colData(sce)
#' @param discard_not_perturbed_hoods Boolean specifying whether perform prior hood selection using Augur based classifier. Note for big datasets it might take a while so recommended if a large portion of hoods is expected to be unperturbed.
#' @param covariates Vector specifying if additional covariates should be passed into experimental design. Default = NULL (no covariates)
#' @param min_n_cells_per_sample positive integer specifying the minimun number of cells per replica to be included in testing. Default = 2
#' @param gene_selection In {"all" , "per_hood"}. If == "per_hood", it will further discard genes not relevant for this particular hood.
#' If == "all", it will perform overall selection and then test each hood across same genes.
#' @param min.count Positive integer, specifying min.count for gene selection. Default = 5
#'
#' @return
#' @export
#' @importFrom SingleCellExperiment colData SingleCellExperiment counts
#' @importFrom SummarizedExperiment assay
#' @import edgeR
#' @importFrom tibble rownames_to_column
#' @importFrom scuttle summarizeAssayByGroup
#'
#' @examples
#' require(SingleCellExperiment)
#' n_row = 500
#' n_col = 100
#' n_latent = 5
#' sce = SingleCellExperiment(assays = list(counts = floor(matrix(rnorm(n_row*n_col), ncol=n_col)) + 4))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' sce$sample = floor(runif(n = n_col , min = 1 , max = 5))
#' sce$type = ifelse(sce$sample %in% c(1,2) , "ref" , "query")
#' reducedDim(sce , "reduced_dim") = matrix(rnorm(n_col*n_latent), ncol=n_latent)
#' sce = assign_hoods(sce, reducedDim.name = "reduced_dim")
#' de_stat = de_test_all_hoods(sce , condition.id = "type" , discard_not_perturbed_hoods = FALSE)
de_test_all_hoods = function(sce ,
                             genes = rownames(sce) ,
                             sample.id = "sample",
                             condition.id ,
                             discard_not_perturbed_hoods = c(TRUE,FALSE),
                             covariates = NULL ,
                             min_n_cells_per_sample = 2,
                             gene_selection = "all" ,
                             min.count = 2){

  sce = sce[genes , ]
  # assign condition and sample ids
  coldata <- as.data.frame(colData(sce))
  sce$condition.id <- as.factor( coldata[, condition.id] )
  sce$sample.id <- as.factor( coldata[, sample.id] )
  if (discard_not_perturbed_hoods){
    sce = .filter_not_perturbed_hoods(sce)
  }
  if (!isFALSE(sce)){
    nhoods_sce = nhoods(sce)

    if (gene_selection == "all"){
      summed = summarizeAssayByGroup(counts(sce), colData(sce)[,c("condition.id", "sample.id")])
      summed = SingleCellExperiment(list(counts=assay(summed, "sum")), colData=colData(summed))
      y <- DGEList(counts(summed), samples=colData(summed), lib.size = colSums(counts(summed)))
      keep <- filterByExpr(y, group=summed$sample.id , min.count = min.count , min.total.count = round(min.count * 1.5))
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

    unq.genes = unique(de_stat$gene)

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
  else {
    message("No hoods were tested.")
    return(F)
  }
}



#' de_test_single_hood
#'
#' Tests single hood for DE. Not intended to be used by itself (however possible), but rather as a part of `de_test_all_hoods`
#' @param sce Milo object
#' @param genes Subset of rownames(sce) for which we will perform DE testing. Default = rownames(sce)
#' @param nhoods_sce Can be extracted from sce as nhoods(sce) prior to running the function. The reason the argument is passed by itself is to avoid calculating it every time.
#' @param hood.id Character specifying for which hood we should perform testing. Should be in colnames(nhoods_sce)
#' @param sample.id Character specifying which variable should be used as a sample/replica id. Should be in colData(sce)
#' @param condition.id Character specifying which variable should be used as a condition id. Should be in colData(sce)
#' @param covariates Vector specifying if additional covariates should be passed into experimental design. Default = NULL (no covariates)
#' @param min_n_cells_per_sample positive integer specifying the minimun number of cells per replica to be included in testing. Default = 2
#' @param gene_selection In {"all" , "per_hood"}. If = "per_hood", we will further discard genes not relevant for this particular hood.
#' @param min.count Positive integer, specifying min.count for gene selection. Default = 5
#'
#' @return
#' @export
#' @importFrom SingleCellExperiment colData SingleCellExperiment counts
#' @importFrom SummarizedExperiment assay
#' @import edgeR
#' @importFrom tibble rownames_to_column
#' @importFrom scuttle summarizeAssayByGroup
#' @examples
#' require(SingleCellExperiment)
#' require(miloR)
#' n_row = 500
#' n_col = 100
#' n_latent = 5
#' sce = SingleCellExperiment(assays = list(counts = floor(matrix(rnorm(n_row*n_col), ncol=n_col)) + 4))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' sce$sample = floor(runif(n = n_col , min = 1 , max = 5))
#' sce$type = ifelse(sce$sample %in% c(1,2) , "ref" , "query")
#' reducedDim(sce , "reduced_dim") = matrix(rnorm(n_col*n_latent), ncol=n_latent)
#' sce = assign_hoods(sce, reducedDim.name = "reduced_dim")
#' nhoods_sce = nhoods(sce)
#' de_stat = de_test_single_hood(sce , nhoods_sce = nhoods_sce, hood.id = colnames(nhoods_sce)[1] , sample.id = "sample" , condition.id = "type")
de_test_single_hood = function(sce , genes = rownames(sce) , nhoods_sce , hood.id , sample.id , condition.id , covariates = NULL,
                               min_n_cells_per_sample = 2 , gene_selection = "all" , min.count = 5){
  sce = sce[genes , ]
  coldata <- as.data.frame(colData(sce))

  # assign condition and sample ids
  sce$condition.id <- as.factor( coldata[, condition.id] )
  sce$sample.id <- as.factor( coldata[, sample.id] )

  # select cells
  current.cells = which(nhoods_sce[,hood.id] == 1)
  current.sce = sce[,current.cells]
  current.sce = .filter_samples_with_low_n_cells_in_hood(current.sce , min_n_cells_per_sample = min_n_cells_per_sample)

  tab = table(current.sce$sample.id , current.sce$condition.id)

  if (ncol(tab) < 2){
    stat = .return_null_stat(genes)
  }
  else if (sum(tab[,1] > 0) == 0 | sum(tab[,2] > 0) == 0 | sum(tab > 0) <= 3){
    stat = .return_null_stat(genes)
  }
  else {
    summed = summarizeAssayByGroup(counts(current.sce), colData(current.sce)[,c("condition.id", "sample.id" , covariates)])
    summed = SingleCellExperiment(list(counts=assay(summed, "sum")), colData=colData(summed))
    y <- DGEList(counts(summed), samples=colData(summed), lib.size = colSums(counts(summed)))
    # make sure covariates have contrasts; delete covariates that don't
    meta_y = as.data.frame(y$samples)

    if (!is.null(covariates)){
      covariates_2_keep = sapply(covariates , function(covariate){
        tab = table(meta_y[, covariate] , meta_y$condition.id)
        if (nrow(tab) == 1 | ncol(tab) == 1 | max(rowMins(tab)) == 0){
          return(F)
        }
        else {
          return(T)
        }
      })
      covariates_2_keep = names(covariates_2_keep)[covariates_2_keep]
      if (length(covariates_2_keep) == 0){
        covariates_2_keep = NULL
      }
    }
    else {
      covariates_2_keep = NULL
    }


    if (gene_selection == "per_hood"){
      keep <- filterByExpr(y, group=summed$sample.id , min.count = min.count , min.total.count = round(min.count * 1.5))
      y <- y[keep,]
    }
    y <- calcNormFactors(y)
    if (!is.null(covariates_2_keep)){
      design <- model.matrix(as.formula( paste("~ ", paste(covariates_2_keep, collapse="+"),sep = "" , " + condition.id") ) , y$samples)
    }
    else {
      design <- model.matrix(~condition.id , y$samples)
    }
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    res <- glmQLFTest(fit, coef=ncol(design))
    stat = as.data.frame(topTags(res, n = Inf ) )
    stat = stat[order(rownames(stat)) , c("logFC" , "logCPM" , "PValue" , "FDR")]
    stat = rownames_to_column(stat , var = "gene")
  }
  stat$Nhood_id = hood.id
  return(stat)
}


#' @importFrom SingleCellExperiment colData
#' @importFrom miloR buildNhoodGraph nhoodIndex nhoods
#'
.filter_not_perturbed_hoods <- function(sce){

  # hardcoded for now - maybe change later
  auc.thresh = 0.5
  min_cells = 3

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
    message(paste0(length(hoods_keep) , " hoods are retained for further analysis."))
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
    stop("All hoods have AUC less than selected threshold. No DE to be performed.")
    return(F)
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


#'
.filter_samples_with_low_n_cells_in_hood = function(sce , min_n_cells_per_sample = 2){
  tab = table(sce$sample.id)
  samples_2_keep = names(tab)[tab > min_n_cells_per_sample]
  sce = sce[, sce$sample.id %in% samples_2_keep]
  return(sce)
}

#'
.return_null_stat = function(genes){
  genes = genes[order(genes)]
  stat = data.frame(logFC = NaN, logCPM = NaN, PValue = NaN, FDR = NaN)
  stat = stat[rep(seq_len(nrow(stat)), each = length(genes)), ]
  stat$gene = genes
  stat = stat[, c("gene", "logFC" , "logCPM" , "PValue" , "FDR")]
  return(stat)
}


