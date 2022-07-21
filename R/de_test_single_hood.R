


#' de_test_single_hood
#'
#' Tests single hood for DE. Not intended to be used by itself (however possible), but rather as a part of `de_test_all_hoods`
#' @param sce SingleCellExperiment object
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
#' @importFrom SingleCellExperiment colData SingleCellExperiment counts assay
#' @import edgeR
#' @importFrom tibble rownames_to_column
#' @importFrom scuttle summarizeAssayByGroup
#' @examples
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
   if (gene_selection == "per_hood"){
     keep <- filterByExpr(y, group=summed$condition.id , min.count = min.count , min.total.count = round(min.count * 1.5))
     y <- y[keep,]
   }
   y <- calcNormFactors(y)
   if (!is.null(covariates)){
    design <- model.matrix(as.formula( paste("~ ", paste(covariates, collapse="+"),sep = "" , " + condition.id") ) , y$samples)
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
