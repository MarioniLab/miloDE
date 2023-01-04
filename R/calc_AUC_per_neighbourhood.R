

#' calc_AUC_per_neighbourhood
#'
#' Returns per neighbourhood AUC from Augur-classifier
#'
#' @param sce Milo object
#'
#' @param sample_id Character specifying which variable should be used as a sample/replica id. Should be in colData(sce).
#' @param genes Character vector specifying genes to be passed for the testing
#' @param condition_id Character specifying which variable should be used as a condition id. Should be in colData(sce).
#' @param min_n_cells_per_sample Positive integer specifying the minimum number of cells per replica to be included in testing. Default = 2.
#' @param n_threads Positive integer specifying the number of cores to be used to calculate AUC. Higher number results in faster calculation, but its feasibility depends on the specs of your machine. Only relevant if BPPARAM = NULL.
#' @param BPPARAM NULL or MulticoreParam object
#'
#' @return
#' @export
#' @importFrom SummarizedExperiment colData assayNames
#' @importFrom miloR nhoods
#' @importFrom BiocParallel bplapply
#'
#' @examples
#' require(SingleCellExperiment)
#' n_row = 500
#' n_col = 100
#' n_latent = 5
#' sce = SingleCellExperiment(assays = list(counts = floor(matrix(rnorm(n_row*n_col), ncol=n_col)) + 4))
#' logcounts(sce) = floor(matrix(rnorm(n_row*n_col), ncol=n_col)) + 4
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' sce$sample = floor(runif(n = n_col , min = 1 , max = 5))
#' sce$type = ifelse(sce$sample %in% c(1,2) , "ref" , "query")
#' reducedDim(sce , "reduced_dim") = matrix(rnorm(n_col*n_latent), ncol=n_latent)
#' sce = assign_neighbourhoods(sce, reducedDim_name = "reduced_dim")
#' sce = calc_AUC_per_neighbourhood(sce, sample_id = "sample" , condition_id = "type")
#'
calc_AUC_per_neighbourhood <- function(sce , genes = rownames(sce) , sample_id = "sample" , condition_id , min_n_cells_per_sample = 1, n_threads = 2 , BPPARAM = NULL){

  out = .check_argument_correct(sce, .check_sce, "Check sce - something is wrong (gene names unique? reducedDim.name is not present?)") &
    .check_sce_milo(sce) &
    .check_argument_correct(sample_id, is.character, "Check sample_id - should be character vector") &
    .check_argument_correct(condition_id, is.character, "Check condition_id - should be character vector") &
    .check_argument_correct(min_n_cells_per_sample, .check_positive_integer, "Check min_n_cells_per_sample - should be positive integer") &
    .check_argument_correct(n_threads , .check_positive_integer , "Check n_threads - should be positive integer") &
    .check_var_in_coldata_sce(sce , sample_id , "sample_id") & .check_condition_in_coldata_sce(sce , condition_id) &
    .check_sample_and_condition_id_valid(sce , condition_id , sample_id) &
    .check_genes_in_sce(sce , genes)

  if (!"logcounts" %in% assayNames(sce)){
    stop("Please calculate log-normalised counts first if you want to calculate AUC per neighbourhood.")
  }

  sce = sce[genes, ]
  # assign condition and sample ids
  coldata <- as.data.frame(colData(sce))
  sce$condition_id <- as.factor( coldata[, condition_id] )
  sce$sample_id <- as.factor( coldata[, sample_id] )

  nhoods_sce = nhoods(sce)

  if (is.null(BPPARAM)){
    auc_stat = lapply(colnames(nhoods_sce) , function(hood_id){
      out = .get_auc_single_hood(sce , nhoods_sce , hood_id , min_cells = 3 , min_n_cells_per_sample = min_n_cells_per_sample , n_threads = n_threads)
      return(out)
    })
  }
  else {
    auc_stat = bplapply(colnames(nhoods_sce) , function(hood_id){
      out = .get_auc_single_hood(sce , nhoods_sce , hood_id , min_cells = 3 , min_n_cells_per_sample = min_n_cells_per_sample , n_threads = 1)
      return(out)
    } , BPPARAM = BPPARAM)
  }
  auc_stat = do.call(rbind , auc_stat)
  # add Nhood
  meta_nhoods = data.frame(Nhood = 1:ncol(nhoods_sce) , Nhood_center = colnames(nhoods_sce))
  auc_stat = merge(auc_stat , meta_nhoods , all.x = TRUE , all.y = FALSE)
  auc_stat = auc_stat[, c("Nhood", "Nhood_center" , "auc" , "auc_calculated")]
  auc_stat = auc_stat[order(auc_stat$Nhood) , ]
  return(auc_stat)
}



#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment colData
#' @import Augur
.get_auc_single_hood = function(sce , nhoods_sce , hood_id , min_cells = 3 , min_n_cells_per_sample = 1 , n_threads = 2){

  out = .check_argument_correct(min_cells, .check_positive_integer, "Check min_cells - should be positive integer")
  # select cells
  current.cells = which(nhoods_sce[,hood_id] == 1)
  current.sce = sce[,current.cells]
  current.sce = .filter_samples_with_low_n_cells_in_hood(current.sce , min_n_cells_per_sample = min_n_cells_per_sample)

  current.sce$celltype.dummy = "dummy"
  meta = as.data.frame(colData(current.sce))

  if (ncol(current.sce) > 0){
    tab = table(current.sce$condition_id)
    if (length(tab) == 2 & tab[1] >= min_cells & tab[2] >= min_cells){
      auc = calculate_auc(logcounts(current.sce), meta, cell_type_col = "celltype.dummy",
                          label_col = "condition_id" , n_subsamples = 0 ,
                          subsample_size = min_cells , min_cells = min_cells ,
                          feature_perc = 1 , n_threads = n_threads , show_progress = FALSE)
      out = as.data.frame(auc$AUC)
      out$auc_calculated = TRUE
    }
    else {
      out = data.frame(cell_type = "dummy" , auc = NaN , auc_calculated = FALSE)
    }
  }
  else {
    out = data.frame(cell_type = "dummy" , auc = NaN , auc_calculated = FALSE)
  }
  out$Nhood_center = hood_id
  return(out)
}

