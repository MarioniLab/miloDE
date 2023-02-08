

#' calc_AUC_per_neighbourhood
#'
#' Returns per neighbourhood AUC from Augur based (RF) classifier
#'
#' @param x A \code{\linkS4class{Milo}} object.
#' @param genes Character vector specifying genes to be passed for the testing. Default \code{genes = rownames(x)}.
#' @param sample_id Character specifying which variable should be used as a replicate ID.
#' Should be in \code{colnames(colData(x))}. Default \code{sample_id = "sample"}.
#' @param condition_id Character specifying which variable should be used as a condition ID.
#' Should be in \code{colnames(colData(x))}.
#' @param conditions In case of multiple comparable groups, character vector specifying which conditions should be tested for separation.
#' Default \code{conditions = NULL} and assumes that only 2 different conditions are present.
#' @param min_n_cells_per_sample Positive integer specifying the minimum number of cells per replicate to be included in testing.
#' Default \code{min_n_cells_per_sample = 3}.
#' @param n_threads Positive integer specifying the number of cores to be used to calculate AUC.
#' Higher number results in faster calculation, but its feasibility depends on the specs of your machine.
#' Only relevant if \code{BPPARAM = NULL}. Default \code{n_threads = 2}.
#' @param BPPARAM NULL or \code{\link{MulticoreParam}} object. Default \code{BPPARAM = NULL} assuming no parallelisation.
#' @details
#' This function calculates for each neighbourhood whether cells between 2 conditions can be separated
#' with Random Forest based classifiers (adapted from \code{\link[Augur]{calculate_auc}}).
#' Accordingly, AUCs of the classifiers represent how well we can separate 2 conditions.
#'
#' We suggest that neighbourhoods with AUC > 0.5 suggest a certain degree of perturbation between 2 conditions that can further be examined
#' with DE testing. You also can set your own AUC threshold if desired as well as use AUCs to rank neighbourhoods.
#'
#' \emph{Note that this function is only relevant for \dQuote{simple} models (i.e. no interactions.)}.
#'
#' @return \code{data.frame} object, with AUC calculated for each neighbourhood.
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
#' sce = SingleCellExperiment(assays = list(counts =
#' floor(matrix(rnorm(n_row*n_col), ncol=n_col)) + 4))
#' logcounts(sce) = floor(matrix(rnorm(n_row*n_col), ncol=n_col)) + 4
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' sce$sample = floor(runif(n = n_col , min = 1 , max = 5))
#' sce$type = ifelse(sce$sample %in% c(1,2) , "ref" , "query")
#' reducedDim(sce , "reduced_dim") = matrix(rnorm(n_col*n_latent),
#' ncol=n_latent)
#' sce = assign_neighbourhoods(sce, reducedDim_name = "reduced_dim")
#' sce = calc_AUC_per_neighbourhood(sce, condition_id = "type")
calc_AUC_per_neighbourhood <- function(x , genes = rownames(x) , sample_id = "sample" ,
                                       condition_id , conditions = NULL,
                                       min_n_cells_per_sample = 3, n_threads = 2 , BPPARAM = NULL){

  out = .check_argument_correct(x, .check_sce, "Check x - something is wrong (are gene names unique?)") &
    .check_sce_milo(x) &
    .check_argument_correct(sample_id, is.character, "Check sample_id - should be character vector") &
    .check_argument_correct(condition_id, is.character, "Check condition_id - should be character vector") &
    .check_argument_correct(min_n_cells_per_sample, .check_positive_integer, "Check min_n_cells_per_sample - should be positive integer") &
    .check_argument_correct(n_threads , .check_positive_integer , "Check n_threads - should be positive integer") &
    .check_var_in_coldata_sce(x , sample_id , "sample_id") & .check_condition_in_coldata_sce(x , condition_id) &
    .check_sample_and_condition_id_valid(x , condition_id , sample_id) &
    .check_genes_in_sce(x , genes)

  if (!"logcounts" %in% assayNames(x)){
    stop("Please calculate log-normalised counts first if you want to calculate AUC per neighbourhood.")
  }

  x = x[genes, ]

  if (is.null(colnames(x))){
    colnames(x) = as.character(c(1:ncol(x)))
  }
  # assign condition and sample ids
  coldata <- as.data.frame(colData(x))
  x$milo_sample_id <- as.factor( coldata[, sample_id] )
  x$milo_condition_id <- as.factor( coldata[, condition_id] )

  if (is.null(conditions)){
    tab = table(x$milo_condition_id)
    if (!length(tab) == 2){
      stop("If conditions == NULL, there should be exactly two levels for tested conditions.")
    }
  } else {
    if (mean(conditions %in% unique(x$milo_condition_id)) < 1){
      stop("All specified conditions should be present.")
    }
    if (!length(conditions) == 2){
      stop("Conditions should be a vector of 2 elements.")
    }
    x = x[, x$milo_condition_id %in% conditions]
  }

  nhoods_sce = nhoods(x)
  # filter out for relevant cells
  current_cols = colnames(nhoods_sce)
  nhoods_sce = as.matrix( nhoods_sce[rownames(nhoods_sce) %in% colnames(x), ] )
  colnames(nhoods_sce) = current_cols
  # delete neighbourhoods that contain 0 cells (possible due to the upstream filtering)
  idx = which(colSums(nhoods_sce) > 0)
  current_cols = colnames(nhoods_sce)[idx]
  nhoods_sce = as.matrix( nhoods_sce[, idx] )
  colnames(nhoods_sce) = current_cols

  if (is.null(BPPARAM)){
    auc_stat = lapply(colnames(nhoods_sce) , function(hood_id){
      out = .get_auc_single_hood(x , nhoods_sce , hood_id , min_cells = 3 , min_n_cells_per_sample = min_n_cells_per_sample , n_threads = n_threads)
      return(out)
    })
  }
  else {
    auc_stat = bplapply(colnames(nhoods_sce) , function(hood_id){
      out = .get_auc_single_hood(x , nhoods_sce , hood_id , min_cells = 3 , min_n_cells_per_sample = min_n_cells_per_sample , n_threads = 1)
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
.get_auc_single_hood = function(x , nhoods_sce , hood_id , min_cells = 3 , min_n_cells_per_sample = 3 , n_threads = 2){
  out = .check_argument_correct(min_cells, .check_positive_integer, "Check min_cells - should be positive integer")
  # select cells
  current.cells = which(nhoods_sce[,hood_id] == 1)
  current.cells = rownames(nhoods_sce)[current.cells]
  current.sce = x[,colnames(x) %in% current.cells]
  current.sce = .filter_samples_with_low_n_cells_in_hood(current.sce , min_n_cells_per_sample = min_n_cells_per_sample)

  if (ncol(current.sce) > 0){
    current.sce$celltype.dummy = "dummy"
    meta = as.data.frame(colData(current.sce))
    tab = table(as.character(current.sce$milo_condition_id))
    if (length(tab) == 2 & tab[1] >= min_cells & tab[2] >= min_cells){
      auc = calculate_auc(logcounts(current.sce), meta, cell_type_col = "celltype.dummy",
                          label_col = "milo_condition_id" , n_subsamples = 0 ,
                          subsample_size = min_cells , min_cells = min_cells ,
                          feature_perc = 1 , n_threads = n_threads , show_progress = FALSE)
      out = as.data.frame(auc$AUC)
      out$auc_calculated = TRUE

    } else {
      out = data.frame(cell_type = "dummy" , auc = NaN , auc_calculated = FALSE)
    }
  } else {
    out = data.frame(cell_type = "dummy" , auc = NaN , auc_calculated = FALSE)
  }
  out$Nhood_center = hood_id
  return(out)
}




