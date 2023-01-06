

#' spatial_pval_adjustment
#'
#' Performs p-values multiple testing correction, with accounting for the overlap. This is achieved by using weighted version of BH correction.
#' @param nhoods_sce nhoods(sce)
#' @param pvalues vector of pvalues
#'
#' @return
#' @export
#'
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
#' reducedDim(sce , "reduced_dim") = matrix(rnorm(n_col*n_latent), ncol=n_latent)
#' sce = assign_neighbourhoods(sce, reducedDim_name = "reduced_dim" , k = 10 , order = 1)
#' nhoods_sce = nhoods(sce)
#' pvalues = runif(n = ncol(nhoods_sce) , min = 0 , max = 1)
#' out = spatial_pval_adjustment(nhoods_sce = nhoods_sce, pvalues = pvalues)
spatial_pval_adjustment = function(nhoods_sce , pvalues){

  # we will only calculate weights for neighbourhoods in which we are testing
  idx_not_nan = which(!is.na(pvalues))
  if (length(idx_not_nan) > 1){
    weights = .get_weights(as.matrix(nhoods_sce[,idx_not_nan]))

    out = .check_weights_and_pvals(weights , pvalues[idx_not_nan] , as.matrix(nhoods_sce[,idx_not_nan])) & .check_nhoods_matrix(nhoods_sce)


    n_comparisons = length(pvalues)
    #n_nans = sum(is.na(pvalues))
    pvalues = pvalues[idx_not_nan]

    # calc correction
    o <- order(pvalues)
    pvalues <- pvalues[o]
    weights <- weights[o]
    adjp <- numeric(length(o))
    adjp[o] <- rev(cummin(rev(sum(weights)*pvalues/cumsum(weights))))
    adjp <- pmin(adjp, 1)

    # add nans back
    adjp_total = rep(NaN, 1,n_comparisons)
    adjp_total[idx_not_nan] = adjp
  } else {
    adjp_total = pvalues
  }
  return(adjp_total)
}



.get_weights = function(nhoods_sce){

  out = .check_nhoods_matrix(nhoods_sce)

  intersect_mat <- crossprod(nhoods_sce)
  t.connect <- unname(rowSums(intersect_mat))
  weights<- 1/unlist(t.connect)
  return(weights)
}


