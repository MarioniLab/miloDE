

#' spatial_pval_adjustment
#'
#' Performs p-values multiple testing correction across neighbourhoods, with accounting for the overlap
#' @param nhoods_x Should be extracted from x as \code{nhoods(x)}.
#' @param pvalues Vector of p-values.
#' @details
#' This function is not intended to be run by itself, but can be useful if the user wants to perform \sQuote{spatially aware} multiple testing correction.
#'
#' Under the hood it performs weighted version of BH correction, where weights are reciprocal to the local desnity of a neighbourhood.
#' Accordingly, big and/or highly connected neighbourhoods will have lower weights.
#'
#' @return Vector with \sQuote{spatially} (i.e. across neighbourhoods) adjusted p-values
#' @export
#' @examples
#' require(SingleCellExperiment)
#' require(miloR)
#' n_row = 500
#' n_col = 100
#' n_latent = 5
#' sce = SingleCellExperiment(assays = list(counts =
#' floor(matrix(rnorm(n_row*n_col), ncol=n_col)) + 4))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' reducedDim(sce , "reduced_dim") =
#' matrix(rnorm(n_col*n_latent), ncol=n_latent)
#' sce = assign_neighbourhoods(sce,
#' reducedDim_name = "reduced_dim" , k = 10 , order = 1)
#' nhoods_x = nhoods(sce)
#' pvalues = runif(n = ncol(nhoods_x) , min = 0 , max = 1)
#' out = spatial_pval_adjustment(nhoods_x, pvalues = pvalues)
spatial_pval_adjustment = function(nhoods_x , pvalues){

  # we will only calculate weights for neighbourhoods in which we are testing
  idx_not_nan = which(!is.na(pvalues))
  if (length(idx_not_nan) > 1){
    weights = .get_weights(as.matrix(nhoods_x[,idx_not_nan]))

    out = .check_weights_and_pvals(weights , pvalues[idx_not_nan] , as.matrix(nhoods_x[,idx_not_nan])) & .check_nhoods_matrix(nhoods_x)

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



.get_weights = function(nhoods_x){
  out = .check_nhoods_matrix(nhoods_x)
  intersect_mat <- crossprod(nhoods_x)
  t.connect <- unname(rowSums(intersect_mat))
  weights<- 1/unlist(t.connect)
  return(weights)
}


