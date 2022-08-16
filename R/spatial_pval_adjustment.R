

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
spatial_pval_adjustment = function(weights = NULL , nhoods_sce , pvalues){
  # get weights
  if (is.null(weights)){
    intersect_mat <- crossprod(nhoods_sce)
    t.connect <- unname(rowSums(intersect_mat))
    weights<- 1/unlist(t.connect)
  }

  n_comparisons = length(pvalues)
  # lets take the not NaNs
  idx_not_nan = which(!is.na(pvalues))
  #n_nans = sum(is.na(pvalues))
  pvalues = pvalues[idx_not_nan]
  weights = weights[idx_not_nan]

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
  return(adjp_total)
}


