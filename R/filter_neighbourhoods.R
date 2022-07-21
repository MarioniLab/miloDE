

#' filter_neighbourhoods
#'
#' Filtering redundant hoods, using the greedy approach to set cover problem
#' @param sce_milo A Milo object (SCE + assigned hoods)
#'
#' @return
#' @export
#' @importFrom RcppGreedySetCover greedySetCover
#' @importFrom miloR buildNhoodGraph nhoodIndex nhoods
#' @examples
#'
filter_neighbourhoods = function(sce_milo){
  require(RcppGreedySetCover)
  nhoods_sce = nhoods(sce_milo)
  stat_hoods = lapply(1:ncol(nhoods_sce) , function(i){
    current.cells = which(nhoods_sce[,i] == 1)
    out = data.frame(set = colnames(nhoods_sce)[i],
                     element = names(current.cells))
    return(out)
  })
  stat_hoods = do.call(rbind, stat_hoods)
  stat_filtered = greedySetCover(stat_hoods)
  hoods_filtered = unique(stat_filtered$set)
  if (length(hoods_filtered) > 1){
    nhoods_sce = nhoods_sce[, colnames(nhoods_sce) %in% hoods_filtered]
  }
  else {
    nhoods_sce = as.matrix(nhoods_sce[, colnames(nhoods_sce) %in% hoods_filtered])
    colnames(nhoods_sce) = hoods_filtered
  }
  nhoods(sce_milo) = nhoods_sce
  sce_milo = buildNhoodGraph(sce_milo)
  nhoodIndex(sce_milo) = nhoodIndex(sce_milo)[match(colnames(nhoods(sce_milo)) , nhoodIndex(sce_milo))]
  return(sce_milo)
}
