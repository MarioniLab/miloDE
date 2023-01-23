

#' filter_neighbourhoods
#'
#' Filtering redundant hoods, using the greedy approach to set cover problem
#' @param sce_milo A Milo object (SCE + assigned hoods)
#'
#' @return
#' @export
#' @importFrom RcppGreedySetCover greedySetCover
#' @importFrom miloR buildNhoodGraph nhoodIndex nhoods graph graph<-
#' @examples
#' require(SingleCellExperiment)
#' n_row = 500
#' n_col = 100
#' n_latent = 5
#' sce = SingleCellExperiment(assays = list(counts = floor(matrix(rnorm(n_row*n_col), ncol=n_col)) + 4))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' reducedDim(sce , "reduced_dim") = matrix(rnorm(n_col*n_latent), ncol=n_latent)
#' sce = assign_neighbourhoods(sce, reducedDim_name = "reduced_dim" , k = 10 , order = 1)
#' sce = filter_neighbourhoods(sce)
filter_neighbourhoods = function(sce_milo){

  #args = c(as.list(environment()))
  #out = .general_check_arguments(args)
  out = .check_argument_correct(sce_milo, .check_sce_milo, "Check sce_milo - something is wrong. Calculate 'assign_neighbourhoods' first.)")

  nhoods_sce = nhoods(sce_milo)
  stat_hoods = lapply(1:ncol(nhoods_sce) , function(i){
    current.cells = which(nhoods_sce[,i] == 1)
    out = data.frame(set = colnames(nhoods_sce)[i],
                     element = names(current.cells))
    return(out)
  })
  stat_hoods = do.call(rbind, stat_hoods)
  stat_filtered = suppressMessages( greedySetCover(stat_hoods) )
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
