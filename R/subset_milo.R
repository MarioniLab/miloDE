

#' subset_milo
#'
#' Function properly subsets data from sce-milo object (including nhoods(sce-milo))
#' @param sce_milo Milo object
#' @param colnames_2_retain Character string of colnames to retain
#'
#' @return
#' @export
#' @importFrom miloR buildNhoodGraph nhoodIndex<- nhoods nhoods<-
#' @importFrom igraph V V<- induced_subgraph
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
#' reducedDim(sce , "reduced_dim") = matrix(rnorm(n_col*n_latent), ncol=n_latent)
#' sce = assign_neighbourhoods(sce, reducedDim.name = "reduced_dim" , k = 10 , order = 1)
#' sce = subset_milo(sce , colnames_2_retain = c("1", "2", "3"))
subset_milo = function(sce_milo , colnames_2_retain){
  sce_milo_filtered = sce_milo[, colnames(sce_milo) %in% colnames_2_retain]

  # get subgraph
  V(miloR::graph(sce_milo))$name = colnames(sce_milo)
  miloR::graph(sce_milo_filtered) = induced_subgraph(miloR::graph(sce_milo) , colnames_2_retain)

  # subset rows from nhoods
  nhoods_sce = nhoods(sce_milo_filtered)
  nhoods_sce = nhoods_sce[rownames(nhoods_sce) %in% colnames_2_retain , ]

  # check if we have empty hoods - exclude them
  hood_sizes = colSums(nhoods_sce)
  idx_2_keep = which(hood_sizes > 0)
  nhoods_sce = nhoods_sce[, idx_2_keep]

  # update graph and hood index
  nhoods(sce_milo_filtered) = nhoods_sce
  sce_milo_filtered = buildNhoodGraph(sce_milo_filtered)
  nhoodIndex(sce_milo_filtered) = nhoodIndex(sce_milo_filtered)[match(colnames(nhoods(sce_milo_filtered)) , nhoodIndex(sce_milo_filtered))]
  return(sce_milo_filtered)
}

