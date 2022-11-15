#' estimate_hood_sizes
#'
#' For a grid of k (order fixed), return hood size distribution; that will help a user to select appropriate k
#' @param sce SCE object
#' @param k.grid Vector of positive integers, defines how many neighbours to use for the hood assignment
#' @param prop Numerical, between 0 and 1, defines fraction of cells from SCE to use for the hoods. Default = 0.1.
#' @param order In {1,2}, defines which order of neighbours to use
#' @param filtering In {T,F}, defines whether to filter hoods (reduces computing time greatly). Default = TRUE
#' @param reducedDim.name defines the slot in reducedDim(sce) to use as embedding for graph construction
#' @param k_init Positive integer, defines how many neighbours to use for identifying anchor cells
#' @param d Positive integer, defines how many dimensions from reducedDim(sce) to use
#' @param quantile.vec Vector of distribution quantiles to be estimated. Default = seq(0,1,0.25)
#'
#' @return
#' @export
#' @importFrom miloR Milo buildGraph graph<- graph nhoods<- nhoodIndex<- buildNhoodGraph
#' @importFrom tibble rownames_to_column
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
#' out = estimate_hood_sizes(sce, k.grid = c(5,10), reducedDim.name = "reduced_dim")
#'
estimate_hood_sizes = function(sce, k.grid = seq(10,100,10) , order = 2, prop = 0.1 , filtering = T,
                               reducedDim.name , k_init = 50 , d = 30 , quantile.vec = seq(0 , 1 , 0.25)){

  args = c(as.list(environment()))
  out = .general_check_arguments(args) & .check_reducedDim_in_sce(sce , reducedDim.name)

  # check that k.grid reasonable -- at least 2 values, the the highest is smaller than 1000;
  # otherwise warn
  k.grid = sort(unique(k.grid[k.grid > 2]))
  if (length(k.grid) == 0){
    stop("No acceptable k values were entered. Input appropriate k.grid (should be positive integers)")
    return(F)
  }
  else {
    message("Accepted k values:\n" )
    message(c( paste0( sapply(k.grid[1:length(k.grid) - 1] , function(x) paste0(x , ", "))) , k.grid[length(k.grid)]))

    if (length(k.grid) == 1){
      warning("You only selected one value for k. If it is intended, we recommend to run directly 'assign_hoods'")
    }
    if (max(k.grid) >= 1000){
      warning("The highest selected value is > 1000. It is gonna cost computationally, and we generally do not recommend
              such high k. Consider reducing.")
    }
    stat = lapply(k.grid , function(k){
      sce_milo = assign_hoods(sce , k = k , prop = prop , order = order , filtering = filtering,
                                       reducedDim.name = reducedDim.name , k_init = k_init , d = d)
      out = .get_stat_single_coverage(nhoods(sce_milo) , quantile.vec)
      return(out)
    })

    stat = as.data.frame( do.call(rbind , stat) )
    rownames(stat) = k.grid
    stat = rownames_to_column(stat , var = "k")
    message(paste0("Finished the estimation of hood sizes ~ k dependancy (order = " , order , ")."))
    print(stat)
    return(stat)
  }
}


#' @importFrom stats quantile
.get_stat_single_coverage = function(nhoods_sce , quantile.vec){
  vec = colSums(nhoods_sce)
  out = quantile(vec , probs = quantile.vec)
  return(out)
}

