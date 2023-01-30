#' estimate_neighbourhood_sizes
#'
#' For a grid of k (order fixed), return hood size distribution; that will help a user to select appropriate k
#' @param sce SCE object
#' @param reducedDim_name defines the slot in reducedDim(sce) to use as embedding for graph construction
#' @param k_grid Vector of positive integers, defines how many neighbours to use for the hood assignment
#' @param prop Numerical, between 0 and 1, defines fraction of cells from SCE to use for the hoods. Default = 0.1.
#' @param order In {1,2}, defines which order of neighbours to use
#' @param filtering In {TRUE,FALSE}, defines whether to filter hoods (reduces computing time greatly). Default = TRUE
#' @param reducedDim_name defines the slot in reducedDim(sce) to use as embedding for graph construction
#' @param k_init Positive integer, defines how many neighbours to use for identifying anchor cells
#' @param d Positive integer, defines how many dimensions from reducedDim(sce) to use
#' @param cluster_id Character specifying which field in colData(sce) to use for 'localised' neighbourhood size estimation. This might be useful in case dataset is
#' rather big which will result in excessive running time. In case cluster_id is provided, we will calculate neighbourhood size distribution within individual clusters and aggregate results
#' across clusters in order to speed up the process (note that it might result in slightly biased estimates). Default is NULL, in which case neighbourhood sizes will be estimated for the whole dataset.
#' @param plot_stat Boolean specifying whether to plot the stat
#' @param verbose Boolean specifying whether to print intermediate output messages. Default = FALSE.
#'
#' @return
#' @export
#' @importFrom miloR Milo buildGraph graph<- graph nhoods<- nhoodIndex<- buildNhoodGraph
#' @importFrom tibble rownames_to_column
#' @importFrom SummarizedExperiment colData
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom stats quantile
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
#' out = estimate_neighbourhood_sizes(sce, k_grid = c(5,10), reducedDim_name = "reduced_dim")
#'
estimate_neighbourhood_sizes = function(sce, reducedDim_name , k_grid = seq(10,100,10) , order = 2, prop = 0.1 , filtering = TRUE,
                                        k_init = 50 , d = 30 , cluster_id = NULL, plot_stat = TRUE , verbose = FALSE){

  out = .check_argument_correct(sce, .check_sce, "Check sce - something is wrong (gene names unique? reducedDim.name is not present?)") &
    .check_argument_correct(k_grid, is.numeric, "Check k_grid - should be numeric vector") &
    .check_argument_correct(prop, .check_prop, "Check prop - should be positive number between 0 and 1") &
    .check_argument_correct(order, function(x) .check_arg_within_options(x, c(1,2)),
                            "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)") &
    .check_argument_correct(filtering, .check_boolean, "Check filtering - should be either TRUE or FALSE") &
    .check_argument_correct(reducedDim_name, is.character, "Check reducedDim_name - should be character vector") &
    .check_argument_correct(k_init, .check_positive_integer, "Check k_init - should be positive integer") &
    .check_argument_correct(d, .check_positive_integer, "Check d - should be positive integer") &
    .check_reducedDim_in_sce(sce , reducedDim_name) & .check_k_grid(k_grid)

  # check that cluster_id is in colData(sce)
  if (!is.null(cluster_id)){
    if (!cluster_id %in% colnames(colData(sce))){
      stop("If cluster_id not NULL, it should be in colnames(colData(sce))")
    }
  }

  quantile_vec = seq(0,1,0.25)
  # check that k_grid reasonable -- at least 2 values, the the highest is smaller than 1000;
  # otherwise warn
  k_grid = sort(unique(k_grid))
  message(paste0("Running for next k values:\n" , paste(k_grid , collapse = ", ")))

  if (is.null(cluster_id)){
    stat = lapply(k_grid , function(k){
      sce_milo = assign_neighbourhoods(sce , k = k , prop = prop , order = order , filtering = filtering,
                                       reducedDim_name = reducedDim_name , k_init = k_init , d = d , verbose = FALSE)
      out = quantile(colSums(nhoods(sce_milo)) , probs = quantile_vec)
      return(out)
    })
  } else {
    meta = as.data.frame(colData(sce))
    clusters = table( meta[, cluster_id] )
    # select only big clusters
    clusters = names(clusters)[clusters > 2*max(k_grid)]
    if (length(clusters) == 0){
      stop("All specified clusters have # cells < 2*max(k). We recommed to provide higher clustering resolution, decreasing max(k) or set cluster_id = NULL.")
    }
    else {
      stat = lapply(k_grid , function(k){
        stat_per_k = sapply(clusters , function(cluster){
          idx = which(meta[, cluster_id] == cluster)
          sce_milo = assign_neighbourhoods(sce[,idx] , k = k , prop = prop , order = order , filtering = filtering,
                                           reducedDim_name = reducedDim_name , k_init = k_init , d = d , verbose = FALSE)
          out = colSums(nhoods(sce_milo))
          return(out)
        })
        stat_per_k = unlist(stat_per_k)
        out = quantile(stat_per_k , probs = quantile_vec)
        return(out)
      })
    }
  }

  stat = as.data.frame( do.call(rbind , stat) )
  rownames(stat) = k_grid
  stat = rownames_to_column(stat , var = "k")
  message(paste0("Finished the estimation of neighbourhood sizes ~ k dependancy (order = " , order , ")."))
  colnames(stat) = c("k" , "min" , "q25" , "med" , "q75" , "max")
  stat$k = factor(stat$k , levels = k_grid)

  if (plot_stat){
    p_stat = tryCatch(
      {
        p = ggplot(stat, aes(k)) +
          geom_boxplot( aes(ymin = min, lower = q25, middle = med, upper = q75, ymax = max , fill = k), stat = "identity") +
          theme_bw() +
          scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(length(k_grid))) +
          labs( y = "Neighbourhood size")
        p
      },
      error=function(err){
        warning("Can not return plot")
        return(NULL)
      }
    )
    print(p_stat)
  }
  return(stat)
}

