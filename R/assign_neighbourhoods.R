

#' assign_neighbourhoods
#'
#' Assign neighbourhoods to \code{SingleCellExperiment} object
#' @param x A \code{\linkS4class{SingleCellExperiment}} object.
#' @param reducedDim_name Defines the assay in \code{reducedDim(x)} to use as the embedding for graph construction.
#' @param k Positive integer, defines how many neighbours to use for the neighbourhood assignment. Default \code{k = 25}.
#' @param prop Numerical, between 0 and 1, defines which fraction of cells to use as neighbourhood centres. Default \code{prop = 0.2}.
#' @param order In \code{c(1,2)}, defines which order of neighbours to use. Default \code{order = 2}.
#' @param filtering In \code{c(TRUE,FALSE)}, defines whether to refine the assignment. Default \code{filtering = TRUE}.
#' @param k_init Positive integer, defines how many neighbours to use for identifying anchor cells (for this step we use 1st-order kNN). Default \code{k_init = 50}.
#' @param d Positive integer, defines how many dimensions from \code{reducedDim(x)} to use. Default \code{d = 30}.
#' @param verbose Boolean specifying whether to print intermediate output messages. Default \code{verbose = TRUE}.
#' @details
#' This function assigns neighbourhoods to single-cell data. This includes assigning graph representation, selecting \sQuote{index} cells and, finally, for each index cell, assigning it along with its neighbourhoors to one neighbourhood.
#'
#' Specifically, algorithm goes as follows:
#' 1. Assigning \sQuote{loose} graph (i.e. ~low k, 1st-order kNN) to select index cells for the selected \code{prop} (greatly reduces computational time to look for \sQuote{index} cells in a loose graph).
#' 2. Reassigning graph following entered by the user \code{order} and \code{k}.
#' 3. Assigning neighbourhoods.
#' 4. (Optional but recommended) Refining the neighbourhood assignment (check \code{\link{filter_neighbourhoods}}).
#'
#' @return Milo object containing cell-neighbourhood matrix in \code{nhoods(out)} slot.
#' @export
#' @importFrom miloR Milo buildGraph graph<- graph nhoods<- nhoodIndex<- buildNhoodGraph
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assay
#' @import Matrix
#' @importFrom igraph connect V neighborhood
#' @importFrom methods as
#' @examples
#' require(SingleCellExperiment)
#' n_row = 500
#' n_col = 100
#' n_latent = 5
#' sce = SingleCellExperiment(assays =
#' list(counts = floor(matrix(rnorm(n_row*n_col), ncol=n_col)) + 4))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' reducedDim(sce , "reduced_dim") =
#' matrix(rnorm(n_col*n_latent), ncol=n_latent)
#' out = assign_neighbourhoods(sce, reducedDim_name = "reduced_dim")
assign_neighbourhoods = function(x , reducedDim_name , k = 25, prop = 0.2, order = 2, filtering = TRUE, k_init = 50, d = 30, verbose = TRUE){

  #args = c(as.list(environment()))
  #out = .general_check_arguments(args) & .check_reducedDim_in_sce(sce , reducedDim_name)
  out = .check_argument_correct(x, .check_sce, "Check x - something is wrong (gene names unique?)") &
        .check_argument_correct(k, .check_positive_integer, "Check k - should be positive integer") &
        .check_argument_correct(prop, .check_prop, "Check prop - should be positive number between 0 and 1") &
        .check_argument_correct(order, function(x) .check_arg_within_options(x, c(1,2)),
                                      "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)") &
        .check_argument_correct(filtering, .check_boolean, "Check filtering - should be either TRUE or FALSE") &
        .check_argument_correct(reducedDim_name, is.character, "Check reducedDim_name - should be character vector") &
        .check_argument_correct(k_init, .check_positive_integer, "Check k_init - should be positive integer") &
        .check_argument_correct(d, .check_positive_integer, "Check d - should be positive integer") &
        .check_argument_correct(verbose, .check_boolean, "Check verbose - should be either TRUE or FALSE") &
        .check_reducedDim_in_sce(x , reducedDim_name)


  d <- min(d , ncol(reducedDim(x , reducedDim_name)))
  k_init <- min(k , k_init)

  if (is.null(colnames(x))){
    colnames(x) = as.character(c(1:ncol(x)))
  }
  if (is(x , "SingleCellExperiment")){
    x = Milo(x)
    # build 1st order to sample vertices
    k_init <- min(k , k_init)
    x <- suppressMessages(buildGraph(x, k = k_init, d = d, reduced.dim = reducedDim_name))
  }
  else {
    message("SCE is Milo object. Checking if graph is already constructed.")
    if (length(miloR::graph(x)) == 0){
      message("Graph is not constructed yet. Building now.")
      x <- suppressMessages(buildGraph(x, k = k_init, d = d, reduced.dim = reducedDim_name))
    }
  }
  # find anchor cells
  sampled_vertices <- .get_graph_refined_sampling(graph(x), prop)

  # rebuild to the actual graph, with parameters specified by user
  if (!k == k_init){
    x <- suppressMessages(buildGraph(x, k = k, d = d, reduced.dim = reducedDim_name))
  }
  # if order == 2 -- reassign edges
  if (order == 2){
    graph(x) = connect(graph(x),order)
  }

  # create nhoods
  nh_mat <- Matrix(data = 0, nrow=ncol(x), ncol=length(sampled_vertices), sparse = TRUE)
  v.class <- V(graph(x))$name
  rownames(nh_mat) <- colnames(x)
  for (X in seq_len(length(sampled_vertices))){
    nh_mat[unlist(neighborhood(graph(x), order = 1, nodes = sampled_vertices[X])), X] <- 1
  }
  colnames(nh_mat) <- as.character(sampled_vertices)
  nhoodIndex(x) <- as(sampled_vertices, "list")
  nhoods(x) <- nh_mat

  # filter
  if (!filtering){
    x = suppressMessages(buildNhoodGraph(x))
  }
  else {
    if (verbose){
      message("Filtering redundant neighbourhoods.")
    }
    x = suppressMessages(filter_neighbourhoods(x))
  }

  if (verbose){
    stat_print =.calc_quick_stat(x , nhoods(x))
    message(paste0("Finished successfully.\nNumber of neighbourhoods assigned: ", stat_print$n_hoods ,
                 ";\naverage neighbourhood size: ", stat_print$avg_hood_size ,
                 ";\nnumber of unassigned cells: ", stat_print$n_cells_unocovered))
  }
  return(x)
}


#'
#'
#' @importFrom igraph V set_vertex_attr induced_subgraph count_triangles neighborhood
#' @importFrom miloR graph nhoodIndex nhoods<-
#' @importFrom dplyr %>%
.get_graph_refined_sampling <- function(graph, prop){
  random_vertices <- sample(V(graph), size=floor(prop*length(V(graph))))
  random_vertices <- as.vector(random_vertices)
  X_graph <- set_vertex_attr(graph, "name", value = 1:length(V(graph)))
  refined_vertices <- lapply(seq_along(random_vertices), function(i){
    target_vertices <- unlist(neighborhood(X_graph, order = 1, nodes = random_vertices[i])) #get neighborhood of random vertex
    target_vertices <- target_vertices[-1] #remove first entry which is the random vertex itself
    rv_induced_subgraph <- induced_subgraph(graph = X_graph, vids = target_vertices)
    triangles <- count_triangles(rv_induced_subgraph)
    max_triangles <- max(triangles)
    max_triangles_indices <- which(triangles == max_triangles)
    #note - take first max_ego_index in the next line of code
    resulting_vertices <- V(rv_induced_subgraph)[max_triangles_indices]$name[1]
    return(resulting_vertices)
  }) %>% unlist() %>% as.integer()
  refined_vertices = unique(refined_vertices)
  return(refined_vertices)
}

#'
.calc_quick_stat = function(x , nhoods_sce){
  n_hoods = ncol(nhoods_sce)
  avg_hood_size = round(mean(colSums(nhoods_sce)))
  n_cells_unocovered = ncol(x) - sum(rowSums(nhoods_sce) > 0)
  out = list(n_hoods = n_hoods ,
             avg_hood_size = avg_hood_size ,
             n_cells_unocovered = n_cells_unocovered)
  return(out)
}
