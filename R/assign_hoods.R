

#' assign_hoods
#'
#' Assign hoods to SCE object
#' @param sce SCE object
#' @param k Positive integer, defines how many neighbours to use for the hood assignment
#' @param prop Numerical, between 0 and 1, defines fraction of cells from SCE to use for the hoods
#' @param order In {1,2}, defines which order of neighbours to use
#' @param filtering In {TRUE,FALSE}, defines whether to filter hoods (reduces computing time greatly). Default = TRUE
#' @param reducedDim_name defines the slot in reducedDim(sce) to use as embedding for graph construction
#' @param k_init Positive integer, defines how many neighbours to use for identifying anchor cells
#' @param d Positive integer, defines how many dimensions from reducedDim(sce) to use
#'
#' @return
#' @export
#' @importFrom miloR Milo buildGraph graph<- graph nhoods<- nhoodIndex<- buildNhoodGraph
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assay
#' @import Matrix
#' @importFrom igraph connect V neighborhood
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
#' out = assign_hoods(sce, reducedDim_name = "reduced_dim")
assign_hoods = function(sce , k = 25, prop = 0.2, order = 2, filtering = TRUE, reducedDim_name , k_init = 50, d = 30){

  args = c(as.list(environment()))
  out = .general_check_arguments(args) & .check_reducedDim_in_sce(sce , reducedDim_name)

  d <- min(d , ncol(reducedDim(sce , reducedDim_name)))
  k_init <- min(k , k_init)
  if (is(sce , "SingleCellExperiment")){
    sce = Milo(sce)
    # build 1st order to sample vertices
    k_init <- min(k , k_init)
    sce <- buildGraph(sce, k = k_init, d = d, reduced.dim = reducedDim_name)
  }
  else {
    message("SCE is Milo object. Checking if graph is already constructed.")
    if (length(miloR::graph(sce)) == 0){
      message("Graph is not constructed yet. Building now.")
      sce <- buildGraph(sce, k = k_init, d = d, reduced.dim = reducedDim_name)
    }
  }
  # find anchor cells
  sampled_vertices <- .get_graph_refined_sampling(graph(sce), prop)

  # rebuild to the actual graph, with parameters specified by user
  if (!k == k_init){
    sce <- buildGraph(sce, k = k, d = d, reduced.dim = reducedDim_name)
  }
  # if order == 2 -- reassign edges
  if (order == 2){
    graph(sce) = connect(graph(sce),order)
  }

  # create hoods
  nh_mat <- Matrix(data = 0, nrow=ncol(sce), ncol=length(sampled_vertices), sparse = TRUE)
  v.class <- V(graph(sce))$name
  rownames(nh_mat) <- colnames(sce)
  for (X in seq_len(length(sampled_vertices))){
    nh_mat[unlist(neighborhood(graph(sce), order = 1, nodes = sampled_vertices[X])), X] <- 1
  }
  colnames(nh_mat) <- as.character(sampled_vertices)
  nhoodIndex(sce) <- as(sampled_vertices, "list")
  nhoods(sce) <- nh_mat

  # filter
  if (!filtering){
    sce = buildNhoodGraph(sce)
  }
  else {
    message("Filtering redundant hoods.")
    sce = suppressMessages(filter_hoods(sce))
  }

  stat_print =.calc_quick_stat(sce , nhoods(sce))
  message(paste0("Finished successfully.\nNumber of hoods assigned: ", stat_print$n_hoods ,
                 ",\naverage hood size: ", stat_print$avg_hood_size ,
                 ",\nnumber of unassigned cells: ", stat_print$n_cells_unocovered))
  return(sce)
}


#'
#'
#' @importFrom igraph V set_vertex_attr induced_subgraph count_triangles neighborhood
#' @importFrom miloR graph nhoodIndex nhoods<-
#' @importFrom dplyr %>%
.get_graph_refined_sampling <- function(graph, prop){
  random_vertices <- sample(V(graph), size=floor(prop*length(V(graph))))
  message("Running refined sampling with graph")
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
.calc_quick_stat = function(sce , nhoods_sce){
  n_hoods = ncol(nhoods_sce)
  avg_hood_size = round(mean(colSums(nhoods_sce)))
  n_cells_unocovered = ncol(sce) - sum(rowSums(nhoods_sce) > 0)
  out = list(n_hoods = n_hoods ,
             avg_hood_size = avg_hood_size ,
             n_cells_unocovered = n_cells_unocovered)
  return(out)
}




