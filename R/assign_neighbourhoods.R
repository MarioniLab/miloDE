

#' assign_neighbourhoods
#'
#' Assign hoods to SCE object
#' @param sce SCE object
#' @param k Positive integer, defines how many neighbours to use for the hood assignment
#' @param prop Numerical, between 0 and 1, defines fraction of cells from SCE to use for the hoods
#' @param order In {1,2}, defines which order of neighbours to use
#' @param filtering In {T,F}, defines whether to filter hoods (reduces computing time greatly). Default = TRUE
#' @param reducedDim.name defines the slot in reducedDim(sce) to use as embedding for graph construction
#' @param k_init Positive integer, defines how many neighbours to use for identifying anchor cells
#' @param d Positive integer, defines how many dimensions from reducedDim(sce) to use
#'
#' @return
#' @export
#' @importFrom miloR buildGraph Milo graph
#' @importFrom SingleCellExperiment reducedDim
#' @import Matrix
#' @import igraph
#' @examples
assign_neighbourhoods = function(sce , k = 25, prop = 0.2, order = 2, filtering = T, reducedDim.name , k_init = 50, d = 30){
  sce_milo <- Milo(sce)
  d <- min(30 , ncol(reducedDim(sce , reducedDim.name)))

  # build 1st order to sample vertices
  k_init <- min(k , k_init)
  sce_milo <- buildGraph(sce_milo, k = k_init, d = d, reduced.dim = reducedDim.name)

  # find anchor cells
  sampled_vertices <- .get_graph_refined_sampling(miloR::graph(sce_milo), prop)

  # rebuild to the actual graph, with parameters specified by user
  sce_milo <- buildGraph(sce_milo, k = k, d = d, reduced.dim = reducedDim.name)
  # if order > 1 -- reassign
  if (order > 1){
    graph(sce_milo) = connect(miloR::graph(sce_milo),order)
  }

  # create hoods
  nh_mat <- Matrix(data = 0, nrow=ncol(sce_milo), ncol=length(sampled_vertices), sparse = TRUE)
  v.class <- V(miloR::graph(sce_milo))$name
  rownames(nh_mat) <- colnames(sce_milo)
  for (X in seq_len(length(sampled_vertices))){
    nh_mat[unlist(neighborhood(miloR::graph(sce_milo), order = 1, nodes = sampled_vertices[X])), X] <- 1
  }
  colnames(nh_mat) <- as.character(sampled_vertices)
  nhoodIndex(sce_milo) <- as(sampled_vertices, "list")
  nhoods(sce_milo) <- nh_mat

  # filter
  if (!filtering){
    sce_milo = buildNhoodGraph(sce_milo)
  }
  else {
    sce_milo = filter_neighbourhoods(sce_milo)
  }
  return(sce_milo)
}


#'
#'
#' @import igraph
#' @importFrom miloR graph
.get_graph_refined_sampling <- function(graph, prop){
  require(igraph)
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





