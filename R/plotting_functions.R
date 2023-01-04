
# plotting functions



#' plot_milo_by_single_metric
#'
#' Return Milo-graph plot; each node is coloured by colour_by column from nhood_stat, if significance_by is smaller than alpha. Vetices are ordered by order_by column
#' @param sce \code{\linkS4class{Milo}} object
#' @param nhood_stat data.frame object, containing columns 'Nhood' (should correspond to neighbourhoods from nhoodGraph(sce))
#' @param colour_by A colname from nhood_stat - nodes will be coloured by the values from this column
#' @param significance_by A colname from nhood_stat (or NULL) - if values for this column exceed alpha, colour_by will be set to 0. Default = NULL meaning that we will not use any column to define the significance (i.e. no correction for colour_by)
#' @param order_by A colname from nhood_stat (or NULL) specifying by which column we will order neighbourhoods for plotting. Default = NULL meaning we will use size of the neighbourhoods.
#' @param order_direction Boolean specifying the direction of ordering neighbourhoods
#' @param size_by NULL or a character vector - colname from nhood_stat - specifying metric by which we define size for neighbourhoods
#' @param alpha A scalar (between 0 and 1) specifying the significance level used
#' @param layout A character indicating the name of the \code{reducedDim} slot in the \code{\linkS4class{Milo}} object to use for layout (default: 'UMAP')
#' @param subset_nhoods A vector (or NULL) specifying which neighbourhoods will be plotted. Default = NULL meaning that all neighbourhoods will be plotted
#' @param size_range A numeric vector indicating the range (min and max) of node sizes to use for plotting (to avoid overplotting in the graph)
#' @param node_stroke A numeric indicating the desired thickness of the border around each node
#' @param edge_width A numeric vector indicating the range (min and max) of edge widths to use for plotting
#' @param edge_weight.thresh A numeric (or NULL) specifying a threshold for minimum cells in common (between neighbourhoods) required for an edge to be plotted
#'
#' @return
#' @importFrom SingleCellExperiment reducedDims reducedDim
#' @importFrom SummarizedExperiment colData<- colData
#' @importFrom miloR nhoodIndex nhoodIndex<- nhoodGraph nhoodGraph<-
#' @importFrom igraph induced_subgraph vertex_attr vertex_attr<- permute delete.edges simplify V<- V E<- E
#' @importFrom ggraph ggraph geom_edge_link geom_edge_link0 geom_node_point scale_edge_width
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
#' @importFrom methods as
#' @export
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
#' sce$sample = floor(runif(n = n_col , min = 1 , max = 5))
#' sce$type = ifelse(sce$sample %in% c(1,2) , "ref" , "query")
#' reducedDim(sce , "reduced_dim") = matrix(rnorm(n_col*n_latent), ncol=n_latent)
#' sce = assign_neighbourhoods(sce, reducedDim_name = "reduced_dim")
#' de_stat = de_test_neighbourhoods(sce , condition_id = "type" , gene_selection = "none")
#' de_stat = de_stat[de_stat$gene == "1", ]
#' umaps = as.data.frame(matrix(rnorm(n_col*2), ncol=2))
#' colnames(umaps) = c("V1" , "V2")
#' reducedDim(sce , "UMAP") = umaps
#' p = plot_milo_by_single_metric(sce, de_stat)
#'
plot_milo_by_single_metric = function(sce, nhood_stat, colour_by = "logFC" , significance_by = NULL , order_by = NULL , order_direction = TRUE,
                                      size_by = NULL ,
                                      alpha = 0.1, layout = "UMAP" , subset_nhoods = NULL , size_range = c(1,3) ,
                                      node_stroke = 0.3, edge_width = c(0.2,0.5), edge_weight.thresh = NULL){

  # checks
  out = .check_argument_correct(sce, .check_sce, "Check sce - something is wrong (gene names unique? reducedDim.name is not present?)") &
    .check_sce_milo(sce) &
    .check_argument_correct(order_direction, .check_boolean, "Check order_direction - should be either TRUE or FALSE") &
    .check_reducedDim_in_sce(sce , layout) & .check_nhood_stat(nhood_stat , sce)

  if (!colour_by %in% colnames(nhood_stat)){
    stop("colour_by should be in colnames(nhood_stat)")
  }
  if (!is.null(significance_by)){
    if (!significance_by %in% colnames(nhood_stat)){
      stop("significance_by should be NULL or in colnames(nhood_stat)")
    }
  }
  if (!is.null(order_by)){
    if (!order_by %in% colnames(nhood_stat)){
      stop("order_by should be NULL or in colnames(nhood_stat)")
    }
  }
  if (!is.null(size_by)){
    if (!size_by %in% colnames(nhood_stat)){
      stop("size_by should be NULL or in colnames(nhood_stat)")
    }
  }

  if (!is.null(subset_nhoods)){
    nhoods_sce = nhoods(sce)
    out = .check_subset_nhoods(subset_nhoods , nhoods_sce)
  }


  if(!.valid_graph(nhoodGraph(sce))){
    stop("Not a valid Milo object - neighbourhood graph is missing. Please run buildNhoodGraph() first.")
  }
  if (is.character(layout)) {
    if (!layout %in% names(reducedDims(sce))) {
      stop(layout, "isn't in readucedDim(sce) - choose a different layout")
    }
  }


  nhood_stat = nhood_stat[order(nhood_stat$Nhood) , ]
  # for hoods that exceed alpha for significance_by -- set colour_by to 0
  if (!is.null(significance_by)){
    nhood_stat[nhood_stat[, significance_by] > alpha, colour_by] <- 0
  }

  # assign colour_by to sce
  colData(sce)[colour_by] <- NA
  colData(sce)[unlist(nhoodIndex(sce)[nhood_stat$Nhood]),colour_by] <- nhood_stat[,colour_by]

  # pull the graph from sce
  nh_graph <- nhoodGraph(sce)

  ## subset for hoods we will plot
  if (!is.null(subset_nhoods)) {
    nh_graph <- induced_subgraph(nh_graph, vids = which(as.numeric(V(nh_graph)$name) %in% unlist(nhoodIndex(sce)[subset_nhoods])))
  }

  if (!is.null(size_by)){
    vertex_attr(nh_graph)$size = rep(0 , 1 , length(vertex_attr(nh_graph)$name) )
    vertex_attr(nh_graph)$size[nhood_stat$Nhood] = nhood_stat[,size_by]
  }


  # assign attributes to vertices
  #vertex_attr(nh_graph)[, significance_by] = nhood_stat[, significance_by]
  if (!is.null(order_by)){
    vertex_attr(nh_graph)$order_by = rep(NA , 1 , length(vertex_attr(nh_graph)$name) )
    vertex_attr(nh_graph)$order_by[nhood_stat$Nhood] = nhood_stat[,order_by]
    vertex_attr(nh_graph)$order_by_rearrange = order(vertex_attr(nh_graph)$order_by , decreasing = order_direction, na.last = FALSE)
  }
  else {
    vertex_attr(nh_graph)$order_by_rearrange = order(vertex_attr(nh_graph)$size , decreasing = FALSE , na.last = FALSE)
  }
  nh_graph <- permute(nh_graph, match( 1:length(vertex_attr(nh_graph)$order_by_rearrange) ,
                                                 vertex_attr(nh_graph)$order_by_rearrange))


  # assign edges lower than some thresh to 0
  if (!is.null(edge_weight.thresh)){
    nh_graph <- delete.edges(nh_graph, which(E(nh_graph)$weight <= edge_weight.thresh)-1)
  }

  ## define layout
  if (is.character(layout)) {
    redDim <- layout
    layout <- reducedDim(sce, redDim)[as.numeric(vertex_attr(nh_graph)$name),]
    # make sure this is a matrix!
    if(!any(class(layout) %in% c("matrix"))){
      warning("Coercing layout to matrix format")
      layout <- as(layout, "matrix")
    }
  }


  ## Define node color
  if (colour_by %in% colnames(colData(sce))) {

    col_vals <- colData(sce)[as.numeric(vertex_attr(nh_graph)$name), colour_by]
    if (!is.numeric(col_vals)) {
      col_vals <- as.character(col_vals)
    }
    V(nh_graph)$colour_by <- col_vals
  } else {
    stop(colour_by, "is not a column in colData(x)")
  }



  pl <- ggraph(simplify(nh_graph), layout = layout) +
    geom_edge_link0(aes(width = weight), edge_colour = "grey66", edge_alpha=0.2) +
    geom_node_point(aes(fill = colour_by, size = size), shape=21, stroke=node_stroke) +
    scale_size(range = size_range, name="Nhood size") +
    scale_edge_width(range = c(0.2,3), name="overlap size") +
    theme_classic(base_size=14) +
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank()) +
    guides(width="none" ,edge_width="none")

  if (is.numeric(V(nh_graph)$colour_by)) {
    pl <- pl + scale_fill_gradient2(name=colour_by)
  } else {
    mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(V(nh_graph)$colour_by)))
    pl <- pl + scale_fill_manual(values=mycolors, name=colour_by, na.value="white")
  }
  return(pl)
}




#
#' plot_DE_single_gene
#'
#' Return Milo-graph plot; each node is coloured by logFC, if p-value corrected across nhoods < alpha
#' @param sce Milo object
#' @param de_stat milo-DE stat (output of 'de_stat_all_neighbourhoods')
#' @param gene A character specifying gene ID
#' @param alpha A scalar (between 0 and 1) specifying the significance level used
#' @param layout A character indicating the name of the \code{reducedDim} slot in the \code{\linkS4class{Milo}} object to use for layout (default: 'UMAP')
#' @param subset_nhoods A vector (or NULL) specifying which neighbourhoods will be plotted. Default = NULL meaning that all neighbourhoods will be plotted
#' @param ... Arguments to pass to \code{plot_milo_by_single_metric} (e.g. size_range, node_stroke etc))
#'
#' @return
#' @export
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
#' sce$sample = floor(runif(n = n_col , min = 1 , max = 5))
#' sce$type = ifelse(sce$sample %in% c(1,2) , "ref" , "query")
#' reducedDim(sce , "reduced_dim") = matrix(rnorm(n_col*n_latent), ncol=n_latent)
#' sce = assign_neighbourhoods(sce, reducedDim_name = "reduced_dim")
#' de_stat = de_test_neighbourhoods(sce , condition_id = "type" , gene_selection = "none")
#' umaps = as.data.frame(matrix(rnorm(n_col*2), ncol=2))
#' colnames(umaps) = c("V1" , "V2")
#' reducedDim(sce , "UMAP") = umaps
#' p = plot_DE_single_gene(sce, de_stat , gene = "1")
#'
plot_DE_single_gene = function(sce, de_stat , gene , alpha = 0.1, layout = "UMAP" , subset_nhoods = NULL , ...){

  if (!gene %in% rownames(sce)){
    stop("gene should be in rownames(sce)")
  }

  out = .check_de_stat_valid(de_stat ,
                             assay_names = c("logFC" , "pval" , "pval_corrected_across_nhoods" , "pval_corrected_across_genes") ,
                             coldata_names = c("Nhood" , "Nhood_center"))

  if (class(de_stat) == "SingleCellExperiment"){
    de_stat = de_stat[gene , ]
    de_stat = convert_de_stat(de_stat, assay_names = c("logFC" , "pval" , "pval_corrected_across_nhoods" , "pval_corrected_across_genes") ,
                              coldata_names = c("Nhood" , "Nhood_center" , "sufficient_n_samples" , "design_matrix_suitable"))
  }

  nhood_stat = de_stat[de_stat$gene == gene & de_stat$design_matrix_suitable == TRUE, ]

  p = plot_milo_by_single_metric(sce, nhood_stat = nhood_stat, colour_by = "logFC" , significance_by = "pval_corrected_across_nhoods" ,
                             order_by = "pval_corrected_across_nhoods" , order_direction = TRUE, size_by = NULL,
                             alpha = alpha, layout = layout , subset_nhoods = subset_nhoods  , ...)
  return(p)

}





#' plot_DE_gene_set
#'
#' Returns Milo plot, in which colour of nodes correposnd to average logFC across selected genes; size corresponds to how many genes show significant DE in the neighbourhood (based on pval_corrected_across_nhoods)
#' @param sce Milo object
#' @param de_stat milo-DE stat (output of 'de_stat_neighbourhoods')
#' @param genes Character vector, each element corresponds to gene ID
#' @param logFC_correction Boolean specifying whether to perfrom logFC correction. If TRUE (default), logFC will be set to 0 if corrected pvalue (defined by 'correction_by' variable) < alpha
#' @param correction_by Character - colname from de_stat - specifying which column to use to decide on significance for logFC. Relevant only if logFC_correction = TRUE.
#' @param alpha A numeric between 0 and 1 specifying the threshold for correction
#' @param layout A character indicating the name of the \code{reducedDim} slot in the \code{\linkS4class{Milo}} object to use for layout (default: 'UMAP')
#' @param subset_nhoods A vector (or NULL) specifying which neighbourhoods will be plotted. Default = NULL meaning that all neighbourhoods will be plotted
#' @param ... Arguments to pass to \code{plot_milo_by_single_metric} (e.g. size_range, node_stroke etc))
#'
#' @return
#' @importFrom SummarizedExperiment assay assay<- colData
#' @export
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
#' sce$sample = floor(runif(n = n_col , min = 1 , max = 5))
#' sce$type = ifelse(sce$sample %in% c(1,2) , "ref" , "query")
#' reducedDim(sce , "reduced_dim") = matrix(rnorm(n_col*n_latent), ncol=n_latent)
#' sce = assign_neighbourhoods(sce, reducedDim_name = "reduced_dim")
#' de_stat = de_test_neighbourhoods(sce , condition_id = "type" , gene_selection = "none")
#' umaps = as.data.frame(matrix(rnorm(n_col*2), ncol=2))
#' colnames(umaps) = c("V1" , "V2")
#' reducedDim(sce , "UMAP") = umaps
#' genes = c("1","2")
#' p = plot_DE_gene_set(sce, de_stat , genes = c("1","2"))
#'
plot_DE_gene_set = function(sce, de_stat , genes ,
                            logFC_correction = TRUE , correction_by = "pval_corrected_across_nhoods", alpha = 0.1,
                            layout = "UMAP" , subset_nhoods = NULL , ...){
  # checks
  if (mean(genes %in% rownames(sce)) < 1){
    stop("genes should be in rownames(sce)")
  }
  out = .check_de_stat_valid(de_stat ,
                             assay_names = c("logFC" , "pval" , "pval_corrected_across_genes" , "pval_corrected_across_nhoods") ,
                             coldata_names = c("Nhood" , "Nhood_center"))

  if (class(de_stat) == "data.frame"){
    if (mean(genes %in% unique(de_stat$gene)) < 1){
      stop("All genes should be in de_stat$gene")
    }
    else {
      de_stat = de_stat[de_stat$gene %in% genes , ]
      de_stat = convert_de_stat(de_stat ,
                              assay_names = c("logFC" , "pval"  , "pval_corrected_across_genes" , "pval_corrected_across_nhoods") ,
                              coldata_names = c("Nhood" , "Nhood_center"))
    }
  }
  else {
    if (mean(genes %in% rownames(de_stat)) < 1){
      stop("All genes should be in rownames(de_stat)")
    }
    else {
      de_stat = de_stat[genes, ]
    }
  }

  if (logFC_correction){
    assay_correction_by = assay(de_stat , correction_by)
    idx_sig = which(assay_correction_by < alpha)
    idx_not_sig = which(assay_correction_by >= alpha)
    assay_correction_by[idx_sig] = 1
    assay_correction_by[idx_not_sig] = 0

    assay_logFC = assay(de_stat , "logFC")
    assay(de_stat , "logFC_corrected") = assay_logFC * assay_correction_by
  }
  else {
    assay(de_stat , "logFC_corrected") = assay(de_stat , "logFC")
  }


  # get stat - average logFC and fraction of genes for which this neighbourhood is significant
  nhood_stat = as.data.frame(colData(de_stat))
  nhood_stat$avg_logFC = colMeans(assay(de_stat , "logFC_corrected"))
  nhood_stat$frac_DE_genes = colMeans(assay(de_stat , "pval_corrected_across_nhoods") < alpha)

  p = plot_milo_by_single_metric(sce, nhood_stat, colour_by = "avg_logFC" , significance_by = NULL ,
                                 order_by = "frac_DE_genes" , order_direction = FALSE, size_by = "frac_DE_genes",
                                           layout = layout , subset_nhoods = subset_nhoods , ...)
  return(p)
}






#' plot_beeswarm_single_gene
#'
#' Returns beeswarm plot for a single gene
#' @param de_stat milo-DE stat (output of 'de_stat_neighbourhoods')
#' @param gene A character specifying gene ID
#' @param nhoodGroup A character specifying which column to use for neighbourhood grouping.
#' @param alpha A numeric between 0 and1 specifying the significance threshold. Default = 0.1.
#' @param subset_nhoods NULL or numeric vector specifying which Nhoods to use
#'
#' @return
#' @importFrom dplyr mutate %>% arrange
#' @importFrom ggbeeswarm geom_quasirandom
#' @import ggplot2
#' @export
#'
#' @examples
#'
#' require(SingleCellExperiment)
#' n_row = 500
#' n_col = 100
#' n_latent = 5
#' sce = SingleCellExperiment(assays = list(counts = floor(matrix(rnorm(n_row*n_col), ncol=n_col)) + 4))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' sce$sample = floor(runif(n = n_col , min = 1 , max = 5))
#' sce$type = ifelse(sce$sample %in% c(1,2) , "ref" , "query")
#' reducedDim(sce , "reduced_dim") = matrix(rnorm(n_col*n_latent), ncol=n_latent)
#' sce = assign_neighbourhoods(sce, reducedDim_name = "reduced_dim")
#' de_stat = de_test_neighbourhoods(sce , condition_id = "type" , gene_selection = "none")
#' de_stat$celltype = 1
#' p = plot_beeswarm_single_gene(de_stat , gene = "1" , nhoodGroup = "celltype")
#'
plot_beeswarm_single_gene = function(de_stat , gene , nhoodGroup , alpha = 0.1 , subset_nhoods = NULL ){

  out = .check_de_stat_valid(de_stat , assay_names = NULL, coldata_names = nhoodGroup)

  if (class(de_stat) == "SingleCellExperiment"){
    if (!gene %in% rownames(de_stat)){
      stop("gene should be in rownames(de_stat)")
    }
    else {
      de_stat = de_stat[gene , ]
      de_stat = convert_de_stat(de_stat , assay_names = NULL, coldata_names = nhoodGroup)
    }
  } else {
    if (!gene %in% unique(de_stat$gene)){
      stop("gene should be in de_stat$gene")
    }
    else {
      de_stat = de_stat[de_stat$gene == gene , ]
    }
  }

  nhood_stat = de_stat
  nhood_stat = nhood_stat[order(nhood_stat$Nhood) , ]
  nhood_stat[, nhoodGroup] = as.factor(nhood_stat[, nhoodGroup])
  nhood_stat = mutate(nhood_stat, group_by = nhood_stat[,nhoodGroup])

  if (!is.null(subset_nhoods)) {
    nhood_stat <- nhood_stat[nhood_stat$Nhood %in% subset_nhoods,]
  }

  p = nhood_stat %>%
    mutate(is_signif = ifelse(pval_corrected_across_nhoods < alpha, 1, 0)) %>%
    mutate(logFC_color = ifelse(is_signif==1, logFC, NA)) %>%
    arrange(group_by) %>%
    mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
    ggplot(aes(group_by, logFC, color=logFC_color)) +
    scale_color_gradient2() +
    guides(color="none") +
    xlab(nhoodGroup) + ylab("Log Fold Change") +
    geom_quasirandom(alpha=1) +
    coord_flip() +
    theme_bw(base_size=22) +
    theme(strip.text.y =  element_text(angle=0))

  return(p)

}



#' plot_beeswarm_gene_set
#'
#' Returns beeswarm plot for many genes
#' @param de_stat milo-DE stat (output of 'de_stat_neighbourhoods')
#' @param genes A character specifying genes ID
#' @param nhoodGroup A character specifying which column to use for neighbourhood grouping.
#' @param logFC_correction Boolean specifying whether to perfrom logFC correction. If TRUE (default), logFC will be set to 0 if corrected pvalue (defined by 'correction_by' variable) < alpha
#' @param correction_by Character - colname from de_stat - specifying which column to use to decide on significance for logFC. Relevant only if logFC_correction = TRUE.
#' @param alpha A numeric between 0 and1 specifying the significance threshold. Default = 0.1.
#' @param subset_nhoods NULL or numeric vector specifying which Nhoods to use
#'
#' @return
#' @importFrom dplyr mutate %>% arrange
#' @importFrom ggbeeswarm geom_quasirandom
#' @import ggplot2
#' @importFrom SummarizedExperiment assay assay<-
#' @export
#'
#' @examples
#'
#' require(SingleCellExperiment)
#' n_row = 500
#' n_col = 100
#' n_latent = 5
#' sce = SingleCellExperiment(assays = list(counts = floor(matrix(rnorm(n_row*n_col), ncol=n_col)) + 4))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' sce$sample = floor(runif(n = n_col , min = 1 , max = 5))
#' sce$type = ifelse(sce$sample %in% c(1,2) , "ref" , "query")
#' reducedDim(sce , "reduced_dim") = matrix(rnorm(n_col*n_latent), ncol=n_latent)
#' sce = assign_neighbourhoods(sce, reducedDim_name = "reduced_dim")
#' de_stat = de_test_neighbourhoods(sce , condition_id = "type" , gene_selection = "none")
#' de_stat$celltype = 1
#' p = plot_beeswarm_gene_set(de_stat , genes = c("1","2") , nhoodGroup = "celltype")
#'
plot_beeswarm_gene_set = function(de_stat , genes , nhoodGroup , logFC_correction = TRUE ,
                                  correction_by = "pval_corrected_across_nhoods" , alpha = 0.1 , subset_nhoods = NULL ){

  out = .check_de_stat_valid(de_stat , assay_names = NULL, coldata_names = nhoodGroup)


  if (class(de_stat) == "data.frame"){
    if (mean(genes %in% unique(de_stat$gene)) < 1){
      stop("genes should be in de_stat$gene")
    }
    else {
      de_stat = convert_de_stat(de_stat , assay_names = NULL, coldata_names = nhoodGroup)
    }
  } else {
    if (mean(genes %in% rownames(de_stat)) < 1){
      stop("genes should be in rownames(de_stat)")
    }
  }
  de_stat = de_stat[genes , ]

  if (logFC_correction){
    assay_correction_by = assay(de_stat , correction_by)
    idx_sig = which(assay_correction_by < alpha)
    idx_not_sig = which(assay_correction_by >= alpha)
    assay_correction_by[idx_sig] = 1
    assay_correction_by[idx_not_sig] = 0

    assay_logFC = assay(de_stat , "logFC")
    assay(de_stat , "logFC_corrected") = assay_logFC * assay_correction_by
  }
  else {
    assay(de_stat , "logFC_corrected") = assay(de_stat , "logFC")
  }

  nhood_stat = as.data.frame(colData(de_stat))
  nhood_stat = nhood_stat[order(nhood_stat$Nhood) , ]
  nhood_stat$avg_logFC = colMeans(assay(de_stat , "logFC_corrected"))
  nhood_stat$frac_DE_genes = colMeans(assay(de_stat , "pval_corrected_across_nhoods") < alpha)
  nhood_stat = mutate(nhood_stat, group_by = nhood_stat[,nhoodGroup])

  if (!is.null(subset_nhoods)) {
    nhood_stat <- nhood_stat[nhood_stat$Nhood %in% subset_nhoods,]
  }

  p = nhood_stat %>%
    arrange(group_by) %>%
    mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
    ggplot(aes(group_by, frac_DE_genes, color=avg_logFC)) +
    scale_color_gradient2(name = "Average logFC") +
    guides(color="none") +
    xlab(nhoodGroup) + ylab("How often gene is DE") +
    geom_quasirandom(alpha=1) +
    coord_flip() +
    theme_bw(base_size=22) +
    theme(strip.text.y =  element_text(angle=0))

  return(p)

}



