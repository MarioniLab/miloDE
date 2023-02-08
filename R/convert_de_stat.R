

#' convert_de_stat
#'
#' Converts output of miloDE between \code{\link[base]{data.frame}} and \code{\linkS4class{SingleCellExperiment}} formats
#' @param de_stat miloDE results, output of \code{\link{de_test_neighbourhoods}}; either in \code{data.frame} or \code{SingleCellExperiment}.
#' @param assay_names Character string specifying which fields are treated as assays.
#' Note that \code{logFC}, \code{pval}, \code{pval_corrected_across_genes} and \code{pval_corrected_across_nhoods} are hard-coded to be included.
#' @param coldata_names Character string specifying which fields are treated as \code{Nhood} metadata.
#' Note that \code{Nhood}, \code{Nhood_center}, \code{test_performed} are hard-coded to be included.
#'
#' \emph{Please note that \code{coldata_names} have to be the attributes of neighbourhoods (i.e. same across different genes for the same neighbourhood).}
#' @details
#' This function converts results of \code{\link{de_test_neighbourhoods}} between \code{data.frame} object and \code{SingleCellExperiment}.
#' \code{data.frame} object is more commonly used and might be easier to navigate, however, if total number of tests (i.e. gene x neighboourhoods)
#' is overwhelmingly large, \code{SingleCellExperiment} might be more suitable.
#' @return A \code{\linkS4class{SingleCellExperiment}} object or \code{data.frame} object, containing miloDE results.
#' @export
#' @examples
#' de_stat = expand.grid(gene = paste0("gene" , c(1:5)) , Nhood = c(1:10))
#' de_stat$Nhood_center = paste0("nhood_" , de_stat$Nhood)
#' de_stat$logFC = sample(seq(-2,2,1) , nrow(de_stat) , 1)
#' de_stat$pval = sample(c(0,1),nrow(de_stat),1)
#' de_stat$pval_corrected_across_genes = sample(c(0,1),nrow(de_stat),1)
#' de_stat$pval_corrected_across_nhoods = sample(c(0,1),nrow(de_stat),1)
#' de_stat$test_performed = TRUE
#' de_stat = convert_de_stat(de_stat)
#' de_stat = convert_de_stat(de_stat)
#'
convert_de_stat = function(de_stat ,
                           assay_names = NULL,
                           coldata_names = NULL){

  assay_names = unique( c("logFC" , "pval" , "pval_corrected_across_genes" , "pval_corrected_across_nhoods", assay_names))
  coldata_names = unique( c("Nhood" , "Nhood_center" , "test_performed"  , coldata_names))
  out = .check_de_stat_valid(de_stat , assay_names , coldata_names)

  if (is(de_stat , "SingleCellExperiment")){
    message("Converting de_stat to 'data.frame' format")
    de_stat = .convert_from_sce(de_stat , assay_names = assay_names , coldata_names = coldata_names)
  } else if (is(de_stat , "data.frame")){
    message("Converting de_stat to 'SingleCellExperiment' format")
    de_stat = .convert_from_df(de_stat , assay_names = assay_names , coldata_names = coldata_names)
  }
  return(de_stat)
}


#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
.convert_from_df = function(de_stat , assay_names , coldata_names){

  de_stat = de_stat[order(de_stat$Nhood) , ]
  #  convert assays
  de_assays = lapply(assay_names , function(assay_name){
    return(.convert_from_df_one_var(de_stat , assay_name))
  })
  names(de_assays) = assay_names

  # convert coldata
  meta_nhoods = unique(de_stat[, coldata_names])
  meta_nhoods = meta_nhoods[order(meta_nhoods$Nhood) , ]

  # combine
  de_stat = SingleCellExperiment(de_assays, colData = DataFrame(meta_nhoods))

  colnames(de_stat) = meta_nhoods$Nhood
  return(de_stat)
}



#' @importFrom reshape2 dcast
.convert_from_df_one_var = function(de_stat , var){
  df = dcast(de_stat , formula = gene ~ Nhood , value.var = var)
  df$gene = as.character(df$gene)
  n_hoods = ncol(df)-1
  n_genes = nrow(df)
  genes = df$gene
  df = df[, 2:ncol(df)]
  df = as.matrix(df , ncol = n_hoods , nrow = n_genes)
  rownames(df) = genes
  return(df)
}


#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
.convert_from_sce = function(de_stat , assay_names , coldata_names){

  de_stat = de_stat[ , order(de_stat$Nhood)]

  # combine assays
  de_assays = lapply(assay_names , function(assay_name){
    return(.convert_from_sce_one_var(de_stat , assay_name))
  })
  de_assays = as.data.frame( do.call(cbind , de_assays) )
  colnames(de_assays) = assay_names

  # add gene and Nhood
  df = as.data.frame( assay(de_stat , assay_names[1]) )
  df = rownames_to_column(df , var = "gene")
  df = melt(df , id = "gene")
  colnames(df) = c("gene" , "Nhood" , "var")
  df$Nhood = as.numeric(as.character(df$Nhood))
  df$Nhood = as.integer(df$Nhood)
  df = df[order(df$Nhood) , ]
  df = df[, c("gene" , "Nhood")]

  de_assays = cbind(df , de_assays)
  meta_nhoods = as.data.frame(colData(de_stat))
  de_assays = merge(de_assays , meta_nhoods , by = c("Nhood") , all.x = TRUE)

  return(de_assays)
}


#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
.convert_from_sce_one_var = function(de_stat , var){
  df = as.data.frame( assay(de_stat , var) )
  df = rownames_to_column(df , var = "gene")
  df = melt(df , id = "gene")
  colnames(df) = c("gene" , "Nhood" , "var")
  df$Nhood = as.numeric(as.character(df$Nhood))
  df$Nhood = as.integer(df$Nhood)
  df = df[order(df$Nhood) , ]
  return(df$var)
}


