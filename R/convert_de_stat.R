

#' convert_de_stat
#'
#' Converts output of milo-DE between \code{\link[base]{data.frame}} and \code{\linkS4class{SingleCellExperiment}} formats
#' @param de_stat Output of milo-DE (from \code{de_test_neighbourhoods}), either in \code{data.frame} or \code{SingleCellExperiment}
#'
#' @return
#' @export
#' @examples
#' de_stat = expand.grid(gene = paste0("gene" , c(1:5)) , Nhood = c(1:10))
#' de_stat$Nhood_id = paste0("nhood_" , de_stat$Nhood)
#' de_stat$logFC = sample(seq(-2,2,1) , nrow(de_stat) , 1)
#' de_stat$pval = sample(c(0,1),nrow(de_stat),1)
#' de_stat$pval_corrected_across_genes = sample(c(0,1),nrow(de_stat),1)
#' de_stat$pval_corrected_across_nhoods = sample(c(0,1),nrow(de_stat),1)
#' de_stat$sufficient_n_samples = TRUE
#' de_stat$design_matrix_suitable = TRUE
#' de_stat = convert_de_stat(de_stat)
#' de_stat = convert_de_stat(de_stat)
#'
convert_de_stat = function(de_stat){

  out = .check_de_stat_valid(de_stat)

  if (class(de_stat) == "SingleCellExperiment"){
    message("Converting to 'data.frame' format")
    de_stat = .convert_from_sce(de_stat)
  } else if (class(de_stat) == "data.frame"){
    message("Converting to 'SingleCellExperiment' format")
    de_stat = .convert_from_df(de_stat)
  }

  return(de_stat)

}


#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
.convert_from_df = function(de_stat){

  df_logFC = .convert_from_df_one_var(de_stat , "logFC")
  df_pval = .convert_from_df_one_var(de_stat , "pval")
  df_pval_corrected_across_genes = .convert_from_df_one_var(de_stat , "pval_corrected_across_genes")
  df_pval_corrected_across_nhoods = .convert_from_df_one_var(de_stat , "pval_corrected_across_nhoods")

  meta_nhoods = unique(de_stat[, c("Nhood" , "Nhood_id" , "sufficient_n_samples" , "design_matrix_suitable")])
  meta_nhoods = meta_nhoods[order(meta_nhoods$Nhood) , ]
  de_stat = SingleCellExperiment(list(logFC = df_logFC ,
                                      pval = df_pval ,
                                      pval_corrected_across_genes = df_pval_corrected_across_genes ,
                                      pval_corrected_across_nhoods = df_pval_corrected_across_nhoods ) ,
                                 colData = DataFrame(meta_nhoods))
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


#' @importFrom SingleCellExperiment SingleCellExperiment colData
.convert_from_sce = function(de_stat){

  df_logFC = .convert_from_sce_one_var(de_stat , "logFC")
  df_pval = .convert_from_sce_one_var(de_stat , "pval")
  df_pval_corrected_across_genes = .convert_from_sce_one_var(de_stat , "pval_corrected_across_genes")
  df_pval_corrected_across_nhoods = .convert_from_sce_one_var(de_stat , "pval_corrected_across_nhoods")
  meta_nhoods = as.data.frame(colData(de_stat))

  df = cbind(df_logFC , df_pval[, "pval"] , df_pval_corrected_across_genes[, "pval_corrected_across_genes"] ,
             df_pval_corrected_across_nhoods[, "pval_corrected_across_nhoods"] ,
             meta_nhoods[,c("Nhood_id" , "sufficient_n_samples" , "design_matrix_suitable")])

  colnames(df) = c("gene" , "Nhood" , "logFC", "pval" , "pval_corrected_across_genes" ,
                   "pval_corrected_across_nhoods" , "Nhood_id" , "sufficient_n_samples" , "design_matrix_suitable")
  df = df[, c("gene" , "Nhood", "Nhood_id" , "logFC" , "pval" , "pval_corrected_across_genes" , "pval_corrected_across_nhoods" ,
                                                       "sufficient_n_samples" , "design_matrix_suitable")]
  return(df)
}


#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
.convert_from_sce_one_var = function(de_stat , var){
  df = as.data.frame( assay(de_stat , var) )
  df = rownames_to_column(df , var = "gene")
  df = melt(df , id = "gene")
  colnames(df) = c("gene" , "Nhood" , var)
  df = df[order(df$Nhood) , ]
  return(df)
}


