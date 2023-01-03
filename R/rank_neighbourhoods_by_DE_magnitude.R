

#' rank_neighbourhoods_by_DE_magnitude
#'
#' Ranks neighbourhoods by the magnitude of DE: number of DE genes and number of 'specifically' DE genes. Note that for this analysis we set NaN p-values (raw and corrected) to 1 - interpret accordingly.
#' @param de_stat Output of milo-DE, either in data.frame format or SCE
#' @param pval.thresh A scalar specifying which threshold to use for deciding on significance for gene being DE in a hood. Default = 0.1.
#' @param z.thresh A scalar specifying which threshold to use for deciding on which z-normalised p-values are going to be considered specifically DE. Default = -3.
#'
#'
#' @return
#' @export
#' @importFrom SummarizedExperiment colData
#' @importFrom stats sd
#' @examples
#' de_stat = expand.grid(gene = paste0("gene" , c(1:5)) , Nhood = c(1:10))
#' de_stat$Nhood_center = paste0("nhood_" , de_stat$Nhood)
#' de_stat$logFC = sample(seq(-2,2,1) , nrow(de_stat) , 1)
#' de_stat$pval = sample(c(0,1),nrow(de_stat),1)
#' de_stat$pval_corrected_across_genes = sample(c(0,1),nrow(de_stat),1)
#' de_stat$pval_corrected_across_nhoods = sample(c(0,1),nrow(de_stat),1)
#' de_stat$sufficient_n_samples = TRUE
#' de_stat$design_matrix_suitable = TRUE
#' out = rank_neighbourhoods_by_DE_magnitude(de_stat)
#'
rank_neighbourhoods_by_DE_magnitude = function(de_stat, pval.thresh = 0.1, z.thresh = -3 ){

  out = .check_de_stat_valid(de_stat , assay_names = c("logFC" , "pval" , "pval_corrected_across_nhoods" , "pval_corrected_across_genes") , coldata_names = c("Nhood" , "Nhood_center")) &
    .check_pval_thresh(pval.thresh) & .check_z_thresh(z.thresh)

  if (class(de_stat) == "data.frame"){
    de_stat = convert_de_stat(de_stat ,
                              assay_names = c("logFC" , "pval" , "pval_corrected_across_nhoods" , "pval_corrected_across_genes") ,
                              coldata_names = c("Nhood" , "Nhood_center"))
    de_stat = de_stat[ , order(de_stat$Nhood)]
  }

  # calculate number of DE genes
  assay_pval_corrected_across_genes = assay(de_stat , "pval_corrected_across_genes")
  idx = which(is.na(assay_pval_corrected_across_genes))
  assay_pval_corrected_across_genes[idx] = 1
  out_n_genes = colSums(assay_pval_corrected_across_genes < pval.thresh , na.rm = T)

  # calculate number of specifically DE genes
  assay_pval_corrected_across_nhoods = assay(de_stat , "pval_corrected_across_nhoods")
  idx = which(is.na(assay_pval_corrected_across_nhoods))
  assay_pval_corrected_across_nhoods[idx] = 1
  assay_pval_corrected_across_nhoods = apply(assay_pval_corrected_across_nhoods , 1 , function(x){
    return((x - mean(x, na.rm = T))/sd(x, na.rm = T))
  })
  assay_pval_corrected_across_nhoods = t(assay_pval_corrected_across_nhoods)
  out_n_specific_genes = colSums(assay_pval_corrected_across_nhoods < z.thresh , na.rm = T)

  out = cbind(out_n_genes , out_n_specific_genes)
  colnames(out) = c("n_DE_genes" , "n_specific_DE_genes")

  meta_nhoods = as.data.frame(colData(de_stat))
  meta_nhoods = cbind(meta_nhoods , out)

  return(meta_nhoods)

}




