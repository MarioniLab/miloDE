

#' de_test_neighbourhoods
#'
#' Tests neighbourhoods for DE + performs per gene correction across neighbourhoods
#' @param sce Milo object
#' @param sample_id Character specifying which variable should be used as a sample/replica id. Should be in colData(sce)
#' @param condition_id Character specifying which variable should be used as a condition id. Should be in colData(sce)
#' @param subset_nhoods NULL or character vector specifying the set of neighbourhoods that will be tested for DE
#' @param covariates Vector specifying if additional covariates should be passed into experimental design. Default = NULL (no covariates)
#' @param min_n_cells_per_sample Positive integer specifying the minimum number of cells per replica to be included in testing. Default = 2
#' @param gene_selection In {"none" , "per_neighbourhood"}. If == "per_neighbourhood", it will further discard genes not relevant for this particular neighbourhood.
#' If == "all", it will perform overall selection and then test each hood across same genes (default). If "none" = no selection based on expression.
#' @param genes_2_test Character vector of genes for which we will perform DE testing. Default = rownames(sce). Only applicable if gene_selection == "none".
#' @param min_count Positive integer, specifying min.count for gene selection. Default = 3.
#' @param output_type In {"data.frame","SCE"} Specifying the output format - either in data.frame or sce (with assays corresponding to logFC, pvals (raw and corrected)); columns correspond to neighbourhoods. Note that default is data.frame, but if number of genes x neighbourhoods combinations is > 10^8 -- the output will be sce.
#' @param plot_summary_stat Boolean specifying if we plot Milo neighbourhood plot summarising per neighbourhood whether testing was performed
#' @param layout A character indicating the name of the \code{reducedDim} slot in the \code{\linkS4class{Milo}} object to use for layout (default: 'UMAP'). Only relevant if plot_summary_stat == TRUE.
#' @param BPPARAM NULL or BPPARAM object
#' @return
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment counts
#' @importFrom SummarizedExperiment assay colData
#' @importFrom edgeR DGEList filterByExpr
#' @importFrom tibble rownames_to_column
#' @importFrom scuttle summarizeAssayByGroup
#' @importFrom reshape2 dcast
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel bplapply
#' @importFrom ggpubr ggarrange
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
#' de_stat = de_test_neighbourhoods(sce , condition_id = "type" , gene_selection = "none" , plot_summary_stat = FALSE)
#' de_stat = convert_de_stat(de_stat)
#' de_stat = convert_de_stat(de_stat)
de_test_neighbourhoods = function(sce ,
                                  sample_id = "sample",
                                  condition_id ,
                                  subset_nhoods = NULL,
                                  covariates = NULL,
                                  min_n_cells_per_sample = 2,
                                  gene_selection = c("none","per_neighbourhood") ,
                                  genes_2_test = rownames(sce) ,
                                  min_count = 3 ,
                                  output_type = "data.frame" ,
                                  plot_summary_stat = FALSE,
                                  layout = "UMAP",
                                  BPPARAM = NULL){


  out = .check_argument_correct(sce, .check_sce, "Check sce - something is wrong (gene names unique? reducedDim.name is not present?)") &
    .check_sce_milo(sce) &
    .check_argument_correct(sample_id, is.character, "Check sample_id - should be character vector") &
    .check_argument_correct(condition_id, is.character, "Check condition_id - should be character vector") &
    .check_argument_correct(covariates, .check_string_or_null, "Check covariates - should be NULL or character vector") &
    .check_argument_correct(min_n_cells_per_sample, .check_positive_integer, "Check min_n_cells_per_sample - should be positive integer") &
    .check_argument_correct(gene_selection, function(x) .check_arg_within_options(x, c("none", "per_neighbourhood")),
                            "Check gene_selection - should be either 'none' or 'per_neighbourhood'") &
    .check_argument_correct(genes_2_test, .check_string_or_null, "Check genes_2_test - should be NULL or character vector") &
    .check_argument_correct(min_count, .check_positive_integer, "Check min_count - should be positive integer") &
    .check_argument_correct(output_type, function(x) .check_arg_within_options(x, c("data.frame", "SCE")),
                            "Check output_type - should be either 'data.frame' or 'SCE'") &
    .check_argument_correct(plot_summary_stat, .check_boolean, "Check plot_summary_stat - should be either TRUE or FALSE") &
    .check_var_in_coldata_sce(sce , sample_id , "sample_id") & .check_condition_in_coldata_sce(sce , condition_id) &
    .check_genes_in_sce(sce , genes_2_test) &
    .check_covariates_in_coldata_sce(sce , covariates) & .check_sample_and_condition_id_valid(sce , condition_id , sample_id)

  if (plot_summary_stat){
    out = .check_reducedDim_in_sce(sce , layout)
  }

  nhoods_sce = nhoods(sce)

  if (!is.null(subset_nhoods)){
    out = .check_subset_nhoods(subset_nhoods , nhoods_sce)
  }

  if (gene_selection == "none" & is.null(genes_2_test)){
    stop("If 'gene_selection' == none, genes_2_test should contain at least 1 gene.")
  }

  # assign condition and sample ids
  coldata <- as.data.frame(colData(sce))
  sce$condition_id <- as.character( coldata[, condition_id] )
  sce$sample_id <- as.character( coldata[, sample_id] )

  # delete sample_id and condition_id from covariates
  if (!is.null(covariates)){
    if ("sample_id" %in% covariates){
      warning("Discarding 'sample_id' from covariates since 'sample_id' can not be a covariate name. If in fact you wish to pass a covariate 'sample_id', please rename it first.")
      covariates = setdiff(covariates , "sample_id")
    }
    if ("condition_id" %in% covariates){
      warning("Discarding 'condition_id' from covariates since 'condition_id' can not be a covariate name. If in fact you wish to pass a covariate 'condition_id', please rename it first.")
      covariates = setdiff(covariates , "condition_id")
    }
  }

  # select hoods for testing
  if (is.null(subset_nhoods)){
    subset_nhoods = c(1:ncol(nhoods_sce))
  }
  else {
    if (is.logical(subset_nhoods)){
      subset_nhoods = which(subset_nhoods)
    }
  }
  subset_nhoods = sort(subset_nhoods)


  # get DE for each hood
  if (length(subset_nhoods) == 1){
    warning("You are testing only one neighbourhood. If this is really the intended case, we recommend to run 'de_test_single_neighbourhood' directly.")
    out = de_test_single_neighbourhood(sce , nhoods_sce = nhoods_sce, hood_id = subset_nhoods,
                                       sample_id = sample_id, condition_id = condition_id, covariates = covariates,
                                       min_n_cells_per_sample = min_n_cells_per_sample ,
                                       gene_selection = gene_selection , genes_2_test = genes_2_test, min_count = min_count , run_separately = F)
    out$pval_corrected_across_nhoods = out$pval
    return(out)
  }
  else {
    if (is.null(BPPARAM)) {
      de_stat = lapply(seq(length(subset_nhoods)) , function(i){
        hood_id = subset_nhoods[i]
        out = de_test_single_neighbourhood(sce , nhoods_sce = nhoods_sce, hood_id = hood_id,
                                             sample_id = sample_id, condition_id = condition_id, covariates = covariates,
                                             min_n_cells_per_sample = min_n_cells_per_sample ,
                                             gene_selection = gene_selection , genes_2_test = genes_2_test, min_count = min_count , run_separately = F)
        return(out)
      })
    }
    else {
      de_stat = bplapply(seq(length(subset_nhoods)) , function(i){
        hood_id = subset_nhoods[i]
        out = de_test_single_neighbourhood(sce , nhoods_sce = nhoods_sce, hood_id = hood_id,
                                             sample_id = sample_id, condition_id = condition_id, covariates = covariates,
                                             min_n_cells_per_sample = min_n_cells_per_sample ,
                                             gene_selection = gene_selection , genes_2_test = genes_2_test, min_count = min_count , run_separately = F)
        return(out)
      } , BPPARAM = BPPARAM)
    }

    # put it together in SCE format
    de_stat_sce = SingleCellExperiment(list(logFC = matrix(NA, nrow = nrow(sce), ncol = length(subset_nhoods)) ,
                                              pval = matrix(NA, nrow = nrow(sce), ncol = length(subset_nhoods)) ,
                                              pval_corrected_across_genes = matrix(NA, nrow = nrow(sce), ncol = length(subset_nhoods))))
    rownames(de_stat_sce) = rownames(sce)
    for (i in seq(length(de_stat))){
      current.de_stat = de_stat[[i]]
      assay(de_stat_sce , "logFC")[current.de_stat$gene , i] = current.de_stat$logFC
      assay(de_stat_sce , "pval")[current.de_stat$gene , i] = current.de_stat$pval
      assay(de_stat_sce , "pval_corrected_across_genes")[current.de_stat$gene , i] = current.de_stat$pval_corrected_across_genes

    }

    # add coldata on nhoods
    meta_nhoods = lapply(seq(length(de_stat)) , function(i){
      current.de_stat = de_stat[[i]]
      current.meta = unique(current.de_stat[, c("Nhood" , "Nhood_center" , "sufficient_n_samples" , "design_matrix_suitable")])
      return(current.meta)
    })
    meta_nhoods = do.call(rbind , meta_nhoods)
    #meta_nhoods = meta_nhoods[order(meta_nhoods$Nhood) , ]
    colData(de_stat_sce) = DataFrame(meta_nhoods)
    colnames(de_stat_sce) = meta_nhoods$Nhood

    # get rid of genes that are never used - really only relevant for hood option
    idx = which( rowMeans(is.na(assay(de_stat_sce , "pval"))) < 1)
    if (length(idx) == 0){
      return(NULL)
    } else {
      de_stat_sce = de_stat_sce[idx , ]


      # add pval-corrected-across-nhoods
      ## get weights from nhoods_sce
      #weights = get_weights(nhoods_sce = as.matrix(nhoods_sce[,subset_nhoods]))

      ## calc corrected pvals
      if (is.null(BPPARAM)){
        pval_corrected = lapply(rownames(de_stat_sce) , function(gene){
          pvals = assay(de_stat_sce , "pval")[gene , ]
          out = spatial_pval_adjustment(nhoods_sce = as.matrix(nhoods_sce[,subset_nhoods]) , pvalues = pvals)
          return(out)
        })
      } else {
        pval_corrected = bplapply(rownames(de_stat_sce) , function(gene){
          pvals = assay(de_stat_sce , "pval")[gene , ]
          out = spatial_pval_adjustment(nhoods_sce = as.matrix(nhoods_sce[,subset_nhoods]) , pvalues = pvals)
          return(out)
        }, BPPARAM = BPPARAM)
      }

      ## add corrected to assay
      assay_pval_corrected_across_nhoods = do.call(rbind , pval_corrected)
      rownames(assay_pval_corrected_across_nhoods) = rownames(de_stat_sce)
      colnames(assay_pval_corrected_across_nhoods) = colnames(de_stat_sce)
      assay(de_stat_sce , "pval_corrected_across_nhoods") = assay_pval_corrected_across_nhoods


      # print output message
      message(paste0("DE testing is done.\n",
                     sum(meta_nhoods$sufficient_n_samples), " neighbourhoods (out of ", nrow(meta_nhoods) , ") can be tested wo/ included covariates.\n",
                     sum(meta_nhoods$design_matrix_suitable), " neighbourhoods (out of ", nrow(meta_nhoods) , ") can be tested w/ included covariates."))

      # plot output
      if (plot_summary_stat){
        cols = c("brown2" , "goldenrod1", "seagreen3")
        names(cols) = c("Testing not possible" , "Testing possible wo/ covariates" , "Testing possible w/ covariates")
        meta_nhoods = meta_nhoods[order(meta_nhoods$Nhood) , ]
        meta_nhoods$class = NaN
        meta_nhoods$class[!meta_nhoods$sufficient_n_samples] = "Testing not possible"
        meta_nhoods$class[meta_nhoods$sufficient_n_samples & !meta_nhoods$design_matrix_suitable] = "Testing possible wo/ covariates"
        meta_nhoods$class[meta_nhoods$sufficient_n_samples & meta_nhoods$design_matrix_suitable] = "Testing possible w/ covariates"

        p = plot_milo_by_single_metric(sce, meta_nhoods, colour_by = "class" , layout = layout , size_range = c(1.5,3) , edge_width = c(0.2,0.5)) +
          scale_fill_manual(values = cols, name = "Is testing possible?")
        print(p)
      }


      # return output
      if (output_type == "data.frame"){
        out = convert_de_stat(de_stat_sce)
      } else {
        out = de_stat_sce
      }
      return(out)
    }
  }
}



#' de_test_single_neighbourhood
#'
#' Tests single hood for DE. Not intended to be used by itself (however possible), but rather as a part of `de_test_all_hoods`
#' @param sce Milo object
#' @param nhoods_sce Can be extracted from sce as nhoods(sce) prior to running the function. The reason the argument is passed by itself is to avoid calculating it every time.
#' @param hood_id Character specifying for which hood we should perform testing. Should be in colnames(nhoods_sce)
#' @param sample_id Character specifying which variable should be used as a sample/replica id. Should be in colData(sce)
#' @param condition_id Character specifying which variable should be used as a condition id. Should be in colData(sce)
#' @param covariates Vector specifying if additional covariates should be passed into experimental design. Default = NULL (no covariates)
#' @param min_n_cells_per_sample positive integer specifying the minimun number of cells per replica to be included in testing. Default = 2
#' @param gene_selection In {"per_neighbourhood" , "none"}. If = "per_neighbourhood", we will further discard genes not relevant for this particular neighbourhood.
#' @param genes_2_test Character vector of genes for which we will perform DE testing. Default = rownames(sce). Only applicable if gene_selection = "none"
#' @param min_count Positive integer, specifying min.count for gene selection. Default = 5
#' @param run_separately A boolean parameter specifying whether the function is to be run as a part of 'de_test_all_hoods'(F) or as a stand-alone run (T)
#'
#' @return
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment counts
#' @importFrom SummarizedExperiment assay colData
#' @importFrom edgeR DGEList glmQLFit glmQLFTest calcNormFactors filterByExpr estimateDisp topTags
#' @importFrom tibble rownames_to_column
#' @importFrom scuttle summarizeAssayByGroup
#' @importFrom stats model.matrix p.adjust as.formula
#' @examples
#' require(SingleCellExperiment)
#' require(miloR)
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
#' nhoods_sce = nhoods(sce)
#' de_stat = de_test_single_neighbourhood(sce , nhoods_sce = nhoods_sce, hood_id = 1 , sample_id = "sample" , condition_id = "type")
#'
de_test_single_neighbourhood = function(sce , nhoods_sce , hood_id , sample_id , condition_id , covariates = NULL,
                                        min_n_cells_per_sample = 1 , gene_selection = "none" , genes_2_test = rownames(sce) ,
                                        min_count = 3 , run_separately = TRUE ){

  if (!hood_id %in% 1:ncol(nhoods_sce)){
    stop("'hood_id' should be in 1:ncol(nhoods(sce))")
    return(FALSE)
  }
  if (gene_selection == "none" & is.null(genes_2_test)){
    stop("If 'gene_selection' == none, genes_2_test should contain at least 1 gene.")
  }

  # if 'de_test_single_hood' is a part of the the main run - these steps have already been carried out before
  if (run_separately){

    out = .check_argument_correct(sce, .check_sce, "Check sce - something is wrong (gene names unique? reducedDim.name is not present?)") &
      .check_sce_milo(sce) &
      .check_argument_correct(sample_id, is.character, "Check sample_id - should be character vector") &
      .check_argument_correct(condition_id, is.character, "Check condition_id - should be character vector") &
      .check_argument_correct(covariates, .check_string_or_null, "Check covariates - should be NULL or character vector") &
      .check_argument_correct(min_n_cells_per_sample, .check_positive_integer, "Check min_n_cells_per_sample - should be positive integer") &
      .check_argument_correct(gene_selection, function(x) .check_arg_within_options(x, c("none", "per_neighbourhood")),
                              "Check gene_selection - should be either 'none' or 'per_neighbourhood'") &
      .check_argument_correct(genes_2_test, .check_string_or_null, "Check genes_2_test - should be NULL or character vector") &
      .check_argument_correct(min_count, .check_positive_integer, "Check min_count - should be positive integer") &
      .check_argument_correct(run_separately, .check_boolean, "Check run_separately - should be either TRUE or FALSE") &
      .check_var_in_coldata_sce(sce , sample_id , "sample_id") & .check_condition_in_coldata_sce(sce , condition_id) &
      .check_genes_in_sce(sce , genes_2_test) &
      .check_covariates_in_coldata_sce(sce , covariates) & .check_sample_and_condition_id_valid(sce , condition_id , sample_id)


    coldata <- as.data.frame(colData(sce))
    # assign condition and sample ids
    sce$condition_id <- as.character( coldata[, condition_id] )
    sce$sample_id <- as.character( coldata[, sample_id] )

    # delete sample_id and condition_id from covariates
    if (!is.null(covariates)){
      if ("sample_id" %in% covariates){
        warning("Discarding 'sample_id' from covariates since 'sample_id' can not be a covariate name. If in fact you wish to pass a covariate 'sample_id', please rename it first.")
        covariates = setdiff(covariates , "sample_id")
      }
      if ("condition_id" %in% covariates){
        warning("Discarding 'condition_id' from covariates since 'condition_id' can not be a covariate name. If in fact you wish to pass a covariate 'condition_id', please rename it first.")
        covariates = setdiff(covariates , "condition_id")
      }
    }
  }


  # select cells in the hood
  current.cells = which(nhoods_sce[,hood_id] == 1)
  current.sce = sce[,current.cells]
  current.sce = .filter_samples_with_low_n_cells_in_hood(current.sce , min_n_cells_per_sample = min_n_cells_per_sample)


  if (ncol(current.sce) > 0){
    summed = summarizeAssayByGroup(counts(current.sce), colData(current.sce)[,c("condition_id", "sample_id" , covariates)])
    summed = SingleCellExperiment(list(counts=assay(summed, "sum")), colData=colData(summed))
    y <- DGEList(counts(summed), samples=colData(summed))

    # select genes for testing
    keep <- filterByExpr(y, group=summed$sample_id , min.count = min_count , min.total.count = round(min_count * 1.5))
    if (gene_selection == "none"){
      keep[names(keep) %in% genes_2_test] = TRUE
    }

    if (sum(keep) == 0){
      stop(paste0("For hood_id " , hood_id , " 0 genes are selected for testing. Check that 'min.count' is of the appropriate value and possibly decrease it?"))
      return(NULL)
    } else {
      y <- y[keep,]
    }

    y <- calcNormFactors(y)

    # assess sample composition
    tab = table(current.sce$sample_id , current.sce$condition_id)

    if (ncol(tab) < 2){
      stat = .return_null_stat(rownames(sce) , sufficient_n_samples = FALSE, design_matrix_suitable = FALSE)
    }
    else if (sum(tab[,1] > 0) == 0 | sum(tab[,2] > 0) == 0 | sum(tab > 0) <= 3){
      stat = .return_null_stat(rownames(sce) , sufficient_n_samples = FALSE , design_matrix_suitable = FALSE)
    }
    else {
      stat = tryCatch(
        {
          if (!is.null(covariates)){
            design <- model.matrix(as.formula( paste("~ ", paste(covariates, collapse="+"),sep = "" , " + condition_id") ) , y$samples)
          }
          else {
            design <- model.matrix(~condition_id , y$samples)
          }
          y <- estimateDisp(y, design)
          fit <- glmQLFit(y, design, robust=TRUE)

          res <- glmQLFTest(fit, coef=ncol(design))
          out = as.data.frame(topTags(res, n = Inf ) )

          # if gene_selection == "none" - select only relevant entries and recalculate p-value corrected across genes
          if (gene_selection == "none"){
            out = out[rownames(out) %in% genes_2_test , ]
            out$FDR = p.adjust(out$PValue,"fdr")
          }

          out = out[order(rownames(out)) , c("logFC" , "logCPM" , "PValue" , "FDR")]
          out = rownames_to_column(out , var = "gene")
          out$sufficient_n_samples = TRUE
          out$design_matrix_suitable = TRUE
          out
        },
        error=function(err){
          warning(paste0("For hood_id " , hood_id , " test can not be performed with given design matrix. Reconsider the design matrix (included covariates)."))
          out = .return_null_stat(rownames(sce) , sufficient_n_samples = TRUE , design_matrix_suitable = FALSE)
          return(out)
        }
      )
    }
    stat$Nhood = hood_id
    stat$Nhood_center = colnames(nhoods_sce)[hood_id]

    # clean output
    stat = stat[, c("Nhood", "Nhood_center", "gene", "logFC" , "PValue" , "FDR" , "sufficient_n_samples" , "design_matrix_suitable")]
    colnames(stat)[colnames(stat) == "PValue"] = "pval"
    colnames(stat)[colnames(stat) == "FDR"] = "pval_corrected_across_genes"
    stat = stat[, c("gene" , "Nhood" , "Nhood_center", "sufficient_n_samples" , "design_matrix_suitable" ,
                    "logFC" , "pval" , "pval_corrected_across_genes")]
    if (run_separately){
      stat$pval_corrected_across_nhoods = stat$pval
    }
    return(stat)
  }
  else {
    message(paste0("For neighbourhood " , hood_id , " no samples have number of cells higher than specified threshold. Consider reducing 'min_n_cells_per_sample'?"))
    stat = .return_null_stat(rownames(sce) , sufficient_n_samples = FALSE, design_matrix_suitable = FALSE)

    # clean output
    stat$Nhood = hood_id
    stat$Nhood_center = colnames(nhoods_sce)[hood_id]
    stat = stat[, c("Nhood", "Nhood_center", "gene", "logFC" , "PValue" , "FDR" , "sufficient_n_samples" , "design_matrix_suitable")]
    colnames(stat)[colnames(stat) == "PValue"] = "pval"
    colnames(stat)[colnames(stat) == "FDR"] = "pval_corrected_across_genes"
    stat = stat[, c("gene" , "Nhood" , "Nhood_center", "sufficient_n_samples" , "design_matrix_suitable" ,
                    "logFC" , "pval" , "pval_corrected_across_genes")]

    return(stat)
  }
}




#'
.filter_samples_with_low_n_cells_in_hood = function(sce , min_n_cells_per_sample = 1){
  tab = table(sce$sample_id)
  samples_2_keep = names(tab)[tab >= min_n_cells_per_sample]
  sce = sce[, sce$sample_id %in% samples_2_keep]
  return(sce)
}






#'
.return_null_stat = function(genes , sufficient_n_samples = c(TRUE,FALSE) , design_matrix_suitable = c(TRUE,FALSE)){
  genes = genes[order(genes)]
  stat = data.frame(logFC = NaN, logCPM = NaN, PValue = NaN, FDR = NaN ,
                    sufficient_n_samples = sufficient_n_samples,
                    design_matrix_suitable = design_matrix_suitable)
  stat = stat[rep(seq_len(nrow(stat)), each = length(genes)), ]
  stat$gene = genes
  stat = stat[, c("gene", "logFC" , "logCPM" , "PValue" , "FDR" , "sufficient_n_samples" , "design_matrix_suitable")]
  return(stat)
}


