# Contains various check functions to examine whether variables are of the right format



### checks ###

# sce-milo check: that is used in de_test_all_hoods, de_test_single_hood, filter_hoods, subset_milo
#' @importFrom miloR Milo nhoods graph graph<- nhoods<- nhoodIndex<- buildNhoodGraph
#' @importFrom methods is
.check_sce_milo = function(x){
  if (!is(x , "Milo")){
    stop("x should be a Milo object. Please run `assign_neighbourhoods` first.")
    return(FALSE)
  } else if (length(miloR::graph(x)) == 0){
    stop("x should contain non-trivial graph. Please run `assign_neighbourhoods` first.")
    return(FALSE)
  } else if (sum(sum(nhoods(x))) == 0){
    stop("x should have calculated nhoods. Please run 'assign_neighbourhoods' first.")
    return(FALSE)
  } else if ( nrow(nhoods(x)) == 1 & ncol(nhoods(x)) == 1 ){
    stop("x should contain non-trivial nhoods. Please run `assign_neighbourhoods` first.")
    return(FALSE)
  } else if (isEmpty(nhoodIndex(x))){
    stop("x should have calculated nhoodIndex. Please run 'assign_neighbourhoods' first.")
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# sce check: that is used in add_embedding
#' @importFrom SingleCellExperiment reducedDimNames
#' @importFrom SummarizedExperiment assays<- assays
#' @importFrom miloR Milo
#' @importFrom methods is
.check_sce = function(x){
  if (!is(x , "SingleCellExperiment") & !is(x , "Milo")){
    stop("x should be a SingleCellExperiment or Milo object.")
    return(FALSE)
  } else if (!("counts" %in% names(assays(x)))){
    stop("x should contain 'counts' assay that will be used to calculate DE. If counts are stored in different assay, please move them to slot 'counts'.")
    return(FALSE)
  } else if (length(unique(rownames(x))) < nrow(x) ){
    stop("x should have unique rownames.")
    return(FALSE)
  } else if (!is.null(colnames(x)) & length(unique(colnames(x))) < ncol(x)){
    stop("If colnames(x) exist, they should be unique.")
    return(FALSE)
  } else {
    return(TRUE)
  }
}





#' @importFrom igraph is_igraph
.valid_nhood <- function(x){
  # check for a valid nhood slot
  n_neigh <- ncol(nhoods(x))
  is_not_empty <- n_neigh > 0
  if (is_not_empty) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @importFrom igraph is_igraph
.valid_graph <- function(x){
  # check for a valid graph
  if(isTRUE(is_igraph(x))){
    return(TRUE)
  } else{
    return(FALSE)
  }
}



#' @importFrom SummarizedExperiment assayNames
.check_assay_in_sce = function(x , assay_type){
  if (.check_sce(x)){
    if (!assay_type %in% assayNames(x)){
      stop("assay_type should be in assayNames(x)")
      return(FALSE)
    }
    else {
      return(TRUE)
    }
  }
  else {
    return(FALSE)
  }
}


#' @importFrom SummarizedExperiment colData
.check_condition_in_coldata_sce = function(x , condition_id){
  if (.check_sce(x)){
    if (!(condition_id %in% colnames(colData(x)))){
      stop("'condition_id' should be in colData(x)")
      return(FALSE)
    }
    else {
      meta = as.data.frame(colData(x))
      tab = table(meta[, condition_id])
      if (length(tab) < 2){
        stop("x should have at least two levels for tested conditions.")
        return(FALSE)
      } else{
        return(TRUE)
      }
    }
  }
  else {
    return(FALSE)
  }
}

#' @importFrom SummarizedExperiment colData
.check_sample_and_condition_id_valid = function(x , condition_id , sample_id){

  if (condition_id == sample_id){
    stop("'sample_id' and 'condition_id' can not be the same")
    return(FALSE)
  }
  else {
    meta = as.data.frame(colData(x))
    tab = table(meta[, sample_id] , meta[, condition_id])
    tab = sapply(1:nrow(tab) , function(i) sum(tab[i,] > 0))
    if (mean(tab == 1) < 1){
      stop("Each sample_id should be associated with one condition")
      return(FALSE)
    }
    else {
      return(TRUE)
    }
  }
}

#' @importFrom SummarizedExperiment colData
.check_var_in_coldata_sce = function( x , var , var_intended){
  if (.check_sce(x)){
    if (!(var %in% colnames(colData(x)))){
      stop(paste0("'", var_intended, "'", " should be in colData(x)"))
      return(FALSE)
    }
    else {
      return(TRUE)
    }
  }
  else {
    return(FALSE)
  }
}


#' @importFrom miloR nhoods
.check_nhood_stat = function(nhood_stat , x){
  if (!is(nhood_stat , "data.frame")){
    stop("'nhood_stat' should be a data.frame.")
    return(FALSE)
  }
  else {
    if (!"Nhood" %in% colnames(nhood_stat)){
      stop("'nhood_stat' should contain column 'Nhood'.")
      return(FALSE)
    }
    else {
      if (!is.numeric(nhood_stat$Nhood)){
        stop("'nhood_stat$Nhood' should be numeric.")
        return(FALSE)
      }
      else {
        nhoods_sce = nhoods(x)
        if (mean(nhood_stat$Nhood %in% c(1:ncol(nhoods_sce))) < 1){
          stop("'nhood_stat$Nhood' should be within c(1:ncol(nhoods(x))).")
          return(FALSE)
        }
        else {
          return(TRUE)
        }
      }
    }
  }
}


#' @importFrom SummarizedExperiment colData
.check_covariates_in_coldata_sce = function(x , covariates){
  if (.check_sce(x)){
    if (is.null(covariates)){
      return(TRUE)
    }
    else {
      coldata = colnames(colData(x))

      covariates_exist = sapply(covariates , function(covariate){
        out = as.numeric(covariate %in% coldata)
        return(out)
      })
      if (sum(covariates_exist) < length(covariates)){
        stop("All covariates should be colnames of colData(x).")
        return(FALSE)
      }
      else {
        covariates_w_contrast = sapply(covariates , function(covariate){
          tab = table(colData(x)[, covariate])
          out = as.numeric(length(tab) != 1)
          return(out)
        })
        if (sum(covariates_w_contrast) < length(covariates)){
          stop("All covariates should have more than 1 contrast.")
          return(FALSE)
        }
        else {
          return(TRUE)
        }
      }
    }
  }
  else {
    return(FALSE)
  }
}




#' @importFrom SummarizedExperiment colData
.check_cell_id_in_sce = function(x , cell_id){
  if (is.null(cell_id) & is.null(colnames(x))){
    stop("If colnames(x) are NULL, cell_id has to be specified in order to assgin unique cell identifiers.")
    return(FALSE)
  } else {
    if (is.null(colnames(x))){
      if (!cell_id %in% colnames(colData(x))){
        stop("cell_id should be in colData(x)")
        return(FALSE)
      }
      else {
        return(TRUE)
      }
    }
    else {
      return(TRUE)
    }
  }
}



.check_genes_in_sce = function(x, genes){
  if (.check_sce(x)){
    if (!is.null(genes)){
      if (mean(genes %in% rownames(x)) < 1){
        stop("Some gene names are missing from x")
        return(FALSE)
      }
      else {
        return(TRUE)
      }
    }
    else {
      return(TRUE)
    }
  }
  else {
    return(FALSE)
  }
}


#' @importFrom SingleCellExperiment reducedDimNames
.check_reducedDim_in_sce = function(x , reducedDim_name){
  if (.check_sce(x)){
    if (!reducedDim_name %in% reducedDimNames(x)){
      stop("reducedDim_name should be in reducedDimNames(x).")
      return(FALSE)
    }
  }
  else {
    return(FALSE)
  }
}


#' @importFrom S4Vectors isEmpty
.check_cells_ref_and_query = function(cells_sce , cells_ref , cells_query){
  if (mean(cells_ref %in% cells_sce) < 1){
    stop("Some of cells_ref are not present.")
    return(FALSE)
  } else if (mean(cells_query %in% cells_sce) < 1){
    stop("Some of cells_query are not present.")
    return(FALSE)
  } else if (!isEmpty(intersect(cells_ref , cells_query))){
    stop("cells_ref and cells_query can not overlap.")
    return(FALSE)
  } else {
    return(TRUE)
  }
}


#' @importFrom SummarizedExperiment assayNames assayNames<- colData
.check_de_stat_valid = function(de_stat , assay_names , coldata_names){
  if (!class(de_stat) %in% c("data.frame" , "SingleCellExperiment")){
    stop("de_stat should be either data.frame or SingleCellExperiment object.\n
         To get valid de_stat object, please run 'de_test_neighbourhoods.R'")
    return(FALSE)
  } else if (length(intersect(assay_names, coldata_names)) > 0){
    stop("assay_names and coldata_names can not overlap")
    return(FALSE)
  } else if (is(de_stat , "data.frame")){
    cols = colnames(de_stat)
    cols_required = c(assay_names , coldata_names)
    if (mean(cols_required %in% cols) < 1){
      stop("colnames(de_stat) missing some of the assay_names or coldata_names.")
      return(FALSE)
    } else if (!is.numeric(de_stat$Nhood)) {
      stop("Nhood field should be numeric. To get valid de_stat object, please run 'de_test_neighbourhoods.R'")
      return(FALSE)
    }
  } else if (is(de_stat , "SingleCellExperiment")){
    cols = assayNames(de_stat)
    if (mean(assay_names %in% cols) < 1){
      stop("de_stat missing some of the required assays.")
      return(FALSE)
    } else {
      meta_nhoods = as.data.frame(colData(de_stat))
      if (mean(coldata_names %in% colnames(meta_nhoods)) < 1){
        stop("de_stat missing some of the coldata_names")
        return(FALSE)
      } else if (!is.numeric(de_stat$Nhood)) {
        stop("colData field Nhood should be numeric. To get valid de_stat object, please run 'de_test_neighbourhoods.R'")
        return(FALSE)
      }
    }
  } else {
    return(TRUE)
  }
}


#
# .check_argument_correct = function(dots, arg_name , fun , message){
#   if (arg_name %in% names(dots)){
#     arg = dots[[which(names(dots) == arg_name)]]
#     out = fun(arg)
#     if (!out){
#       stop(message)
#     }
#     return(out)
#   }
#   else {
#     return(TRUE)
#   }
# }
#

.check_weights = function(weights){
  if (!is.numeric(weights)){
    stop("weights should be a numeric vector. To get valid weights, run 'get_weights.R'")
    return(FALSE)
  } else {
    if (sum(is.na(weights)) > 0){
      stop("weights can not contain NaNs. To get valid weights, run 'get_weights.R'")
      return(FALSE)
    } else if (sum(weights <= 0) > 0){
      stop("weights should be positive. To get valid weights, run 'get_weights.R'")
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}


.check_weights_and_pvals = function(weights , pvalues , nhoods_sce){
  if (!length(weights) == length(pvalues)){
    stop("weights should be of the same size as pvalues.")
    return(FALSE)
  } else if (!length(weights) == ncol(nhoods_sce)){
    stop("weights should be of the same size as number of columns in nhoods_sce.")
    return(FALSE)
  } else {
    return(TRUE)
  }
}




.check_argument_correct = function(arg , fun , message){
  out = fun(arg)
  if (!out){
    stop(message)
  }
  return(out)
}


.check_quantile_vec = function(quantile_vec){
  if (!is.numeric(quantile_vec)){
    stop("Check quantile_vec - should be numeric vector")
    return(FALSE)
  } else {
    quantile_vec = sort(unique(quantile_vec))
    if (min(quantile_vec) < 0 | max(quantile_vec) > 1){
      stop("quantile_vec should have all its values between 0 and 1. Please enter valid quantile_vec")
      return(FALSE)
    }
    else {
      return(TRUE)
    }
  }
}


.check_k_grid = function(k_grid){
  if (!is.numeric(k_grid)){
    stop("Check k_grid - should be numeric vector")
    return(FALSE)
  } else {
    k_grid = sort(unique(k_grid))
    if (min(k_grid) < 0 | max(k_grid%%1 > 0)){
      stop("Values of k_grid should be positive integers. Please enter valid k_grid.")
      return(FALSE)
    }
    else {
      if (length(k_grid) == 1){
        warning("You only selected one value for k. If it is intended, we recommend to run directly 'assign_neighbourhoods'.")
      }
      if (max(k_grid) >= 1000){
        warning("The highest selected value is > 1000. It is gonna cost computationally, and we generally do not recommend such high k. Consider reducing.")
      }
      return(TRUE)
    }
  }
}


.check_pval_thresh = function(pval.thresh){
  if (!is.numeric(pval.thresh)){
    stop("pval.thresh should be numeric")
    return(FALSE)
  }
  else {
    if (!length(pval.thresh) == 1){
      stop("pval.thresh should be a single number")
      return(FALSE)
    }
    else {
      if (pval.thresh <= 0 | pval.thresh >=1){
        stop("pval.thresh should be between 0 and 1")
        return(FALSE)
      }
      else {
        return(TRUE)
      }
    }
  }
}



.check_z_thresh = function(z.thresh){
  if (!is.numeric(z.thresh)){
    stop("z.thresh should be numeric")
    return(FALSE)
  }
  else {
    if (!length(z.thresh) == 1){
      stop("z.thresh should be a single number")
      return(FALSE)
    }
    else {
      if (z.thresh > 0){
        stop("z.thresh should be not higher than 0")
        return(FALSE)
      }
      else {
        return(TRUE)
      }
    }
  }
}


#
# .general_check_arguments = function(dots){
#   out = TRUE
#   out = .check_argument_correct(dots, "sce", .check_sce, "Check sce - something is wrong (gene names unique? reducedDim.name is not present?)")
#   out = .check_argument_correct(dots, "sce_milo", .check_sce_milo, "Check sce_milo - something is wrong. Calculate 'assign_hoods' first.)")
#   out = .check_argument_correct(dots, "genes", .check_string_or_null, "Check genes - should be NULL or character vector")
#   out = .check_argument_correct(dots, "genes_2_exclude", .check_string_or_null, "Check genes_2_exclude - should be NULL or character vector")
#   out = .check_argument_correct(dots, "n_hvgs", .check_positive_integer, "Check n_hvgs - should be positive integer")
#   out = .check_argument_correct(dots, "assay_type", function(x) .check_arg_within_options(x, c("counts", "logcounts")),
#                                 "Check assay_type - should be either 'counts' or 'logcounts'")
#   out = .check_argument_correct(dots, "reduction_type", function(x) .check_arg_within_options(x, c("Azimuth", "MNN")),
#                                 "Check reduction_type - should be either 'Azimuth' or 'MNN'")
#   out = .check_argument_correct(dots, "reducedDim_name", is.character, "Check reducedDim_name - should be character vector")
#   out = .check_argument_correct(dots, "sample_id", is.character, "Check sample_id - should be character vector")
#   out = .check_argument_correct(dots, "condition_id", is.character, "Check condition_id - should be character vector")
#   out = .check_argument_correct(dots, "cell_id", .check_string_or_null, "Check cell_id - should be NULL or string")
#   out = .check_argument_correct(dots, "d", .check_positive_integer, "Check d - should be positive integer")
#   out = .check_argument_correct(dots, "order", function(x) .check_arg_within_options(x, c(1,2)),
#                                 "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)")
#   out = .check_argument_correct(dots, "k", .check_positive_integer, "Check k - should be positive integer")
#   out = .check_argument_correct(dots, "k_init", .check_positive_integer, "Check k_init - should be positive integer")
#   out = .check_argument_correct(dots, "prop", .check_prop, "Check prop - should be positive number between 0 and 1")
#   out = .check_argument_correct(dots, "filtering", .check_boolean, "Check filtering - should be either TRUE or FALSE")
#   out = .check_argument_correct(dots, "k.grid", is.numeric, "Check k.grid - should be numeric vector")
#   out = .check_argument_correct(dots, "quantile_vec", is.numeric, "Check quantile_vec - should be numeric vector")
#   out = .check_argument_correct(dots, "discard_not_perturbed_hoods", .check_boolean, "Check discard_not_perturbed_hoods - should be either TRUE or FALSE")
#   out = .check_argument_correct(dots, "gene_selection", function(x) .check_arg_within_options(x, c("all", "none", "per_hood")),
#                                 "Check gene_selection - should be either 'all', 'none' or 'per_hood'")
#   out = .check_argument_correct(dots, "min_n_cells_per_sample", .check_positive_integer, "Check min_n_cells_per_sample - should be positive integer")
#   out = .check_argument_correct(dots, "min_count", .check_positive_integer, "Check min_count - should be positive integer")
#   out = .check_argument_correct(dots, "run_separately", .check_boolean, "Check run_separately - should be either TRUE or FALSE")
#   out = .check_argument_correct(dots, "covariates", .check_string_or_null, "Check covariates - should be NULL or character vector")
#   out = .check_argument_correct(dots, "seed", .check_number_or_null, "Check seed - should be NULL or number")
#   return(out)
# }


.check_positive_integer = function(x){
  out = TRUE
  if (!is.numeric(x)){
    out = FALSE
  } else if (!x%%1 == 0 | x <= 0){
    out = FALSE
  }
  return(out)
}


.check_non_negative = function(x){
  out = TRUE
  if (!is.numeric(x)){
    out = FALSE
  } else if (x < 0){
    out = FALSE
  }
  return(out)
}

.check_design = function(x){
  if (is(x , "formula")){
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}


.check_number_or_null = function(x){
  out = TRUE
  if (!is.null(x)){
    if (!is.numeric(x)){
      out = FALSE
    }
  }
  return(out)
}


.check_arg_within_options = function(x , options){
  out = TRUE
  if (is.null(x)){
    out = FALSE
  }
  else  if (!x %in% options){
    out = FALSE
  }
  return(out)
}

.check_string_or_null = function(x){
  out = TRUE
  if (!is.null(x)){
    if (!is.character(x)){
      out = FALSE
    }
  }
  return(out)
}

.check_prop = function(x){
  out = TRUE
  if (!is.numeric(x)){
    out = FALSE
  } else if (x <= 0 | x > 1){
    out = FALSE
  }
  return(out)
}


.check_boolean = function(x){
  out = TRUE
  if (is.null(x)){
    out = FALSE
  }
  else if (!x %in% c(TRUE, FALSE)){
    out = FALSE
  }
  return(out)
}



.check_nhoods_matrix = function(nhoods_sce){
  unq_elements = unique(as.numeric(nhoods_sce))
  if (mean(unq_elements %in% c(0,1)) < 1){
    stop("All elements of nhoods matrix should be either 0 or 1. To get valid matrix, run nhoods(x).")
    return(FALSE)
  }
  else {
    return(TRUE)
  }
}


.check_subset_nhoods = function(subset_nhoods, nhoods_sce){
  if (!is.numeric(subset_nhoods) & !is.logical(subset_nhoods)){
    stop("'subset_nhoods' should be either numeric or logical.")
    return(FALSE)
  }
  else if (is.numeric(subset_nhoods)){
    if (mean(subset_nhoods %in% c(1:ncol(nhoods_sce))) < 1){
      stop("If 'subset_nhoods' is numeric vector, it should lie within c(1:ncol(nhoods(x))).")
      return(FALSE)
    }
    else {
        return(TRUE)
    }
  }
  else {
    if (!length(subset_nhoods) == ncol(nhoods_sce)){
      stop("If 'subset_nhoods' is logical vector, it should be the same size as ncol(nhoods(x)).")
      return(FALSE)
    }
    else {
      return(TRUE)
    }
  }
}

#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr distinct
#' @importFrom stats model.matrix
.check_design_and_covariates_match = function(x , design , sample_id , covariates){
  design.df = as.data.frame(colData(x)[,c(sample_id , covariates)])
  design.df = distinct(design.df)
  out = tryCatch(
    {
      design = model.matrix(design , data = design.df)
      TRUE
    },
      error=function(err){
        stop("Some of the design's arguments are not in the covariate vector.")
        return(FALSE)
    }
  )
  return(out)
}


