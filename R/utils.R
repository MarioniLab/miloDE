# Contains various check functions to examine whether variables are of the right format



### checks ###

# sce-milo check: that is used in de_test_all_hoods, de_test_single_hood, filter_hoods, subset_milo
#' @importFrom miloR nhoods graph graph<- nhoods<-
.check_sce_milo = function(sce){
  if (!is(sce , "Milo")){
    stop("SCE should be a Milo object. Please run `assign_neighbourhoods` first.")
    return(F)
  } else if (isEmpty(graph(sce))){
    stop("SCE should contain non-trivial graph. Please run `assign_neighbourhoods` first.")
    return(F)
  } else if ( nrow(nhoods(sce)) == 1 & ncol(nhoods(sce)) == 1 ){
    stop("SCE should contain non-trivial nhoods. Please run `assign_neighbourhoods` first.")
    return(F)
  } else {
    return(T)
  }
}

# sce check: that is used in add_embedding
#' @importFrom SingleCellExperiment reducedDimNames
#' @importFrom SummarizedExperiment assays<- assays
#' @importFrom miloR Milo
.check_sce = function(sce){
  if (!is(sce , "SingleCellExperiment") & !is(sce , "Milo")){
    stop("sce should be a SingleCellExperiment object. If you are working with different format, please first convert to SingleCellExperiment object.")
    return(F)
  } else if (!("counts" %in% names(assays(sce)))){
    stop("sce should contain 'counts' assay that will be used to calculate DE. If counts are stored in different assay, please move them to slot 'counts'.")
    return(F)
  } else if (length(unique(rownames(sce))) < nrow(sce) ){
    stop("sce should have unique rownames.")
    return(F)
  } else if (!is.null(colnames(sce)) & length(unique(colnames(sce))) < ncol(sce)){
    stop("If colnames(sce) exist, they should be unique.")
    return(F)
  } else {
    return(T)
  }
}


#' @importFrom miloR Milo buildGraph graph<- graph nhoods<- nhoodIndex<- buildNhoodGraph
.check_sce_milo = function(sce_milo){
  if (!is(sce_milo , "Milo")){
    stop("sce_milo should be a Milo object. Run 'assign_hoods' first.")
    return(F)
  } else if (length(miloR::graph(sce_milo)) == 0){
    stop("sce_milo should have calculated graph. Run 'assign_hoods' first.")
    return(F)
  } else if (sum(sum(nhoods(sce_milo))) == 0){
    stop("sce_milo should have calculated nhoods. Run 'assign_hoods' first.")
    return(F)
  } else if (isEmpty(nhoodIndex(sce_milo))){
    stop("sce_milo should have calculated nhoodIndex Run 'assign_hoods' first.")
    return(F)
  } else {
    return(T)
  }
}




#' @importFrom SummarizedExperiment assayNames
.check_assay_in_sce = function(sce , assay.type){
  if (.check_sce(sce)){
    if (!assay.type %in% assayNames(sce)){
      stop("assay.type should be in assayNames(sce)")
      return(F)
    }
    else {
      return(T)
    }
  }
}


#' @importFrom SingleCellExperiment colData
.check_condition_in_coldata_sce = function(sce , condition.id){
  if (.check_sce(sce)){
    if (!(condition.id %in% colnames(colData(sce)))){
      stop("condition.id should be in colData(sce)")
      return(F)
    }
    else {
      meta = as.data.frame(colData(sce))
      tab = table(meta[, condition.id])
      if (!(length(tab) == 2)){
        stop("sce should have exactly 2 levels for the condition.id (control and case).")
        return(F)
      } else{
        return(T)
      }
    }
  }
}



#' @importFrom SingleCellExperiment colData
.check_sample_in_coldata_sce = function(sce , sample.id){
  if (.check_sce(sce)){
    if (!(sample.id %in% colnames(colData(sce)))){
      stop("sample.id should be in colData(sce)")
      return(F)
    }
    else {
      return(T)
    }
  }
}


#' @importFrom SingleCellExperiment colData
.check_covariates_in_coldata_sce = function(sce , covariates){
  if (.check_sce(sce)){
    coldata = colnames(colData(sce))
    for (covariate in covariates){
      if (!covariate %in% coldata){
        stop(paste0(covariate , " covariate not in colnames(colData(sce)). Please ensure that all covariates are present."))
        return(F)
      }
      else {
        tab = table(colData(sce)[, covariate])
        if (length(tab)){
          stop(paste0("Covariate '" , covariate , "' has only 1 contrast. Please exclude it prior to testing."))
          return(F)
        }
        else {
          return(T)
        }
      }
    }
  }
}


.check_cell_id_in_sce = function(sce , cell.id){
  if (is.null(cell.id) & is.null(colnames(sce))){
    stop("If colnames(sce) are NULL, cell.id has to be specified in order to assgin unique cell identifiers.")
    return(F)
  } else {
    if (is.null(colnames(sce))){
      if (!cell.id %in% colnames(colData(sce))){
        stop("cell.id should be in colData(sce)")
        return(F)
      }
      else {
        return(T)
      }
    }
    else {
      return(T)
    }
  }
}



.check_genes_in_sce = function(sce, genes){
  if (.check_sce(sce)){
    if (!is.null(genes)){
      out = mean(genes %in% rownames(sce))
      if (out < 1){
        stop("Some gene names are missing from sce")
        return(F)
      }
      else {
        return(T)
      }
    }
    else {
      return(T)
    }
  }
  else {
    return(F)
  }
}


#' @importFrom SingleCellExperiment reducedDimNames
.check_reducedDim_in_sce = function(sce , reducedDim.name){
  if (.check_sce(sce)){
    if (!reducedDim.name %in% reducedDimNames(sce)){
      stop("reducedDim.name should be in reducedDimNames(sce). If you do not have embedding precalculated, run first 'add_embedding'.")
      return(F)
    }
  }
}



.check_cells_ref_and_query = function(cells_sce , cells_ref , cells_query){
  if (.check_sce(sce)){
    if (mean(cells_ref %in% cells_sce) < 1){
      stop("Some of cells_ref are not present.")
      return(F)
    } else if (mean(cells_query %in% cells_sce) < 1){
      stop("Some of cells_query are not present.")
      return(F)
    } else if (!isEmpty(intersect(cells_ref , cells_query))){
      stop("cells_ref and cells_query can not overlap.")
      return(F)
    } else {
      return(T)
    }
  }
  else {
    return(F)
  }

}



.check_argument_correct = function(dots, arg_name , fun , message){
  if (arg_name %in% names(dots)){
    arg = dots[[which(names(dots) == arg_name)]]
    out = fun(arg)
    if (!out){
      stop(message)
    }
    return(out)
  }
  else {
    return(TRUE)
  }
}


.check_quantile_vec = function(quantile.vec){
  if (!is.numeric(quantile.vec)){
    stop("Check quantile.vec - should be numeric vector")
    return(F)
  } else {
    quantile.vec = sort(unique(quantile.vec))
    if (min(quantile.vec) < 0 | max(quantile.vec) > 1){
      stop("quantile.vec should have all its values between 0 and 1. Please enter valid quantile.vec.")
      return(F)
    }
    else {
      return(T)
    }
  }
}

.check_k_grid = function(k.grid){
  if (!is.numeric(k.grid)){
    stop("Check k.grid - should be numeric vector")
    return(F)
  } else {
    k.grid = sort(unique(k.grid))
    if (min(k.grid) < 0){
      stop("Values of k.grid should be positive integers. Please enter valid quantile.vec.")
      return(F)
    }
    else {
      if (length(k.grid) == 1){
        warning("You only selected one value for k. If it is intended, we recommend to run directly 'assign_hoods'")
      }
      if (max(k.grid) >= 1000){
        warning("The highest selected value is > 1000. It is gonna cost computationally, and we generally do not recommend
              such high k. Consider reducing.")
      }
      return(T)
    }
  }
}




.general_check_arguments = function(dots){
  out = TRUE
  out = .check_argument_correct(dots, "sce", .check_sce, "Check sce - something is wrong (gene names unique? reducedDim.name is not present?)")
  out = .check_argument_correct(dots, "sce_milo", .check_sce_milo, "Check sce_milo - something is wrong. Calculate 'assign_hoods' first.)")
  out = .check_argument_correct(dots, "genes", .check_string_or_null, "Check genes - should be NULL or character vector")
  out = .check_argument_correct(dots, "genes_2_exclude", .check_string_or_null, "Check genes_2_exclude - should be NULL or character vector")
  out = .check_argument_correct(dots, "n_hvgs", .check_positive_integer, "Check n_hvgs - should be positive integer")
  out = .check_argument_correct(dots, "assay.type", function(x) .check_arg_within_options(x, c("counts", "logcounts")),
                                "Check assay.type - should be either 'counts' or 'logcounts'")
  out = .check_argument_correct(dots, "reduction_type", function(x) .check_arg_within_options(x, c("Azimuth", "MNN")),
                                "Check reduction_type - should be either 'Azimuth' or 'MNN'")
  out = .check_argument_correct(dots, "reducedDim.name", is.character, "Check reducedDim.name - should be character vector")
  out = .check_argument_correct(dots, "sample.id", is.character, "Check sample.id - should be character vector")
  out = .check_argument_correct(dots, "condition.id", is.character, "Check sample.id - should be character vector")
  out = .check_argument_correct(dots, "cell.id", .check_string_or_null, "Check cell.id - should be NULL or string")
  out = .check_argument_correct(dots, "d", .check_positive_integer, "Check d - should be positive integer")
  out = .check_argument_correct(dots, "order", function(x) .check_arg_within_options(x, c(1,2)),
                                "Check order - should be either 1 (standard kNN-graph) or 2 (2nd-order kNN-graph)")
  out = .check_argument_correct(dots, "k", .check_positive_integer, "Check k - should be positive integer")
  out = .check_argument_correct(dots, "k_init", .check_positive_integer, "Check k_init - should be positive integer")
  out = .check_argument_correct(dots, "prop", .check_prop, "Check prop - should be positive number between 0 and 1")
  out = .check_argument_correct(dots, "filtering", .check_boolean, "Check filtering - should be either TRUE or FALSE")
  out = .check_argument_correct(dots, "k.grid", is.numeric, "Check k.grid - should be numeric vector")
  out = .check_argument_correct(dots, "quantile.vec", is.numeric, "Check quantile.vec - should be numeric vector")
  out = .check_argument_correct(dots, "discard_not_perturbed_hoods", .check_boolean, "Check discard_not_perturbed_hoods - should be either TRUE or FALSE")
  out = .check_argument_correct(dots, "gene_selection", function(x) .check_arg_within_options(x, c("all", "none", "per_hood")),
                                "Check gene_selection - should be either 'all', 'none' or 'per_hood'")
  out = .check_argument_correct(dots, "min_n_cells_per_sample", .check_positive_integer, "Check min_n_cells_per_sample - should be positive integer")
  out = .check_argument_correct(dots, "min.count", .check_positive_integer, "Check min.count - should be positive integer")
  out = .check_argument_correct(dots, "run_separately", .check_boolean, "Check run_separately - should be either TRUE or FALSE")
  out = .check_argument_correct(dots, "covariates", .check_string_or_null, "Check covariates - should be NULL or character vector")
  return(out)
}


.check_positive_integer = function(x){
  out = TRUE
  if (!is.numeric(x)){
    out = FALSE
  } else if (!x%%1 == 0 | x <= 0){
    out = FALSE
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
