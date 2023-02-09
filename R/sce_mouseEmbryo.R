#' Chimeric mouse embryo, Tal1- (Pijuan-Sala et al., 2019)
#'
#' \code{SingleCellExperiment} object, containing counts
#' (raw are pulled using \code{\link[MouseGastrulationData]{Tal1ChimeraData}}
#' and log-normalised are estimated using \code{\link[scuttle]{logNormCounts}}) matrices
#' for chimeric mouse embryos (both Tal1+ and Tal1-).
#'
#' Additionally, we subselected only 3 cell types (Endothelium , Blood progenitors 2, Neural crest)
#' and 300 HVGs.
#'
#' \code{colData(sce_mouseEmbryo)} contains information about replicate (\code{sample}) and Tal1 status (\code{tomato = TRUE} means Tal1-).
#'
#' @docType data
#' @usage data(sce_mouseEmbryo)
#'
#' @format \code{\linkS4class{SingleCellExperiment}} object
#'
#' @name sce_mouseEmbryo
#' @source MouseGastrulationData package
NULL
