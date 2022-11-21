


#' add_embedding
#'
#' @param sce SCE object
#' @param genes genes that should be used to calculate joint embedding (usually some sort of HVG selection). If = NULL (default) , we will calculate HVGs using scran::getTopHVGs
#' @param n_hvgs number of HVGs to retain for the joint embedding estimation
#' @param assay.type A string specifying the assay in SCE containing the expression values to be used (default = "logcounts")
#' @param reduction_type Either "Azimuth" or "MNN" - specifying which method to use for the embedding estimation
#' @param reducedDim.name String specifying the name under which to save the joint embedding
#' @param sample.id Character specifying which variable should be used as a sample/replica id. Should be in colData(sce)
#' @param cells_ref IDs of cells from reference samples
#' @param cells_query IDs of cells from query samples
#' @param cell.id Character specifying which variable should be used as a cell id (if colnames(sce) are not specified).
#' @param d Positive integer specifying number of dimensions in the joint embedding
#'
#' @return
#' @export
#'
#' @importFrom scran modelGeneVar getTopHVGs
#' @examples
#' require(SingleCellExperiment)
#' n_row = 500
#' n_col = 100
#' sce = SingleCellExperiment(assays = list( counts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' logcounts(sce) = matrix(rnorm(n_row*n_col), ncol=n_col)
#' sce$sample = floor(runif(n = n_col , min = 1 , max = 5))
#' sce$type = ifelse(sce$sample %in% c(1,2) , "ref" , "query")
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' out = add_embedding(sce, reduction_type = "MNN" , reducedDim.name = "MNN", cells_ref = colnames(sce[, sce$type == "ref"]) ,  cells_query = colnames(sce[, sce$type == "query"]) , d = 5)
add_embedding = function(sce ,
                         genes = NULL,
                         n_hvgs = 3000,
                         assay.type = "logcounts" ,
                         reduction_type = c("Azimuth" , "MNN") ,
                         reducedDim.name ,
                         sample.id = "sample",
                         cells_ref ,
                         cells_query ,
                         cell.id = NULL ,
                         d = 30){

  # check that arguments look alright
  args = c(as.list(environment()))
  out = .general_check_arguments(args) & .check_genes_in_sce(sce , genes) & .check_cell_id_in_sce(sce , cell.id) & .check_sample_in_coldata_sce(sce , sample.id)

  # select only relevant genes for the embedding: wither passed with argument genes or in case genes=NULL (default), we will calculate top `n_hvgs` genes
  if (!is.null(genes)){
    sce = sce[genes , ]
  }
  else {
    dec.sce = modelGeneVar(sce , assay.type = assay.type)
    hvgs = getTopHVGs(dec.sce, n = n_hvgs)
    sce = sce[hvgs , ]
  }

  # if sce doenst have colnames, assign
  if (is.null(colnames)){
    message("Assigning colnames based on provided 'cell.id'.")
    meta = as.data.frame(colData(sce))
    colnames(sce) = meta[, cell.id]
  }
  # check that cells_ref and cells_query belong to colnames(sce) and do not overlap
  out = .check_cells_ref_and_query(colnames(sce) , cells_ref , cells_query)

  # if assay.type == "counts" and reduction_type == "MNN" -> warn the user
  if (assay.type == "counts" & reduction_type == "MNN"){
    warning("For reduction_type == 'MNN', we recommed to use normalised logcounts (either pre-compute and use logcounts or use reduction_type == 'Azimuth').")
  }

  if (reduction_type == "MNN"){
    sce_w_embedding = .add_mnn_based_embedding(sce , assay.type = assay.type , reducedDim.name = reducedDim.name, sample.id = sample.id, cells_ref = cells_ref, cells_query = cells_query, d = d)
  }
  else if (reduction_type == "Azimuth"){
    sce_w_embedding = .add_azimuth_embedding(sce , reducedDim.name = reducedDim.name, sample.id = sample.id, cells_ref = cells_ref, cells_query = cells_query, d = d)
  }
  sce = sce[order(colnames(sce)) , ]
  sce_w_embedding = sce_w_embedding[order(colnames(sce_w_embedding)) , ]
  reducedDim(sce , reducedDim.name) = reducedDim(sce_w_embedding , reducedDim.name)
  return(sce)

}



#' @importFrom batchelor multiBatchPCA reducedMNN
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assays<- assays assay assay<-
.add_mnn_based_embedding = function(sce , assay.type = "logcounts" , reducedDim.name , sample.id = "sample", cells_ref , cells_query , d = 30){
  out = .check_assay_in_sce(sce , assay.type)

  sce_ref = sce[ , colnames(sce) %in% cells_ref]
  sce_query = sce[ , colnames(sce) %in% cells_query]

  meta_ref = as.data.frame(colData(sce_ref))
  meta_query = as.data.frame(colData(sce_query))

  pca_full = multiBatchPCA(sce_ref , batch = factor(meta_ref[, sample.id]), d = d, preserve.single = TRUE, assay.type = assay.type)
  pca_ref <- pca_full[[1]]
  to_proj <- (as.matrix(assay(sce_query , assay.type)) - rowMeans(as.matrix(assay(sce_query , assay.type))))/sqrt(rowVars(as.matrix(assay(sce_query , assay.type))))

  # project
  pca_query <- t(to_proj) %*% metadata(pca_full)$rotation
  pca_joint <- rbind(pca_ref, pca_query)
  rownames(pca_joint) = c(colnames(sce_ref) , colnames(sce_query))

  # batch correct
  pca_joint = reducedMNN(pca_joint , batch = factor(c(meta_ref[, sample.id], meta_query[, sample.id])) )
  pca_joint = pca_joint$corrected

  # reorder to ensure that colnames of sce are aligned with rownames of the embedding
  pca_joint = pca_joint[order(rownames(pca_joint)) , ]
  sce = sce[, order(colnames(sce))]

  reducedDim(sce , reducedDim.name) = pca_joint
  return(sce)
}


#' @importFrom Seurat CreateSeuratObject AddMetaData SplitObject VariableFeatures<- NormalizeData FindIntegrationAnchors IntegrateData ScaleData DefaultAssay DefaultAssay<- Embeddings RunPCA FindTransferAnchors TransferData CreateDimReducObject as.SingleCellExperiment
#' @importFrom SingleCellExperiment reducedDim
.add_azimuth_embedding = function(sce , reducedDim.name , sample.id = "sample", cells_ref , cells_query , d = 30){
  sce_ref = sce[ , colnames(sce) %in% cells_ref]
  sce_query = sce[ , colnames(sce) %in% cells_query]

  meta_ref = as.data.frame(colData(sce_ref))
  meta_query = as.data.frame(colData(sce_query))

  samples_ref = unique(meta_ref[ , sample.id])
  samples_query = unique(meta_query[ , sample.id])

  sce_seurat <- CreateSeuratObject(counts = counts(sce))
  sce_seurat = AddMetaData(sce_seurat, as.data.frame(colData(sce)), col.name = NULL)
  sce_seurat.list <- SplitObject(sce_seurat, split.by = sample.id)
  # normalise
  for (i in 1:length(sce_seurat.list)) {
    sce_seurat.list[[i]] <- NormalizeData(sce_seurat.list[[i]], verbose = FALSE)
    VariableFeatures(sce_seurat.list[[i]]) <- rownames(sce_seurat.list[[i]] )
  }
  sce_reference.list = sce_seurat.list[samples_ref]

  # integrate reference
  ref_anchors <- FindIntegrationAnchors(object.list = sce_reference.list, anchor.features = nrow(sce), dims = 1:d)
  sce_reference.integrated <- IntegrateData(anchorset = ref_anchors, dims = 1:d)
  DefaultAssay(sce_reference.integrated) <- "integrated"
  sce_reference.integrated <- ScaleData(sce_reference.integrated, verbose = FALSE)
  sce_reference.integrated <- RunPCA(sce_reference.integrated, npcs = d, verbose = FALSE)

  # map query
  pca_proj_query = lapply(samples_query , function(current.sample){
    current.sce = sce_seurat.list[[which(names(sce_seurat.list) == current.sample)]]
    query_anchors <- FindTransferAnchors(reference = sce_reference.integrated, query = current.sce,
                                         dims = 1:d, reference.reduction = "pca")
    predictions <- TransferData(anchorset = query_anchors, refdata = t(Embeddings(sce_reference.integrated[['pca']])), dims = 1:d)
    predictions = predictions[1:nrow(predictions), 1:ncol(predictions)]
    current.sce[["pca"]] <- CreateDimReducObject(embeddings = t(as.matrix(predictions)), key = "PC_", assay = DefaultAssay(current.sce))
    current.sce = as.SingleCellExperiment(current.sce)
    out = reducedDim(current.sce , "PCA")
    colnames(out) = paste0("PC_" , c(1:d))
    return(out)
  })

  # combine
  pca_proj_query = do.call(rbind , pca_proj_query)
  sce_reference.integrated = as.SingleCellExperiment(sce_reference.integrated)
  pca_ref = reducedDim(sce_reference.integrated , "PCA")
  colnames(pca_ref) = paste0("PC_" , c(1:d))
  pca_joint = rbind(pca_proj_query , pca_ref)
  pca_joint = pca_joint[order(rownames(pca_joint)) , ]
  sce = sce[, order(colnames(sce))]
  reducedDim(sce , reducedDim.name) = pca_joint
  return(sce)
}




