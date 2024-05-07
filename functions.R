norm_function <- function(sce){
  
  clusters <- quickCluster(sce, use.ranks = FALSE, min.size=50)
  sce <- computeSumFactors(sce, cluster=clusters)
  
  rm(clusters)
  
  sce <- logNormCounts(sce)
  
  rownames(sce@assays@data$logcounts) = rownames(sce)
  colnames(sce@assays@data$logcounts) = colnames(sce)
  
  return(sce)
}

hvg_function <- function(sce, n=1000){
  
  dec <- modelGeneVar(sce)
  
  plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
  curve(metadata(dec)$trend(x), col="blue", add=TRUE)
  
  hvg <- getTopHVGs(dec, n=n)
  print(length(hvg))
  
  return(hvg)
}

make_seurat_function <- function(sce, hvg){
  
  if(is.null(reducedDims(sce)) == FALSE){
    reducedDims(sce) = NULL
  }
  
  Seurat <- as.Seurat(sce,
                      counts = "counts",
                      data = "logcounts")
  VariableFeatures(Seurat) <- hvg
  
  all.genes <- rownames(Seurat)
  Seurat <- ScaleData(Seurat, features = all.genes)
  
  Seurat <- RunPCA(Seurat,
                   features = VariableFeatures(Seurat),
                   reduction.key = "pca_",
                   verbose = FALSE)
  print(plot(Seurat@reductions$pca@stdev,
             xlab = "PC",
             ylab = "Eigenvalue"))
  
  return(Seurat)
}
seurat_cluster_function <- function(Seurat, PCA=15, resolution=0.8, k.param=20, reduction="pca", assay="originalexp"){
  
  Seurat <- FindNeighbors(Seurat, dims=1:PCA, k.param=k.param, reduction=reduction, assay=assay,
                          features = VariableFeatures(Seurat))
  Seurat <- FindClusters(Seurat,
                         resolution = resolution)
  
  Seurat <- RunTSNE(Seurat,
                    dims = 1:PCA,
                    assay=assay,
                    features = VariableFeatures(Seurat),
                    check_duplicates = FALSE,
                    reduction=reduction)
  Seurat <- RunUMAP(Seurat,
                    dims = 1:PCA,
                    assay=assay,
                    reduction=reduction)
  
  return(Seurat)
}

