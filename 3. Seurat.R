library(scran)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(scran)
library(plyr)
library(dplyr)
library(palettetown)

source('functions.R')

#
r2flpe_seurat <- make_seurat_function(r2flpe_sce.norm, r2flpe_hvg.norm)
r2flpe_seurat <- seurat_cluster_function(r2flpe_seurat, PCA=50, resolution = 2.0)

#
r2flpe_seurat.2 <- subset(r2flpe_seurat, idents = c(6), invert=TRUE)

r2flpe_sce.2 <- as.SingleCellExperiment(r2flpe_seurat.2)
r2flpe_sce.norm.2 <- norm_function(r2flpe_sce.2)
r2flpe_hvg.norm.2 <- hvg_function(r2flpe_sce.norm.2)

r2flpe_seurat.2 <- make_seurat_function(r2flpe_sce.norm.2, r2flpe_hvg.norm.2[1:500])
r2flpe_seurat.2 <- seurat_cluster_function(r2flpe_seurat.2, PCA=50, resolution=2.0, k.param = 10)

#
r2flpe_seurat.2 <- CellCycleScoring(r2flpe_seurat.2, s.genes, g2m.genes)

#
r2flpe_seurat.2$celltype = NA

Basal = c(7,8,5,10,0,6,4,2)
Suprabasal = c(1,3,9)

types = c("Basal", "Suprabasal")
for(t in types){
  cl = eval(parse(text = t))
  cells = r2flpe_seurat.2[, r2flpe_seurat.2$seurat_clusters %in% cl] %>% colnames()
  r2flpe_seurat.2$celltype[cells] = rep(t, length(cells))
}

UMAPPlot(r2flpe_seurat.2, group.by = "color", label=TRUE, label.size=8)

# celltype2
r2flpe_seurat.2$celltype2 = NA

Basal.Prolif = c(7,8,5,10)
Basal.Quies = c(0,6,4,2)
Suprabasal = c(1,3,9)

types = c("Basal.Prolif", "Basal.Quies", "Suprabasal")
for(t in types){
  cl = eval(parse(text = t))
  cells = r2flpe_seurat.2[, r2flpe_seurat.2$seurat_clusters %in% cl] %>% colnames()
  r2flpe_seurat.2$celltype2[cells] = rep(t, length(cells))
}

UMAPPlot(r2flpe_seurat.2, group.by = "celltype2", label=TRUE, label.size=8)
