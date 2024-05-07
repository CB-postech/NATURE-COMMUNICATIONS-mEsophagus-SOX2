library(scran)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(scran)
library(plyr)
library(dplyr)
library(palettetown)

source('functions.R')

total_counts = rbind.fill.matrix(t(s20_sce.qc@assays@data$counts), 
                                 t(s21_sce.qc@assays@data$counts))
rownames(total_counts) = c(colnames(s20_sce.qc), colnames(s21_sce.qc))
total_counts = t(total_counts)
total_counts[is.na(total_counts)] = 0
total_metadata = rbind.fill(as.data.frame(s20_sce.qc@colData), 
                            as.data.frame(s21_sce.qc@colData))
r2flpe_sce <- SingleCellExperiment(assays = list(counts = total_counts),
                                      colData = total_metadata)
head(r2flpe_sce@colData)

ERCC = rownames(r2flpe_sce)[grep("ERCC", rownames(r2flpe_sce))]
remain = rownames(r2flpe_sce)[!rownames(r2flpe_sce) %in% ERCC]
r2flpe_sce = r2flpe_sce[remain, ]

r2flpe_sce.norm <- norm_function(r2flpe_sce)

r2flpe_hvg.norm <- hvg_function(r2flpe_sce.norm)