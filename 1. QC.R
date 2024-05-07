library(SingleCellExperiment)

# load file and check cell barcodes
cellbarcodes = read.table("Cel-seq2_cell_barcode_whitelist.tsv", sep="\t")
table(cellbarcodes$V1 %in% colnames(count$umicount$exon$all))

# s20
count.2 = readRDS("s20.dgecounts.rds")

s20.count <- count.2$umicount$exon$all
s20.count <- s20.count[, cellbarcodes$V1]
colnames(s20.count) = sprintf("s20_%03d", 1:384)

s20_sce <- SingleCellExperiment(assays = list(counts = s20.count))

s20_color = c(rep("RFP", 48), rep("YFP", 48), rep("RFP", 48), rep("YFP", 48), rep("Empty", 9), rep("RFP", 15), rep("Empty", 168))
s20_sce$color = s20_color
s20_sce$plate = "Plate1"
s20_sce$sample = c(rep("Sample1", 24*8), rep("Empty", 9), rep("Sample1", 15), rep("Empty", 24*7))

# s21
count.3 = readRDS("s21.dgecounts.rds")

s21.count <- count.3$umicount$exon$all
s21.count <- s21.count[, cellbarcodes$V1]
colnames(s21.count) = sprintf("s21_%03d", 1:384)

s21_sce <- SingleCellExperiment(assays = list(counts = s21.count))

s21_color = c(rep("RFP", 48), rep("YFP", 48), rep("RFP", 48), rep("YFP", 48), rep("RFP", 48), rep("YFP", 48), rep("RFP", 48), rep("YFP", 22), rep("Empty", 2), rep("YFP", 22), rep("Empty", 2))
s21_sce$color = s21_color
s21_sce$plate = "Plate2"
s21_sce$sample = c(rep("Sample2", 24*6), rep("Sample3", 24*7), rep("Sample1", 20), rep("Sample3", 4), rep("Sample1", 22), rep("Empty", 2), rep("Sample1", 22), rep("Empty", 2))

# QC
library(scater)
library(RColorBrewer)

qc_hist_function <- function(sce, v=0, column, name, path){
  
  hist_column = sce[[column]]
  
  png(paste0(path, "hist_", column, ".png"))
  hist(
    hist_column,
    breaks = 100,
    xlab = name,
    main="")
  if(v != 0){
    abline(v = v, col="red")
  }
  dev.off()
}
qc_dimplot_function <- function(df, column, name, path){
  
  df = df[order(df[,column]),]
  d <- ggplot(df) +
    geom_point(aes(x=PC1, y=PC2, color = eval(parse(text = column)))) +
    scale_colour_gradientn(colours = brewer.pal(9, "Reds"),
                           limits = c(min(df[,column]), max(df[,column]))) +
    ggtitle(name) +
    theme_bw() +
    theme(plot.title = element_text(size = 15, face="italic", hjust = 0.5, margin = margin(0,0,4,0,"mm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size=15, face="italic"),
          legend.title = element_blank(),
          legend.key = element_blank())
  ggsave(paste0(path, "dimplot_", column, ".png"), plot=d, width=6, height=5)
  
}
qc_function <- function(sce, total_counts_cutoff=0, total_features_cutoff=0, pct_counts_ERCC_cutoff=0,
                        save = FALSE, objectname=NULL){
  
  if(sum(is.null(objectname)) == 1){
    objectname = deparse(substitute(sce))
    objectname = gsub("raw", "", objectname)
  }
  
  path = paste0(plotdir, strsplit(objectname, "_")[[1]][1], "_qc/")
  if(sum(dir.exists(path)) == 0){
    dir.create(path)
  }
  
  is.ERCC = grepl("ERCC", rownames(sce))
  print(sum(is.ERCC))
  
  colData <- perCellQCMetrics(sce,
                              subsets = list(ERCC = is.ERCC))
  sce@colData = cbind(sce@colData, colData)
  sce$log10_sum = log10(sce$sum + 1)
  sce$log10_detected = log10(sce$detected + 1)
  
  ### plot function
  
  coldata = c("log10_sum", "log10_detected", "subsets_ERCC_percent")
  columnName = c("log10 total umi counts", "log10 total feature counts", "pct counts ERCC")
  cutoffs = c(total_counts_cutoff, total_features_cutoff, pct_counts_ERCC_cutoff)
  
  # run pca
  # vector <-  c(unique(colnames(sce@colData)))
  print(colnames(sce@colData))
  vector <- c("sum", "detected", "subsets_ERCC_percent")
  sce <- runColDataPCA(sce, ncomponents=5, variables = vector)
  
  # hist
  for(i in 1:length(coldata)){
    qc_hist_function(sce, v=cutoffs[i], coldata[i], columnName[i], path=path)
  }
  
  df <- as.data.frame(sce@int_colData$reducedDims$PCA_coldata)
  df$log10_sum <- sce$log10_sum
  df$log10_detected <- sce$log10_detected
  df$subsets_ERCC_percent <- sce$subsets_ERCC_percent
  
  # pca plot
  for(i in 1:length(coldata)){
    qc_dimplot_function(df, coldata[i], columnName[i], path=path)
  }
  
  filter_by_total_counts = sce$log10_sum > total_counts_cutoff
  filter_by_feature_counts = sce$log10_detected > total_features_cutoff
  filter_by_pct_counts_ERCC = sce$subsets_ERCC_percent < pct_counts_ERCC_cutoff
  
  sce$use <- (
    filter_by_feature_counts &
      filter_by_total_counts &
      filter_by_pct_counts_ERCC
  )
  print(table(sce$use))
  
  ### save function
  if(sum(save) == 1){
    
    sce.qc <- sce[, sce$use]
    
    p <- plotReducedDim(sce, "PCA_coldata",
                        colour_by = "use",
                        size_by = "detected") +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            axis.title = element_blank(),
            panel.border = element_blank(),
            legend.text=element_text(size=13),
            legend.key=element_blank(),
            axis.text = element_blank())
    ggsave(paste0(path, "dimplot_use.png"), plot=p, width=6, height=5)
    
    print(paste0("sce.qc is generated. There are ", ncol(sce.qc), " cells finally."))
    return(sce.qc)
    
  }else{
    return(sce)
  }
}

objectname = "s20"
s20_sce.qc = s20_sce
s20_sce <- qc_function(s20_sce, total_counts_cutoff = 2.5, pct_counts_ERCC_cutoff = 50, save=FALSE, objectname = objectname)
s20_sce.qc <- qc_function(s20_sce.qc, total_counts_cutoff = 2.5, pct_counts_ERCC_cutoff = 50, save=TRUE, objectname = objectname)
non.empty.cells = colnames(s20_sce[, s20_sce$color != "Empty"])
s20_sce.qc = s20_sce.qc[, colnames(s20_sce.qc) %in% non.empty.cells]

objectname = "s21"
s21_sce.qc = s21_sce
s21_sce <- qc_function(s21_sce, total_counts_cutoff = 2.5, pct_counts_ERCC_cutoff = 50, save=FALSE, objectname = objectname)
s21_sce.qc <- qc_function(s21_sce.qc, total_counts_cutoff = 2.5, pct_counts_ERCC_cutoff = 50, save=TRUE, objectname = objectname)
non.empty.cells = colnames(s21_sce[, s21_sce$color != "Empty"])
s21_sce.qc = s21_sce.qc[, colnames(s21_sce.qc) %in% non.empty.cells]
