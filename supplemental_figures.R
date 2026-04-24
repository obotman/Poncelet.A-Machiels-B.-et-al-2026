rm(list=ls())
# -*-
library(future)
library(tidyverse)
library(dplyr)
library(readr)
library(readxl)
library(msigdbr)
library(ggrepel)
library(scCustomize)
# -*-
library(Seurat)
library(SeuratExtend)
library(AUCell)
library(SCENIC)
library(fgsea)
library(SCPA)
library(dittoSeq)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(ggVennDiagram)
library(ComplexHeatmap)
library(rstatix)
library(speckle)
library(scuttle)
library(muscat)

plan("multisession", workers = 1)
options(future.globals.maxSize= 36 * 1024 * 1024 * 1024)

# FUN
printOverlapGenes4 <- function(x, filename){
  ov1 <- calculate.overlap(x)
  ov1mx <- max(sapply(ov1, length))
  ov1Mat <- as_tibble(matrix(nrow = ov1mx, ncol = 15))
  ov1Fields <- c("a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10", "a11", "a12", "a13", "a14", "a15")
  colnames(ov1Mat) <- ov1Fields
  for(f1 in ov1Fields){
    if(length(ov1[[f1]]) > 0){
      ov1Mat[1:length(ov1[[f1]]), f1] <- ov1[[f1]]
    }
  }
  write_csv2(as_tibble(ov1Mat), filename, na = "")
}
# FUN

# FUN
get_deg_names <- function(clus, condition, day, direction) {
  # Reconstruct the original list name (e.g., "clus2_MuHV4_d8")
  list_name <- paste0(clus, "_", condition, "_", day)
  
  # Get the DEG table for this comparison
  deg_table <- degs_filtered[[list_name]] # change "degs_filtered" with the good table if necessary
  
  # Filter genes based on direction (up/down)
  if (direction == "Upregulated") {
    genes <- rownames(deg_table[deg_table$avg_log2FC > 0, ])
  } else {
    genes <- rownames(deg_table[deg_table$avg_log2FC < 0, ])
  }
  
  # Collapse into a comma-separated string (or use another separator)
  paste(genes, collapse = ", ")
}
# FUN

# FUN
get_names_conditions_day <- function(condition, day, direction) { # clus,
  # Reconstruct the original list name (e.g., "clus2_MuHV4_d8")
  list_name <- paste0(condition, "_", day) # clus, "_", 
  
  # Get the DEG table for this comparison
  deg_table <- degs_filtered[[list_name]] # change "degs_filtered" with the good table if necessary
  
  # Filter genes based on direction (up/down)
  if (direction == "Upregulated") {
    genes <- rownames(deg_table[deg_table$avg_log2FC > 0, ])
  } else {
    genes <- rownames(deg_table[deg_table$avg_log2FC < 0, ])
  }
  
  # Collapse into a comma-separated string (or use another separator)
  paste(genes, collapse = ", ")
}
# FUN

# FUN
waterfall_plot <- function (fsgea_results, graph_title) {
  fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>% # removes 'HALLMARK_' from the pathway title 
    ggplot( aes(reorder(short_name,NES), NES)) +
      geom_bar(stat= "identity", aes(fill = padj < 0.05))+
      coord_flip()+
      labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
      theme(axis.text.y = element_text(size = 7), 
            plot.title = element_text(hjust = 1))
}
# FUN

# FUN
volcanoCustom <- function(degs, main, submain = NULL){
	p <- ggplot(degs, aes(x = avg_log2FC, y = -log10(plot_pval))) + # -log10(plot_pval)
	geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
	geom_point(data = subset(degs, !is_highlight),aes(color = signif),size = 2, alpha = 0.6) +
	geom_point(data = subset(degs, is_highlight),color = "black", size = 2.5) +     
	geom_text_repel(data = subset(degs, is_highlight),aes(label = label_gene),size = 4.5,color = "black",
					box.padding = 0.5,point.padding = 0.3,segment.color = "black",segment.size = 0.3, max.overlaps = Inf) +
	labs(title = main, x = "Log2 fold change", y = "-log10 adjusted p-value", color = "Significance") +
	scale_color_manual(values = c("Significant" = "lightblue", "NS" = "red")) +
	scale_x_continuous(limits = c(-max(abs(degs$avg_log2FC)),max(abs(degs$avg_log2FC)))) +
	scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000)) +
	# scale_y_sqrt() +
	theme_minimal() +
	annotation_logticks(sides = "l", outside = FALSE)
	  
	if (!is.null(submain)) {
		p <- p + labs(subtitle = submain)
	 }
	  
	return(p)
}
# FUN

# FUN
printOverlapGenes4 <- function(x, filename){
  ov1 <- calculate.overlap(x)
  ov1mx <- max(sapply(ov1, length))
  ov1Mat <- as_tibble(matrix(nrow = ov1mx, ncol = 15))
  ov1Fields <- c("a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10", "a11", "a12", "a13", "a14", "a15")
  colnames(ov1Mat) <- ov1Fields
  for(f1 in ov1Fields){
    if(length(ov1[[f1]]) > 0){
      ov1Mat[1:length(ov1[[f1]]), f1] <- ov1[[f1]]
    }
  }
  write_csv2(as_tibble(ov1Mat), filename, na = "")
}
# FUN

# Color palette
brewer_palette <- "RdBu"
ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_palette))
mr <- ramp(256)[256:1]

# ************
# SUPPLEMENTAL
# ************
setwd("~/Projects/Arthur/scRNAseq_Analysis3/Article/Figures/Supplemental/")

# *********************
# - fig S2 B, C, D, E -
# *********************

# # QC fig
so <- readRDS("../../../RawData/so.rds")
# soList <- readRDS("../../RawData/soList.rds")

expName <- c("d8","d38","d120")

for (name in expName){
  so[[name]]$log10GenesPerUMI = log10(so[[name]]$nFeature_RNA) / log10(so[[name]]$nCount_RNA) # Calculate the complexity (amount genes detected per UMI)
  # so[[name]] <- Add_Mito_Ribo(so[[name]], mito_name = "percent.mt", ribo_name = "percent.rb", species="mouse")#
  so[[name]][["percent.mt"]] <- PercentageFeatureSet(so[[name]], pattern = "^mt.")
}


# **********
# QC figures
# **********

for (name in expName){
  
  p <- VlnPlot(so[[name]], features = "nFeature_RNA", group.by = "Sample_Name", pt.size=0.1, alpha=0.25, raster = FALSE) + 
          labs(x="SampleID", title="Nb of genes detected per cell")+
          geom_hline(yintercept = 6500, color = "red", linetype= "dashed", size=0.75) +
		  geom_hline(yintercept = 1000, color = "red", linetype= "dashed", size=0.75)
          # geom_text(aes(2, 14, label = "6500 genes"), color="red", size=5)
	print(p)
  # {png(paste0("S2_QcNfeature_",name,".png"), width = 6, height = 5, units = "in", res = 300);print(p);dev.off()}
  {pdf(paste0("S2_QcNfeature_",name,".pdf"), width = 6, height = 5);print(p);dev.off()}
  
  
  p <- VlnPlot(so[[name]], features = "nCount_RNA", group.by = "Sample_Name", pt.size=0.00000000001, raster = FALSE) + 
          labs(x="SampleID", title="No. of UMIs per cell")
  print(p)
  # {png(paste0("S2_QcNcount_",name,".png"), width = 6, height = 5, units = "in", res = 300);print(p);dev.off()}
  {pdf(paste0("S2_QcNcount_",name,".pdf"), width = 6, height = 5);print(p);dev.off()}
  
  p <- VlnPlot(so[[name]], features = "percent.mt", group.by = "Sample_Name", pt.size=0.1, alpha = 0.25, raster = FALSE) + 
          labs(x="SampleID", title="% of mitochondrial UMIs per cell") +
          geom_hline(yintercept = 10, color = "red", linetype= "dashed", size=0.75) +
          # geom_text(aes(2, 14, label = "10%"), color="red", size=5) + 
		  NoLegend()
  print(p)
  # {png(paste0("S2_QcPercentMt_",name,".png"), width = 6, height = 5, units = "in", res = 300);print(p);dev.off()}
  {pdf(paste0("S2_QcPercentMt_",name,".pdf"), width = 6, height = 5);print(p);dev.off()}

  
  p <- FeatureScatter(so[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt", raster = FALSE)  + 
          labs(x="No. of UMIs per cell", y="% of mitochondrial UMIs per cell", title=NULL) + NoLegend() +
          geom_hline(yintercept = 10, color = "red", linetype= "dashed", size=1) &
          geom_text(aes(150000, 14, label = "10%"), color="red", size=5)
  print(p)
  # {png(paste0("S2_QcPercentMtVsNcount_",name,".png"), width = 6, height = 5, units = "in", res = 300);print(p);dev.off()}
  {pdf(paste0("S2_QcPercentMtVsNcount_",name,".pdf"), width = 6, height = 5);print(p);dev.off()}

  p <-FeatureScatter(so[[name]], feature1 = "nFeature_RNA", feature2 = "percent.mt", raster = FALSE) + 
          labs(x="No. of genes detected per cell", y="% of mitochondrial UMIs per cell", title=NULL) + NoLegend()
  print(p)
  # {png(paste0("S2_QcPercentMtVsNfeature_",name,".png"), width = 6, height = 5, units = "in", res = 300);print(p);dev.off()}
  {pdf(paste0("S2_QcPercentMtVsNfeature_",name,".pdf"), width = 6, height = 5);print(p);dev.off()}

  p <- FeatureScatter(so[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE) + 
          labs(x="No. of UMIs per cell", y="No. of genes detected per cell", title=NULL) + NoLegend()
  print(p)
  # {png(paste0("S2_QcNfeatureVsNcount_",name,".png"), width = 6, height = 5, units = "in", res = 300);print(p);dev.off()}
  {pdf(paste0("S2_QcNfeatureVsNcount_",name,".pdf"), width = 6, height = 5);print(p);dev.off()}
  
  # png(paste0("S2_QcNfeatureHist_",name,".png"), width = 6, height = 5, units = "in", res = 300)
  pdf(paste0("S2_QcNfeatureHist_",name,".pdf"), width = 6, height = 5)
  print(hist(so[[name]]@meta.data$`nFeature_RNA`, 50, xlab = "No. of features", ylab = "Frequency", main = "Histogram of features/cell"))
  dev.off()
  
  # png(paste0("S2_QcNmtHist_",name,".png"), width = 6, height = 5, units = "in", res = 300)
  pdf(paste0("S2_QcNmtHist_",name,".pdf"), width = 6, height = 5)
  print(hist(so[[name]]@meta.data$`percent.mt`,50, xlab = "Percent of mitochondrial genes", ylab = "Frequency", main = "Histogram of mitochondrial genes/cell"))
  dev.off()
  
  # https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

  
  p <- so[[name]]@meta.data %>% ggplot(aes(x=log10GenesPerUMI, color = Sample_Name, fill=Sample_Name)) +
		geom_density(alpha = 0.2) + theme_classic() + geom_vline(xintercept = 0.8)
  print(p)
  # {png(paste0("S2_QcComplexicity_",name,".png"), width = 6, height = 5, units = "in", res = 300);print(p);dev.off}
  {pdf(paste0("S2_QcComplexicity_",name,".pdf"), width = 6, height = 5);print(p);dev.off}
	
  p <- so[[name]]@meta.data %>% ggplot(aes(x=Sample_Name, y=log10(nFeature_RNA), fill=Sample_Name)) + 
		geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
		theme(plot.title = element_text(hjust=0.5, face="bold")) + ggtitle("NCells vs NGenes")
  print(p)
  # {png(paste0("S2_QcBoxPlotNrna_",name,".png"), width = 6, height = 5, units = "in", res = 300);print(p);dev.off()}
  {pdf(paste0("S2_QcBoxPlotNrna_",name,".pdf"), width = 6, height = 5);print(p);dev.off()}
}

# ************
# Filter cells
# ************

for (name in expName){
  minFeatureRna <- 1000
  maxFeatureRna <- 6500
  maxPercentMt <- 10
  so[[name]] <- subset(so[[name]], subset = (nFeature_RNA > minFeatureRna) & (nFeature_RNA < maxFeatureRna) & (percent.mt < maxPercentMt))
}

# *******************
# Merge seurat Object
# *******************

seuObj <- Merge_Seurat_List(so, add.cell.ids = expName)
rm(so)
gc()

# *************************************
# # Density plot and correlation matrix
# *************************************
library(SingleCellExperiment)
library(muscat)
library(DESeq2)

seuObj <- JoinLayers(seuObj)
sce <- as.SingleCellExperiment(seuObj)
sce$Sample_Name <- factor(sce$Sample_Name, levels = c("Mock_r1_d8","Mock_r2_d8",
													"MuHV4_r1_d8","MuHV4_r2_d8","MuHV4_r3_d8",
													"del73_r1_d8","del73_r2_d8","del73_r3_d8",
													"Mock_r1_d38","Mock_r2_d38","Mock_r3_d38",
													"MuHV4_r1_d38","MuHV4_r2_d38","MuHV4_r3_d38",
													"del73_r1_d38","del73_r2_d38","del73_r3_d38",
													"Mock_r1_d120","Mock_r2_d120","Mock_r3_d120",
													"MuHV4_r1_d120","MuHV4_r2_d120","MuHV4_r3_d120",
													"del73_r1_d120","del73_r2_d120","del73_r3_d120"))

rm(seuObj)

colnames(sce) = paste0("cell_",1:ncol(sce))
assayNames(sce) <- "counts"
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
qc <- perCellQCMetrics(sce)
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

sceAgr <- aggregateData(sce,
    assay = "counts", fun = "sum",
    by = "Sample_Name")

assayNames(sceAgr) <- "counts"

countMatrix <- data.frame(as.matrix(counts(sceAgr)))
countMatrix$Genes <- rownames(countMatrix)
countMatrix <- distinct(countMatrix, Genes, .keep_all=TRUE)
rownames(countMatrix) <- countMatrix$Genes
countMatrix$Genes <- NULL

epsilon <- 1 # pseudo-count to avoid problems with log(0)
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMax(countMatrix)

all_counts <- as.numeric(as.matrix(countMatrix))
# png("S2_histDistributionAggr.png")
pdf("S2_histDistributionAggr.pdf")
hist(log2(all_counts + epsilon), breaks=100, col="blue", border="white",
     main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)
dev.off()

# png("S2_sampleDistribution.png",width=10, height=25, units = "in", res = 300)
pdf("S2_sampleDistribution.pdf",width=10, height=25)
boxplot(log2(countMatrix + epsilon), pch=".", 
        horizontal=TRUE, cex.axis=0.5,
        las=1, ylab="Samples", xlab="log2(Counts +1)")
dev.off()

Mo <- grepl("^Mock*",colnames(sceAgr))%>%sum()
del <- grepl("^del73*",colnames(sceAgr))%>%sum() 
Mu <- grepl("^MuHV4*",colnames(sceAgr))%>%sum() 

condition <- c(rep("Mock", 2), rep("Del73", 3), rep("MuHV4", 3),
				rep("Mock", 3), rep("Del73", 3), rep("MuHV4", 3),
				rep("Mock", 3), rep("Del73", 3), rep("MuHV4", 3))
my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(sceAgr)
my_colData
my_colData$group <- rownames(my_colData)

d8 <- grep("*d8", my_colData$group)
for (i in d8){my_colData[i,]$group <- "d8"}
d38 <- grep("*d38", my_colData$group)
for (i in d38){my_colData[i,]$group <- "d38"}
d120 <- grep("*d120", my_colData$group)
for (i in d120){my_colData[i,]$group <- "d120"}

dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = my_colData,
                              design = ~ condition + group + condition:group)
# design(dds) <- ~ condition + group + condition:group
dds <- DESeq(dds)
resultsNames(dds)
normalized_counts <- counts(dds, normalized = T)
head(normalized_counts)
vsd <- vst(dds, blind = TRUE)

library(affy)
library(scales)
myColors <- hue_pal()(4)

# png("S2_densityPlot.png",width=10, height=14, units = "in", res = 300)
pdf("S2_densityPlot.pdf",width=10, height=14)
plotDensity(log2(normalized_counts + 0.1), col=rep(myColors, each=3),
            lty=c(1:ncol(normalized_counts)), xlab='Log2(normalized_counts)',
            main='Expression Distribution')
legend('topright', colnames(normalized_counts), lty=c(1:ncol(normalized_counts)),col=rep(myColors, each=3))
abline(v=-1.5, lwd=1, col='red', lty=2)
dev.off()

count_melt <- reshape2::melt(log2(normalized_counts + epsilon))
head(count_melt)
# png("S2_densityPlot2.png",width=10, height=14, units = "in", res = 300)
pdf("S2_densityPlot2.pdf",width=10, height=14)
ggplot(count_melt, aes(x=value, color = Var2)) + geom_density() #mapping=aes(x=value, color=Var2))
dev.off()

plotFun <- function(x,y){ 
  dns <- densCols(x,y); 
  points(x,y, col=dns, pch=".", panel.first=grid());  
#  abline(a=0, b=1, col="brown")
  }
set.seed(123) # forces the random number generator to produce fixed results. Should generally not be used, except for the sake of demonstration with a particular selection. 
nb.pairs <- 26 # too much
# png("S2_scatter2.png", width=26, height = 26,units = "in", res = 300)
pdf("S2_scatter2.pdf", width=26, height = 26)
# pairs(log2(countMatrix[,sample(ncol(countMatrix), nb.pairs)] + epsilon), 
      # panel=plotFun, lower.panel = NULL)
pairs(log2(countMatrix[,seq(1:26)] + epsilon), 
      panel=plotFun, lower.panel = NULL)
dev.off()

vsd.obj=vsd
plotDists = function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd.obj$condition)
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
			cluster_rows = FALSE,
			cluster_cols = FALSE,
			clustering_distance_rows = sampleDists,
			clustering_distance_cols = sampleDists,
			main = "Correlation based on normalized expression values",
			col = colors)
}
p <- plotDists(vsd)
print(p)
# {png(paste0("S2_QC_correlationPlot.png"), width = 10, height = 10, units = "in", res = 300);print(p);dev.off()}
{pdf(paste0("S2_QC_correlationPlot.pdf"), width = 10, height = 10);print(p);dev.off()}


plotDists = function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd.obj$condition)
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
			cluster_rows = TRUE,
			cluster_cols = TRUE,
			clustering_distance_rows = sampleDists,
			clustering_distance_cols = sampleDists,
			main = "Correlation based on normalized expression values",
			col = colors)
}
p <- plotDists(vsd)
print(p)
# {png(paste0("S2_QC_correlationPlot_clustered.png"), width = 10, height = 10, units = "in", res = 300);print(p);dev.off()}
{pdf(paste0("S2_QC_correlationPlot_clustered.pdf"), width = 10, height = 10);print(p);dev.off()}

plot_PCA = function(vsd.obj){
    pcaData <- plotPCA(vsd.obj,  intgroup = c("condition"), returnData = T, pcsToUse = c(1,2))
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes(PC1, PC2, color=condition)) +
        geom_point(size=3) +
        labs(x = paste0("PC1: ",percentVar[1],"% variance"),
             y = paste0("PC2: ",percentVar[2],"% variance"),
             title = "PCA Plot colored by condition") +
        ggrepel::geom_text_repel(aes(label = name), color = "black")
}
p <- plot_PCA(vsd)
print(p)
# {png(paste0("S2_QC_pca1_2.png"), width = 8, height = 6, units = "in", res = 300);print(p);dev.off()}
{pdf(paste0("S2_QC_pca1_2.pdf"), width = 8, height = 6);print(p);dev.off()}

plot_PCA = function(vsd.obj){
    pcaData <- plotPCA(vsd.obj,  intgroup = c("condition"), returnData = T, pcsToUse = c(1,3))
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes(PC1, PC3, color=condition)) +
        geom_point(size=3) +
        labs(x = paste0("PC1: ",percentVar[1],"% variance"),
             y = paste0("PC3: ",percentVar[2],"% variance"),
             title = "PCA Plot colored by condition") +
        ggrepel::geom_text_repel(aes(label = name), color = "black")
}
p <- plot_PCA(vsd)
print(p)
# {png(paste0("S2_QC_pca1_3.png"), width = 8, height = 6, units = "in", res = 300);print(p);dev.off()}
{pdf(paste0("S2_QC_pca1_3.pdf"), width = 8, height = 6);print(p);dev.off()}

# ************
# - fig S2 F -
# ************
setwd("~/Projects/Arthur/scRNAseq_Analysis3/Article/Figures/Supplemental/")
so <- readRDS("../../../RawData/seuObjRmCycle.rds")

cluster_to_color = c(
'0'='#E69F00',
'1'='#56B4E9',
'2'='#009E73',
'3'='#F0E442',
'4'='#0072B2',
'5'='#D55E00',
'6'='#CC79A7',
'7'='#666666',
'8'='#AD7700',
'9'='#1C91D4',
'10'='#007756',
'11'='#D5C711',
'12'='#005685',
'13'='#A04700',
'14'='#B14380',
'15'='#4D4D4D',
'16'='#FFBE2D'
)

p <- DimPlot(so, reduction = "umap_SCT", group.by = "SCT_clusters_0.6", cols = cluster_to_color, label = TRUE)
# {png("S2F_umap_SCT_cluster_0.6.png", width = 10, height = 8, units = "in", res = 300);print(p);dev.off()}
{pdf("S2F_umap_SCT_cluster_0.6.pdf", width = 10, height = 8);print(p);dev.off()}

Idents(so) <- "SCT_clusters_0.6"
mo <- subset(so, ident = c(0,2,4,9,10,11,14,15,16))
hsc <- subset(so, ident = c(1,3,5,6,7,8,12,13))

mo_markers <- FindAllMarkers(mo, group.by = "SCT_clusters_0.6", assay="RNA")
write.csv2(mo_markers, file = "S2F_markers_clusters_mo.csv", row.names = FALSE)

hsc_markers <- FindAllMarkers(hsc, group.by = "SCT_clusters_0.6", assay="RNA")
write.csv2(hsc_markers, file = "S2F_markers_clusters_hsc.csv", row.names = FALSE)

top_mo <- mo_markers %>% dplyr::group_by(cluster) %>%
						dplyr::filter(avg_log2FC > 0, pct.1 > 0.1, p_val_adj < 0.05) %>% 
						dplyr::top_n(10, avg_log2FC)

geneVector <- unique(top_mo$gene)
widthEst <- 2 + length(unique(top_mo$cluster)) * 0.4
heightEst <- 3 + length(geneVector) * 0.15
	
p <- DotPlot(mo, features = geneVector, group.by = "SCT_clusters_0.6", assay = "RNA") + 
			coord_flip() +
			labs(x="Genes", y="Clusters", title = "top markers in MO clusters")
print(p)
# {png("S2F_dotPlot_mo_cluster_marker.png", width = widthEst, height = heightEst, units = "in", res = 300);plot(p);dev.off()}
{pdf("S2F_dotPlot_mo_cluster_marker.pdf", width = widthEst, height = heightEst);plot(p);dev.off()}

top_hsc <- hsc_markers %>% dplyr::group_by(cluster) %>%
						dplyr::filter(avg_log2FC > 0, pct.1 > 0.1, p_val_adj < 0.05) %>% 
						dplyr::top_n(10, avg_log2FC)

geneVector <- unique(top_hsc$gene)
widthEst <- 2 + length(unique(top_hsc$cluster)) * 0.4
heightEst <- 3 + length(geneVector) * 0.15
	
p <- DotPlot(hsc, features = geneVector, group.by = "SCT_clusters_0.6", assay = "RNA") + 
			coord_flip() +
			labs(x="Genes", y="Clusters", title = "top markers in HSC clusters")
print(p)
# {png("S2F_dotPlot_hsc_cluster_marker.png", width = widthEst, height = heightEst, units = "in", res = 300);plot(p);dev.off()}
{pdf("S2F_dotPlot_hsc_cluster_marker.pdf", width = widthEst, height = heightEst);plot(p);dev.off()}

# **********
# - fig S3 -
# **********

# so <- readRDS(file="~/Projects/Arthur/scRNAseq_Analysis3/RawData/seuObjClustersFiltered.rds")
so <- readRDS(file="../so.rds")
DefaultAssay(so) <- "RNA"
Idents(so) <- "SCT_clusters_0.6"

# so$cell_group <- dplyr::recode(so$SCT_clusters_0.6, "0"="MO", "2"="MO","4"="MO","9"="MO","10"="MO","11"="MO",
								# "1"="HSC","2"="HSC","3"="HSC","5"="HSC","6"="HSC","7"="HSC","8"="HSC","12"="HSC","13"="HSC")

# table(so$cell_group)

# - fig S3A -
# volcanoHMpathway_pMO_0+9_MuHV4_del73_d8
# 10-2025_fig_APO.R
Idents(so) <- "SCT_clusters_0.6"
soSub <- subset(so, ident = c("0","9"))
Idents(soSub) <- "day"
soSub <- subset(soSub, ident = "d8")

resList <- list()
for(virus in c("MuHV4","del73")){
	resting <- seurat_extract(soSub,
							  meta1 = "condition", value_meta1 = "Mock")
	activated <- seurat_extract(soSub,
							  meta1 = "condition", value_meta1 = paste0(virus))

	pathways <- msigdbr(species = "Mus musculus", category = "H") %>%
		 format_pathways()
	
	# pathways <- msigdbr(species = "Mus musculus", category = "C5", subcollection = "GO:BP") %>%
		 # format_pathways()

	resList[[virus]] <- compare_pathways(samples = list(activated,resting), 
								 pathways = pathways, parallel = TRUE, cores = 4)
	resList[[virus]] <- resList[[virus]] %>% mutate(condition = case_when( FC > 0 ~ virus,
																		FC < 0 ~ paste0("Mock_",virus)))

}
rest_act <- bind_rows(resList)
save(resList, file = "D_resListHM_cMO_clus2+4+10_d8.Rdata")
# rest_act$log10Pval <- -log10(res_act$adjPval)

rest_act <- rest_act %>%
mutate(color = case_when(FC > 5 & adjPval < 0.05 ~ "Strong enrichment",
						FC < -5 & adjPval < 0.05 ~ "Strong enrichment",
						between(FC, -5, 5) & adjPval < 0.05 ~ "Mild enrichment",
						TRUE ~ "Not significant"),
		shape_group = case_when(
		grepl("MuHV4", condition) ~ "MuHV4",
		grepl("del73", condition) ~ "del73",
		TRUE ~ NA_character_))


aa_path <- rest_act %>% filter(FC > 0.5 & adjPval < 0.05) %>% arrange(desc(FC)) %>% head(,n=5)
# aa_path <- rest_act %>% filter(-FC > 0.5 & adjPval < 0.05) %>% arrange(FC) %>% head(,n=5) %>% rbind(aa_path)
# filter(grepl(pattern = "Hallmark_interferon", ignore.case = T, x = Pathway))
aa_path$Pathway <- gsub("HALLMARK_","",aa_path$Pathway)
aa_path$Pathway <- gsub("_"," ",aa_path$Pathway)
aa_path$Pathway <- str_to_title(aa_path$Pathway)
# filter(grepl(pattern = "Hallmark_interferon", ignore.case = T, x = Pathway))

# Create a separate column to identify the aa_path points
rest_act$Pathway <- gsub("HALLMARK_","",rest_act$Pathway)
rest_act$Pathway <- gsub("_"," ",rest_act$Pathway)
rest_act$Pathway <- str_to_title(rest_act$Pathway)
rest_act <- rest_act %>%
  mutate(is_aa_path = Pathway %in% aa_path$Pathway)

write.csv2(rest_act, file = "../../Table/S3A_HMpathway_results_pMO_0+9_MuHV4_del73_d8.csv")

ggplot(rest_act, aes(x = FC, y = qval, shape = shape_group, fill = color)) +
  scale_shape_manual(values = c("MuHV4" = 21, "del73" = 24), na.translate = FALSE) +
  scale_fill_manual(
    name = "Significance",
    values = c(
      "Strong enrichment" = "orangered2",
      "Mild enrichment" = "mediumseagreen",
      "Not significant" = "black"
    )
  ) +
  geom_vline(xintercept = c(-5, 0, 5), linetype = "dashed", color = c("red", "black", "red"), size = 0.3) +
  geom_point(size = 3, stroke = 0.3) +  # All points with mapped aesthetics
  guides(fill = guide_legend(override.aes = list(shape = 21, color = NA))) +
  # Highlight aa_path points with a border or different styling
  # geom_point(data = subset(rest_act, is_aa_path), 
             # aes(x = FC, y = qval), 
             # shape = 21, size = 3, color = "darkred", stroke = 0.8) +
  ggrepel::geom_text_repel(data=aa_path, aes(label=Pathway), size = 3, max.overlaps = 10, 
							# arrow = arrow(length = unit(0.02, "npc")),segment.color = "black", # Dark gray lines
  segment.size = 0.4,                 # Slightly thinner
  segment.alpha = 0.6,                # More transparent
  segment.linetype = "dashed",        # Dashed lines
  box.padding = 0.8,                  # More padding
  point.padding = 0.5,
  force = 20         ) +

  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1) +
  labs(title = "pMO (0+9): MuHV4 & del73 vs mock at d8", 
  subtitle = "Strong enrichment: FC > 5 & padj < 0.05 \nMild enrichment: −5 < FC < 5 & padj < 0.05 \nNot significant: padj > 0.05")
# ggsave(paste0("D_volcanoPathway_cMO_0+9_MuHV4_del73_d8.png"), width = 7, height = 7, units = "in", dpi = 300)
ggsave(paste0("S3A_volcanoHMpathway_pMO_0+9_MuHV4_del73_d8.pdf"), width = 7, height = 7, units = "in", dpi = 300)

# - fig S3B - 
# Venn cMO d8d38d120 ; pMO d8d38d120 ; d120 cMO vs pMO
Idents(so) <- "SCT_clusters_0.6"
MO <- subset(so, idents = c(0,2,4,9,10))
MO$MO_cp <- dplyr::recode(MO$SCT_clusters_0.6, "0"="pMO", "9"="pMO", "2"="cMO", 
							"4"="cMO", "10"="cMO")
							
MO$condition_day <- dplyr::recode(MO$condition_day,
								"1_Mock_d8"="Mock_d8", "1_Mock_d38"="Mock_d38", "1_Mock_d120"="Mock_d120", 
								"3_del73_d8"="del73_d8", "3_del73_d38"="del73_d38", "3_del73_d120"="del73_d120", 
								"2_MuHV4_d8"="MuHV4_d8", "2_MuHV4_d38"="MuHV4_d38", "2_MuHV4_d120"="MuHV4_d120")

degs <- list()
condition <- c("MuHV4","del73")
for(clus in c("pMO","cMO")){ 
	Idents(MO) <- "MO_cp"
	soSub <- subset(MO, ident = clus)
	for(con in condition){
		for(d in c("d8","d38","d120")){
			Idents(soSub) <- "condition_day"
			deg_key <- paste0(clus, "_", con, "_", d)			
			deg_result <- FindMarkers(soSub, ident.1 = paste0(con, "_", d), ident.2 = paste0("Mock_", d), only.pos = FALSE)
			deg_result$gene <- rownames(deg_result)
			deg_result$group <- deg_key
			degs[[deg_key]] <- deg_result
		}
	}
}

allDEGs <- bind_rows(degs)
write.csv2(allDEGs, file = "../../Table/S3B_allDEGs_cMO_pMO.csv")

# At d120 only MuHV4 pMO vs cMO
d120List <- list(degs$pMO_MuHV4_d120,degs$cMO_MuHV4_d120)

degs_f <- list()
degs_f <- lapply(d120List, function(x) {x %>% filter(avg_log2FC > 0.5 , p_val_adj < 0.05)})
names(degs_f) <- c("pMO_MuHV4_d120","cMO_MuHV4_d120")

lfc = 0.5
gene_lists <- list(degs_f[[1]]$gene, degs_f[[2]]$gene)
names(gene_lists) <- c("pMO_MuHV4_d120","cMO_MuHV4_d120")
totalGene <- calculate.overlap(gene_lists) %>% sapply(., length) %>% sum

p <- ggvenn::ggvenn(gene_lists, auto_scale=FALSE, text_size = 6, fill_color = c("darkblue","lightblue")) + 
ggtitle(paste0("Log2FC>", lfc,", total genes=", totalGene)) #+ 
		# theme(plot.title = element_text(hjust = 0.5)))
print(p)
# {png(file=paste0("pMOcMO_MuHV4_d120_FC_",lfc,".png"), width=7, height=7, units="in", res=300);print(p);dev.off()}
{pdf(file=paste0("S3B_pMOcMO_MuHV4_d120_FC_",lfc,".pdf"), width=7, height=7);print(p);dev.off()}

printOverlapGenes4(gene_lists, file=paste0("../../Table/S3B_pMOcMO_MuHV4_d120_FC",lfc,".csv"))

# - fig S3C -
# pMO d120 volcano MuHV4 vs Mock
Idents(so) <- "SCT_clusters_0.6"
pMO <- subset(so, ident = c(0,9))
Idents(pMO) <- "day"
pMO <- subset(pMO, ident = "d120")
Idents(pMO) <- "condition"
markers <- FindMarkers(pMO, ident.1 = "MuHV4", ident.2 = "Mock")
markers$gene <- rownames(markers)
write.csv2(markers, file = "../../Table/S3C_DEGs_pMO_MuHV4vsMock_d120.csv")
# markers <- read.csv2(file = "../../Table/S3C_DEGs_pMO_MuHV4vsMock_d120.csv")

# significant_genes <- markers %>% filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.05)
# top_up <- significant_genes %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 20)
# top_down <- significant_genes %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC) %>% slice_head(n = 20)
# highlight_genes <- bind_rows(top_up, top_down) %>% pull(gene)

highlight_genes <- markers %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
  mutate(direction = ifelse(avg_log2FC > 0, "positive", "negative")) %>%
  group_by(direction) %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE) %>%
  slice_head(n = 25) %>%
  ungroup() %>% pull(gene)

# highlight_genes <- c("AW112010","H2.Eb1","H2.Aa","H2.Ab1","Iigp1","Cd74","Ciita","Gbp4","Gbp6","Gbp8","Ifit1","Ifit3","Ifit1bL1","Ifit3b",
# "Clec14a","Tm4sf1","Cyp26b1","Ifi208","Slamf7","Mpo","Slpi","Elane","Ms4a3","Clec2a","Olr1","Tmem119","Tmem178")          

# markers <- markers%>%filter(!str_detect(gene, 'MT-'))
markers$label_gene <- ifelse(markers$gene %in% highlight_genes, markers$gene, NA)
markers$is_highlight <- ifelse(markers$gene %in% highlight_genes, TRUE, FALSE)
markers$signif <- with(markers, ifelse(p_val_adj < 0.05, "Significant", "NS")) # & abs(avg_log2FC) > 0.5
markers$plot_pval <- ifelse(markers$p_val_adj == 0, 1e-300, markers$p_val_adj)

p <- volcanoCustom(markers, "scRNAseq: DEG between pMO MuHV4 vs Mock at d120", "clusters 0+9")
print(p)
# {png(paste0("S3C_volcanoCustom_pMO_MuHV4vsMock_d120.png"), width=10, height=10,units="in", res=300);print(p);dev.off()}
{pdf(paste0("S3C_volcanoCustom_pMO_MuHV4vsMock_d120.pdf"), width=10, height=10);print(p);dev.off()}

# - fig S3D -
# Pyhton language
import os
import h5py
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import anndata
import scvelo as scv
import numpy as np
import scipy.sparse
import gc

os.chdir('/mnt/c/Users/u192308/Mes documents/Projects/Arthur/scRNAseq_Analysis3/Trajectories/Velocyto')
outFolder = '/mnt/c/Users/u192308/Mes documents/Projects/Arthur/scRNAseq_Analysis3/Article/Figures/Supplemental/'
ad = sc.read_h5ad('ad_all_clus.h5ad')
umap = pd.read_csv("umap_Full.csv")
meta = pd.read_csv("meta_Full.csv")

MO_clusters = [0, 2, 4, 9, 10, 11]
ad_mo = ad[
    (ad.obs['SCT_clusters_0.6'].astype(int).isin(MO_clusters)) &
    (ad.obsm['X_umap'][:, 0] < -1.5)
].copy()

cluster_to_color = {
'0':'#E69F00',
'2':'#009E73',
'4':'#0072B2',
'9':'#1C91D4',
'10':'#007756',
'11':'#D5C711'
}
clusters = sorted(ad_mo.obs['SCT_clusters_0.6'].unique())
seurat_palette = [cluster_to_color.get(str(cluster), '#000000') for cluster in clusters]

ad_mo.obs['SCT_clusters_0.6'] = ad_mo.obs['SCT_clusters_0.6'].astype('category')
sc.pl.umap(ad_mo, color='SCT_clusters_0.6', size=20, alpha=0.7, frameon=False, legend_loc='on data',
            title='UMAP - Clusters SCT_clusters_0.6', show=True, palette=seurat_palette)

scv.pp.filter_genes(ad_mo, min_shared_counts=20)
scv.pp.normalize_per_cell(ad_mo, enforce=False)
scv.pp.filter_genes_dispersion(ad_mo, n_top_genes=2000)
sc.pp.pca(ad_mo)
sc.pp.neighbors(ad_mo, n_pcs=30, n_neighbors=30)
scv.pp.moments(ad_mo, n_pcs=None, n_neighbors=None)

scv.tl.velocity(ad_mo, mode='stochastic', vkey='velocity_stochastic')
scv.tl.velocity_graph(ad_mo, vkey='velocity_stochastic')
scv.tl.recover_dynamics(ad_mo, n_top_genes=2000)
scv.tl.velocity(ad_mo, mode='dynamical', vkey='velocity_dynamical')  
scv.tl.velocity_graph(ad_mo, vkey='velocity_dynamical')
scv.tl.velocity_embedding(ad_mo, vkey='velocity_stochastic', basis='umap')
scv.tl.velocity_embedding(ad_mo, vkey='velocity_dynamical', basis='umap')

group='SCT_clusters_0.6'
size=10
alpha=0.5
mods=['stochastic','dynamical']
for mode in mods:
    scv.pl.velocity_embedding_stream(ad_mo, basis='umap', color=group, figsize=(10, 10), 
                                title=f'MO, {mode}', save=f"{outFolder}S3D_umap_MO_traj_{mode}.pdf", 
                                size=size, alpha=alpha, dpi=300, palette=seurat_palette, vkey=f"velocity_{mode}")  
        
    scv.pl.velocity_embedding_grid(ad_mo, basis='umap', color=group, title=f'MO, {mode}', 
                                scale=0.25, figsize=(10, 10), save=f"{outFolder}S3D_umap_MO_gridTraj_{mode}.pdf", 
                                size=size, alpha=alpha, dpi=300, palette=seurat_palette, vkey=f"velocity_{mode}")
        
    scv.pl.velocity_embedding(ad_mo, basis='umap', color=group, figsize=(10, 10),arrow_length=3, arrow_size=2,
                            title=f'MO, {mode}', save=f"{outFolder}S3D_umap_MO_embTraj_{mode}.pdf", 
                            size=size, alpha=alpha, dpi=300, palette=seurat_palette, vkey=f"velocity_{mode}")

# - fig V -
# boxPlot_MHCII_low_clus0+2+4+9+10
# 10-2025_fig_APO.R

# cellAssignAutoMHCII_MO_clus0+2+4+9+10
Idents(so) <- "SCT_clusters_0.6"
soSub <- subset(so, ident = c(0,2,4,9,10))

# MHCII scores
cellList <- list()
cellList[["MHCII_scores"]] <- c("H2-ia-ie-AbSeq", "H2.Aa", "H2.Eb1", "H2.Ab1","Ciita","H2.DMb1","AW112010") #"H2.DMb2",

# DefaultAssay(soSub) <- "SCT"
# soSub <- AddModuleScore(soSub, features = cellList, name = "NES_", nbin = 50, ctrl = 100) # bin = 24, ctrl = 50
# for (i in seq_along(cellList)) {
  # colname <- paste0("NES_", names(cellList)[i])
  # soSub@meta.data[[colname]] <- soSub@meta.data[[paste0("NES_", i)]]
  # soSub@meta.data[[paste0("NES_", i)]] <- NULL
# }

expr_matrix <- as.matrix(GetAssayData(soSub, layer = "data")) # Get normalized counts (log-normalized layer)
gene_rankings <- AUCell_buildRankings(expr_matrix) # Build gene-expression rankings per cell
cells_AUC <- AUCell_calcAUC(geneSets = cellList, rankings = gene_rankings) # Calculate AUC scores
auc_scores <- getAUC(cells_AUC)
MHCII_scores <- auc_scores["MHCII_scores",]

# png(paste0("cellAssignAutoMO.png"), width = 6, height = 6, units = "in", res = 300)
pdf(paste0("M_cellAssignAutoMHCII_MO_clus0+2+4+9+10.pdf"), width = 6, height = 6)   
# par(mfrow=c(1,2))
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
dev.off()
soSub$MHCII_AUC_score <- as.numeric(MHCII_scores[colnames(soSub)])

cells_assignment[["MHCII_scores"]][["aucThr"]]
# $selected
# minimumDens 
   # 0.357437 

# $thresholds
                 # threshold nCells
# tenPercentOfMax 0.09664038  11518
# Global_k1       0.34118769   5952
# minimumDens     0.35743704   5803

ncol(soSub)
26144
sum(MHCII_scores > 0.35743704)
5803

threshold = 0.35743704
MHCII_class <- ifelse(MHCII_scores > threshold, "MHCII high", "MHCII low")

soSub$MHCII_AUC_status <- MHCII_class
soSub$MHCII_AUC_score <- as.numeric(MHCII_scores[colnames(soSub)])

soTemp <- subset(soSub, cells = WhichCells(soSub, expression = umapSCT_1 < 0))
soTemp <- subset(soTemp, cells = WhichCells(soTemp, expression = umapSCT_2 < 5)) 
p <- DimPlot(soTemp, group.by = "MHCII_AUC_status", label = FALSE, reduction = "umap_SCT") +
	scale_color_manual(values = c("MHCII high" = "red","MHCII low" = "gray")) +
	 labs(title = "MO MHCII signature without AbSeq integration")
print(p)
# {png(file = paste0("M_umap_MO_MHCII_sig1.png"), width = 6, height= 3, units = "in", res = 300);plot(p);dev.off()}
{pdf(file = paste0("M_umap_MO_MHCII_sig1.pdf"), width = 8, height= 8);plot(p);dev.off()}
rm(soTemp)

# dotPlot_MO_MHCII_sig1
p <- DotPlot(soSub, assay = "RNA", features = cellList[["MHCII_scores"]], group.by="MHCII_AUC_status", scale = TRUE) +
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        legend.box = "horizontal",           # Horizontal arrangement
        legend.box.just = "center",          # Center alignment
        legend.spacing.x = unit(0.5, 'cm')) + # Spacing between legends
    guides(color = guide_colorbar(title = "Avg Expression", order = 1),
        size = guide_legend(title = "Pct Exp", order = 2)) + 
	labs(title = "MO MHCII signature without AbSeq integration") +
	theme(plot.title = element_text(size = 10))
p[["data"]][["features.plot"]][is.na(p[["data"]][["features.plot"]])] <- "H2-ia-ie-AbSeq"
# {png(file = paste0("M_dotPlot_MO_MHCII_sig1.png"), width = 6, height= 3, units = "in", res = 300);plot(p);dev.off()}
{pdf(file = paste0("M_dotPlot_MO_MHCII_sig1.pdf"), width = 6, height= 3);plot(p);dev.off()}

# # Integrate ADT with RNA scores
soSub$MHCII_RNA <- soSub$MHCII_AUC_score

DefaultAssay(soSub) <- "ADT"
adt_expression <- GetAssayData(soSub, layer = "data")["H2-ia-ie-AbSeq", ]
hist(adt_expression)
# quantile(adt_expression, 0.75)
     # 75% 
# 0.805331 
soSub$MHCII_ADT_binary <- adt_expression > median(adt_expression)

# 3. Create final classification
# soSub$MHCII_final_class <- ifelse(
  # soSub$MHCII_RNA > quantile(soSub$MHCII_RNA, 0.75) & soSub$MHCII_ADT_binary,
  # "High Confidence MHC-II+",
  # ifelse(soSub$MHCII_RNA > quantile(soSub$MHCII_RNA, 0.75),
         # "RNA-based MHC-II+",
         # "MHC-II low/negative")
# )

soSub$MHCII_AbSeq_integration <- ifelse(
    soSub$MHCII_RNA > quantile(soSub$MHCII_RNA, 0.75) & soSub$MHCII_ADT_binary,
    "MHCII high","MHCII low")

soTemp <- subset(soSub, cells = WhichCells(soSub, expression = umapSCT_1 < 0))
soTemp <- subset(soTemp, cells = WhichCells(soTemp, expression = umapSCT_2 < 5))
p <- DimPlot(soTemp, group.by = "MHCII_AbSeq_integration", label = FALSE, reduction = "umap_SCT") +
	scale_color_manual(values = c("MHCII high" = "red", "MHCII low" = "gray")) +
	 labs(title = "MO MHCII signature with AbSeq integration")
print(p)
# {png(file = paste0("M_umap_MO_MHCII_sig2.png"), width = 6, height= 3, units = "in", res = 300);plot(p);dev.off()}
{pdf(file = paste0("M_umap_MO_MHCII_sig2.pdf"), width = 8, height= 8);plot(p);dev.off()}
rm(soTemp)

# dotPlot_MO_MHCII_sig2
p <- DotPlot(soSub, assay = "RNA", features = cellList[["MHCII_scores"]], group.by = "MHCII_AbSeq_integration", scale = TRUE) +
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        legend.box = "horizontal",           # Horizontal arrangement
        legend.box.just = "center",          # Center alignment
        legend.spacing.x = unit(0.5, 'cm')) + # Spacing between legends
    guides(color = guide_colorbar(title = "Avg Expression", order = 1),
        size = guide_legend(title = "Pct Exp", order = 2)) + 
	labs(title = "MO MHCII signature with AbSeq integration") +
	theme(plot.title = element_text(size = 10))
p[["data"]][["features.plot"]][is.na(p[["data"]][["features.plot"]])] <- "H2-ia-ie-AbSeq"
print(p)
# {png(file = paste0("M_dotPlot_MO_MHCII_sig2.png"), width = 6, height= 3, units = "in", res = 300);plot(p);dev.off()}
{pdf(file = paste0("M_dotPlot_MO_MHCII_sig2.pdf"), width = 6, height= 3);plot(p);dev.off()}


soSub$day <- factor(soSub$day, levels = c("d8","d38","d120"))
p <- dittoBarPlot(soSub, "MHCII_AUC_status", group.by = "condition", split.by = "day", sub = "cluster 1+3+5+8")
# {png(paste0("V_cellProp_MHCII_clus0+2+4+9+10.png"), width=6, height=4, units = "in", res = 300);plot(p);dev.off()}
{pdf(paste0("V_cellProp_MHCII_clus0+2+4+9+10.pdf"), width=6, height=4);plot(p);dev.off()}

soSub$Sample_Name <- factor(soSub$Sample_Name, levels = c("1_Mock_r1_d8","1_Mock_r2_d8",
													"2_MuHV4_r1_d8","2_MuHV4_r2_d8","2_MuHV4_r3_d8",
													"3_del73_r1_d8","3_del73_r2_d8","3_del73_r3_d8",
													"1_Mock_r1_d38","1_Mock_r2_d38",
													"2_MuHV4_r1_d38","2_MuHV4_r2_d38","2_MuHV4_r3_d38",
													"3_del73_r1_d38","3_del73_r2_d38","3_del73_r3_d38",
													"1_Mock_r1_d120","1_Mock_r2_d120","1_Mock_r3_d120",
													"2_MuHV4_r1_d120","2_MuHV4_r2_d120","2_MuHV4_r3_d120",
													"3_del73_r1_d120","3_del73_r2_d120","3_del73_r3_d120"))

soSub$MHCII_AUC_status <- factor(soSub$MHCII_AUC_status, levels = c("MHCII low","MHCII high"))

sce <- as.SingleCellExperiment(soSub)
sce$cluster_col <- soSub@meta.data[["MHCII_AUC_status"]]
	
prop <- propeller(sce, clusters = sce$cluster_col, sample = sce$Sample_Name, group = sce$condition_day)
write.csv2(prop, file = paste0("V_propellerPropReplicate_MHCII_clus1+3+5+8.csv"))

propT <- getTransformedProps(clusters = sce$cluster_col, sample = sce$Sample_Name, transform="logit")
write.csv2(propT[["Proportions"]], file = paste0("V_transformedPropReplicate_MHCII_clus1+3+5+8.csv"))

# BoxPlots pariwise comparisons

dataPropT <- as.data.frame(propT[["TransformedProps"]])
dataProp <- as.data.frame(propT[["Proportions"]])
dataProp$Tfreq <- dataPropT$Freq
group <- c(rep("Mock_d8",4), rep("MuHV4_d8",6), rep("del73_d8",6),
			rep("Mock_d38",4), rep("MuHV4_d38",6), rep("del73_d38",6),
			rep("Mock_d120",6), rep("MuHV4_d120",6), rep("del73_d120",6))		
dataProp$group <- group
dataProp$day <- c(rep("d8",16),rep("d38",16), rep("d120",18))
dataProp$condition <- c(rep("Mock",4),rep("MuHV4",6), rep("del73",6),
						rep("Mock",4),rep("MuHV4",6), rep("del73",6),
						rep("Mock",6),rep("MuHV4",6), rep("del73",6))
dataProp$group <- factor(dataProp$group, levels = c("Mock_d8","MuHV4_d8","del73_d8",
													"Mock_d38","MuHV4_d38","del73_d38",
													"Mock_d120","MuHV4_d120","del73_d120"))

cellClusters <- unique(sce$cluster_col)
for(cell in cellClusters){
    temp <- dataProp %>% filter(clusters == cell)
    temp <- temp %>% filter(group %in% c("Mock_d8", "Mock_d38", "Mock_d120", "MuHV4_d8", "del73_d8",
                                        "MuHV4_d38", "del73_d38", "MuHV4_d120", "del73_d120")) 
    temp %>% group_by(group) %>% identify_outliers(Freq)

    stats <- temp %>% pairwise_t_test(Tfreq ~ group, p.adjust.method = "bonferroni")
    stats <- stats %>% mutate(across(where(is.numeric), ~ round(., 3)))
    stats_sig <- stats %>% dplyr::filter(p.adj.signif != "ns")
    write.csv2(stats, file=paste0("test_V_boxPlot_",cell,"_clus0+2+4+9+10_stats.csv"))

    temp <- temp %>% group_by(group)
    temp$group <- factor(temp$group, levels = c("Mock_d8","MuHV4_d8","del73_d8",
                                                "Mock_d38","MuHV4_d38","del73_d38", 
                                                "Mock_d120","MuHV4_d120", "del73_d120"))

    b <- ggboxplot(temp, x = "group", y = "Freq", color = "group", 
                   ylab = "Proportions", xlab = "", outlier.shape = 19)
    
	# ifelse(stats$p.adj.signif == "ns", hide.ns <- F, hide.ns <- T)
	# ifelse(stats$p.adj.signif =="ns", label <- "p={p}", label <- "p.adj.signif")
	
    if(nrow(stats_sig) > 0) {
        b <- b + stat_pvalue_manual(stats_sig, 
                                    label = "p.adj.signif",  # Use the significance stars
                                    hide.ns = FALSE,  # Set to FALSE since we pre-filtered
                                    xmin = "group2", 
                                    xmax = "group1",
                                    y.position = max(temp$Freq) + (max(temp$Freq)/15),
                                    step.increase = 0.1)
    }
    
    b <- b + theme(legend.position = "none") + 
             geom_point(position = position_jitter(width = 0.01)) + 
             ggtitle(paste("cluster", cell)) +
             theme(legend.position = "none",
                   axis.text.x = element_text(angle = 45, hjust = 1))
    print(b)
    pdf(paste0("V_boxPlot_",cell,"_clus0+2+4+9+10.pdf"), width = 5, height = 7)
    print(b)
    dev.off()
}

# # Not in the article
# # Last supplemental for Béné

so <- readRDS(file="~/Projects/Arthur/scRNAseq_Analysis3/RawData/seuObjRmCycle.rds")
# so <- readRDS(file="../so.rds")
DefaultAssay(so) <- "RNA"
Idents(so) <- "SCT_clusters_0.6"
so$SCT_clusters_0.6 <- droplevels(so$SCT_clusters_0.6)

so$cell_group <- dplyr::recode(so$SCT_clusters_0.6, "0"="MO", "2"="MO","4"="MO","9"="MO","10"="MO","11"="MO","14"="MO","15"="MO","16"="MO",
								"1"="HSC","2"="HSC","3"="HSC","5"="HSC","6"="HSC","7"="HSC","8"="HSC","12"="HSC","13"="HSC")

Idents(so) <- "cell_group"
MO <- subset(so, ident = "MO")

umap_coords <- Embeddings(MO, reduction = "umap_SCT")
cells_to_keep <- rownames(umap_coords)[umap_coords[, "umapSCT_1"] < 0] #  & umap_coords[, "umapSCT_2"] < 8
MO <- subset(MO, cells = cells_to_keep)

p <- DimPlot(MO, group.by = "SCT_clusters_0.6", reduction = "umap_SCT")
print(p)
# {png("S2F_dotPlot_mo_cluster_marker.png", width = widthEst, height = heightEst, units = "in", res = 300);plot(p);dev.off()}
{pdf("umap_MOfull.pdf", width = 8, height = 10);plot(p);dev.off()}

Idents(MO) <- "SCT_clusters_0.6"
MO <- subset(MO, idents= c("14","15","16"), invert = TRUE)

umap_coords <- Embeddings(MO, reduction = "umap_SCT")
cells_to_keep <- rownames(umap_coords)[umap_coords[, "umapSCT_2"] < 8] # umap_coords[, "umapSCT_1"] < 0 & 
MO <- subset(MO, cells = cells_to_keep)

p <- DimPlot(MO, group.by = "SCT_clusters_0.6", reduction = "umap_SCT")
print(p)
# {png("S2F_dotPlot_mo_cluster_marker.png", width = widthEst, height = heightEst, units = "in", res = 300);plot(p);dev.off()}
{pdf("umap_MO.pdf", width = 8, height = 10);plot(p);dev.off()}

Idents(MO) <- "day"
d8MO <- subset(MO, ident = "d8")

# https://huayc09.github.io/SeuratExtend/articles/Visualization.html

# DotPlot(MO, features = c("Acod1","Cxcl9","Cxcl10","Gbp5","Gbp6","Gbp2","Irf1"), #"Slpi","Cers6","Ifi203",,"Ccl8","Siglecf"), 
		# group.by = "condition_day", cols = c("Lightblue","Red"), scale = TRUE)

# p <- DotPlot(d8MO, features = c("Acod1","Cxcl9","Cxcl10","Gbp5","Gbp6","Slpi"), #"Gbp2",,"Cers6","Ifi203","Irf1","Ccl8","Siglecf"
		# group.by = "condition", cols = c("Lightblue","Red"), scale = TRUE) +
	# # coord_flip() +
	# theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
	# labs(x="Genes", y="Condition", title = "specific markers in MO clusters at d8")
# print(p)
# # {png("S2F_dotPlot_mo_cluster_marker.png", width = widthEst, height = heightEst, units = "in", res = 300);plot(p);dev.off()}
# {pdf("S2F_dotPlot_mo_cluster_marker.pdf", width = widthEst, height = heightEst);plot(p);dev.off()}

# DimPlot2(d8MO, features = c("Acod1","Cxcl9"), cols = "A", reduction = "umap_SCT")

# DimPlot2(d8MO, features = c("Acod1","Cxcl9"), reduction = "umap_SCT", cols = "Spectral-rev")

p <- FeaturePlot3(d8MO, color = "rgb", feature.1 = "Acod1", feature.2 = "Cxcl9",  pt.size = 1,dark.theme = TRUE, reduction = "umap_SCT")
print(p)
# {png("S2F_dotPlot_mo_cluster_marker.png", width = widthEst, height = heightEst, units = "in", res = 300);plot(p);dev.off()}
{pdf("umapMOd8_Acod1_Cxcl10_dark.pdf", width = 8, height = 10);plot(p);dev.off()}

p <- FeaturePlot3.grid(d8MO, features = c("Acod1","Cxcl9"), pt.size = 1, reduction = "umap_SCT")
print(p)
# {png("S2F_dotPlot_mo_cluster_marker.png", width = widthEst, height = heightEst, units = "in", res = 300);plot(p);dev.off()}
{pdf("umapMOd8_Acod1_Cxcl10.pdf", width = 8, height = 10);plot(p);dev.off()}

p <- DimPlot2(d8MO, features = c("Acod1","Cxcl9"), reduction = "umap_SCT", cols = "OrRd")
print(p)
# {png("S2F_dotPlot_mo_cluster_marker.png", width = widthEst, height = heightEst, units = "in", res = 300);plot(p);dev.off()}
{pdf("umapMOd8_Acod1_Cxcl10_split.pdf", width = 16, height = 8);plot(p);dev.off()}

cells <- colnames(d8MO)[d8MO$condition %in% c("Mock", "MuHV4")]
# VlnPlot2(d8MO, features = c("Acod1","Cxcl9"), violin = FALSE, pt.style = "quasirandom", ncol = 1)

# # MDP/GMP
DefaultAssay(so) <- "RNA"
Idents(MO) <- "SCT_clusters_0.6"
MO <- subset(MO, ident = c(0,2,4,9,10))
umap_coords <- Embeddings(MO, reduction = "umap_SCT")
cells_to_keep <- rownames(umap_coords)[umap_coords[, "umapSCT_2"] < 5]
MO <- subset(MO, cells = cells_to_keep)

# AUC/NES scoring and dotplot for GMP & MDP
MO$day <- dplyr::recode(MO$day, "d8"="d8","d38"="d30","d120"="d120")
MO$day <- factor(MO$day, levels = c("d8","d30","d120"))
MO$condition <- factor(MO$condition, levels = c("Mock","MuHV4","del73"))

MO <- AddModuleScore(MO, features = list(c("Cd177", "Cebpb", "Chil3", "Fpr2", "Hopx", "Lcn2", "Sell", "Slpi", "Padi4")), 
						name = "NES_GMP", nbin = 50, ctrl = 100) # bin = 24, ctrl = 50
{MO$NES_GMP <- MO$NES_GMP1;MO$NES_GMP1 <- NULL}

MO <- AddModuleScore(MO, features = list(c("Batf3", "Btla", "C1qb", "Cd209a", "Cd40", "H2.Ab1", "Hpgd", "Id3", "L1cam", "Slamf7", "Tmem176a")), 
						name = "NES_MDP", nbin = 50, ctrl = 100) # bin = 24, ctrl = 50
{MO$NES_MDP <- MO$NES_MDP1;MO$NES_MDP1 <- NULL}

p <- VlnPlot2(MO, features = c("NES_GMP","NES_MDP"), group.by = "day", split.by = "condition", 
			stat.method="wilcox.test", p.adjust.method = "holm", comparisons=list(c("Mock","MuHV4"), c("Mock","del73")), 
			hide.ns = TRUE)
# {png("S2F_dotPlot_mo_cluster_marker.png", width = widthEst, height = heightEst, units = "in", res = 300);plot(p);dev.off()}
{pdf("vln_MDP_GMP_condition_day.pdf", width = 16, height = 10);plot(p);dev.off()}

MO$MDP_class <- ifelse(MO$NES_MDP > 0.25, "MDP_pos", "MDP_neg")

# # Proportion
p <- ClusterDistrPlot(origin = MO$Sample_Name, cluster = MO$MDP_class, condition = MO$condition_day)
print(p)
# {png("S2F_dotPlot_mo_cluster_marker.png", width = widthEst, height = heightEst, units = "in", res = 300);plot(p);dev.off()}
{pdf("prop_MDP_GMP_condition_day.pdf", width = 16, height = 10);plot(p);dev.off()}

sce <- as.SingleCellExperiment(MO)
sce$cluster_col <- MO@meta.data[["MDP_class"]]
	
prop <- propeller(sce, clusters = sce$cluster_col, sample = sce$Sample_Name, group = sce$condition_day)
# write.csv2(prop, file = paste0("H_propellerPropReplicate_MO.csv"))

propT <- getTransformedProps(clusters = sce$cluster_col, sample = sce$Sample_Name, transform="logit")
# write.csv2(propT[["Proportions"]], file = paste0("H_transformedPropReplicate_MO.csv"))

# BoxPlots pariwise comparisons
dataPropT <- as.data.frame(propT[["TransformedProps"]])
dataProp <- as.data.frame(propT[["Proportions"]])
dataProp$Tfreq <- dataPropT$Freq

group_info <- data.frame(
  sample = sce$Sample_Name,
  day = sce$day,
  condition = sce$condition,
  group = sce$condition_day  # This should be the combined variable
) %>%
  distinct()

# Merge with your proportion data
dataProp <- dataProp %>%
  left_join(group_info, by = c("sample" = "sample"))

dataProp$group <- gsub("^[0-9]_","",dataProp$group)
dataProp$group <- factor(dataProp$group, 
                         levels = c("Mock_d8", "MuHV4_d8", "del73_d8",
                                    "Mock_d38", "MuHV4_d38", "del73_d38",
                                    "Mock_d120", "MuHV4_d120", "del73_d120"))

dataProp$day <- factor(dataProp$day, levels = c("d8", "d30", "d120"))
dataProp$condition <- factor(dataProp$condition, levels = c("Mock", "MuHV4", "del73"))

clusters <- unique(sce$MDP_class)
condition_colors <- c("Mock" = "#E69F00", 
                     "MuHV4" = "#56B4E9", 
                     "del73" = "#009E73")

days <- unique(dataProp$day)
clusters <- unique(dataProp$clusters)

for(d in days){
	temp_day <- dataProp %>% filter(day == d) %>% filter(condition %in% c("Mock", "MuHV4", "del73"))
	all_stats <- list()
	
	stats <- temp_day %>% group_by(clusters) %>% pairwise_t_test(
															Tfreq ~ condition, 
															p.adjust.method = "bonferroni", 
															# ref.group = "Mock",
															detailed = TRUE)
	
	stats <- stats %>% add_xy_position(x = "clusters", dodge = 0.8)
	ifelse(stats$p.adj.signif == "ns", hide.ns <- F, hide.ns <- T)
	ifelse(stats$p.adj.signif =="ns", label <- "p={p}", label <- "p.adj.signif")
	
	y_positions <- temp_day %>%
    group_by(clusters) %>%
    summarise(
        max_freq = max(Freq, na.rm = TRUE),
        y.pos = max_freq + (max_freq/15)
    )

	# Then merge with stats
	stats <- stats %>% left_join(y_positions, by = c("clusters" = "clusters"))
	
	p <- ggbarplot(temp_day, x = "clusters", y = "Freq", add = "mean_sd", fill = "condition", color = "condition", 
				palette = condition_colors, position = position_dodge(0.8)) + labs(title = paste0("MO at ",d), 
																				subtitle = "clusters 0+2+4+9+10")
	
	# p <- ggboxplot(temp_day,x = "clusters", y = "Freq", fill = "condition")
	p <- p + stat_pvalue_manual(stats, label=label, hide.ns=hide.ns, xmin = "xmin", xmax = "xmax", #,hide.ns=T
                     y.position = "y.pos", step.increase = 0.1, step.group.by = "clusters")
	
  
  # Create summary data for bar plot (mean ± SEM)
	summary_data <- temp_day %>%
		group_by(clusters, condition) %>%
		summarise(
		  mean_Freq = mean(Freq, na.rm = TRUE),
		  sd_Freq = sd(Freq, na.rm = TRUE),
		  sem_Freq = sd_Freq / sqrt(n()),
		  n = n(),
		  .groups = 'drop'
		) %>%
		mutate(
		  ymax = mean_Freq + sem_Freq,
		  ymin = mean_Freq - sem_Freq
		)
	
	stats$groups <- NULL
	# write.csv2(stats, file = paste0("H_propMO_stats_", d, ".csv"), row.names = FALSE)
	# write.csv2(summary_data, file = paste0("H_propMO_summary_", d, ".csv"), row.names = FALSE)

	print(p)
	# {png(paste0("H_barPlot_propMO_",d,".png"), width = 10, height = 6, units = "in", res = 300);plot(p);dev.off()}
	{pdf(paste0("barPlot_propMO_",d,".pdf"), width = 6, height = 10);plot(p);dev.off()}
}
