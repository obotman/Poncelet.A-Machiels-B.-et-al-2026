rm(list=ls())
library(Seurat)
library(SeuratData)
library(tidyseurat)
library(SeuratWrappers)
library(clustree)
library(glmGamPoi)
library(future)
library(dplyr)
library(tidyverse)
library(readr)
library(scCustomize)
library(gdata)
library(ComplexHeatmap)
library(ggplot2)
library(dittoSeq)
library(cowplot)
library(scater)
library(RColorBrewer)
library(SingleCellExperiment)
library(SingleR)
library(patchwork)
library(gdata)
library(SCpubr)

plan("multicore", workers = 1)
options(future.globals.maxSize= 39 * 1024 * 1024 * 1024)# Normalization SCT

setwd("~/Projects/Arthur/scRNAseq_Analysis3/Norm_PCA_UMAP/")
# so <- readRDS(file="~/Projects/Arthur/scRNAseq_Analysis2/RawData/seuObjFull.rds")
so <- readRDS(file="~/Projects/Arthur/scRNAseq_Analysis2/RawData/seuObjRmCycle.rds")
# soList <- readRDS(file="~/Projects/Arthur/scRNAseq_Analysis2/RawData/soList.rds")

# **************************************************************************************
# Method 3: SCT & RNA normalization after merging with vars.to.regress = mt + cell cycle
# **************************************************************************************
# https://satijalab.org/seurat/articles/cell_cycle_vignette.html

# var.features <- SelectIntegrationFeatures(object.list = soList) # , nfeatures = 3000
# rm(soList)
# gc()

# so[["RNA"]] <- split(so[["RNA"]], f = so$experiment) <- to redo with ?

DefaultAssay(so) <- "RNA"
so <- ScaleData(so, assay = "RNA", vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")) 
DefaultAssay(so) <- "RNA"
so <- SCTransform(so,method = "glmGamPoi",return.only.var.genes = FALSE, vst.flavor = "v2", vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
gc()
# DefaultAssay(so) <- "SCT"
# VariableFeatures(so) <- var.features

DefaultAssay(so) <- "RNA"
so <- RunPCA(so, assay = "RNA", reduction.name = "pca_RNA")
so <- RunUMAP(so, assay = "RNA", reduction = "pca_RNA", reduction.name = "umap_RNA", dims = 1:30)
so <- FindNeighbors(so, assay = "RNA", reduction = "pca_RNA", dims = 1:30)
clusTrials <- data.frame(Barcode = rownames(so@meta.data))
for(r1 in c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1)){
     so <- FindClusters(so, reduction = "pca_RNA", resolution = r1, cluster.name = paste0("RNA_clusters_",r1)) #
     png(paste0("umapClusteringRNA_res",r1,".png"), width = 10, height = 8, units = "in", res = 300)
     print(DimPlot(so, reduction = "umap_RNA", group.by = paste0("RNA_clusters_",r1), label = TRUE, label.size = 6))
     dev.off()
     # soClus <- FindClusters(soNorm, resolution = r1, future.seed=TRUE) !! seed to fix
     clusColname <- paste0("Res", r1)
     clusTrials[, clusColname] <- so@meta.data$seurat_clusters
}

png("ClusteringTreeRNA10.png", width = 10, height = 12, units = "in", res = 300)
print(clustree(clusTrials, prefix = "Res"))
dev.off()

DefaultAssay(so) <- "SCT"
so <- RunPCA(so, assay = "SCT", reduction.name = "pca_SCT")
so <- RunUMAP(so, assay = "SCT", reduction = "pca_SCT", reduction.name = "umap_SCT", dims = 1:30)
so <- FindNeighbors(so, assay = "SCT", reduction = "pca_SCT", dims = 1:30)
clusTrials <- data.frame(Barcode = rownames(so@meta.data))
for(r1 in c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1)){
     so <- FindClusters(so, reduction = "pca_SCT", resolution = r1, cluster.name = paste0("SCT_clusters_",r1)) #
     png(paste0("umapClusteringSCT_res",r1,".png"), width = 10, height = 8, units = "in", res = 300)
     print(DimPlot(so, reduction = "umap_SCT", group.by = paste0("SCT_clusters_",r1), label = TRUE, label.size = 6))
     dev.off()
     # soClus <- FindClusters(soNorm, resolution = r1, future.seed=TRUE) !! seed to fix
     clusColname <- paste0("Res", r1)
     clusTrials[, clusColname] <- so@meta.data$seurat_clusters
}

png("ClusteringTreeSCT10.png", width = 10, height = 12, units = "in", res = 300)
print(clustree(clusTrials, prefix = "Res"))
dev.off()

DefaultAssay(so) <- "SCT"
r1 = 0.7
so <- FindClusters(so, resolution = r1, cluster.name = paste0("SCT_clusters_",r1)) # reduction = "pca_SCT",

png(paste0("umapClusteringSCT_res",r1,".png"), width = 10, height = 8, units = "in", res = 300)
print(DimPlot(so, reduction = "umap_SCT", group.by = paste0("SCT_clusters_",r1), label = TRUE, label.size = 6))
dev.off()

png("ClusteringTreeSCT.png", width = 10, height = 12, units = "in", res = 300)
print(clustree(clusTrials, prefix = "Res"))
dev.off()

png("pca_RNA_CellCycle.png", width = 10, height = 8, units = "in", res = 300)
DimPlot(so, reduction = "pca_RNA", group.by = "Phase")
dev.off()

png("pca_SCT_CellCycle.png", width = 10, height = 8, units = "in", res = 300)
DimPlot(so, reduction = "pca_SCT", group.by = "Phase")
dev.off()

png("cellCycleGroupSplitSCT.png", width = 24, height = 8, units = "in", res = 300)
dittoDimPlot(so, "Phase", size = 1, reduction.use="umap_SCT", do.ellipse = TRUE, shape.by = "condition", split.by= "day")
dev.off()

png("umap_SCT_CellCycle.png", width = 24, height = 8, units = "in", res = 300)
DimPlot(so, reduction = "umap_SCT", group.by = c("APO_annotation","Phase","SCT_clusters_0.6"), label = TRUE)
dev.off()

so$condition_day <- factor(so$condition_day, levels = c("Mock_d8","Mock_d38","Mock_d120","del73_d8","del73_d38","del73_d120","MuHV4_d8","MuHV4_d38","MuHV4_d120"))  

png("umap_SCT_CellCycle_split.png", width = 24, height = 14, units = "in", res = 300)
DimPlot(so, reduction = "umap_SCT", group.by = "SCT_clusters_0.6", split.by = "condition_day", ncol=3, pt.size = 0.5, label = TRUE)
dev.off()

# HighLight cells per condition and split by d8 and d60

Idents(so) <- "condition"
highlight <- so$condition %>% unique
for (name in highlight) {
	cell <- WhichCells(so, idents = name)
	s <- paste0("umapHighlight_",name,".png")
	png(s, width = 24, height = 8, units = "in", res = 200)
	print(DimPlot(so, reduction = "umap_SCT", cells.highlight=cell, split.by="day", pt.size = 0.2, sizes.highlight = 0.1, cols.highlight = c("darkblue", "darkred"), cols= "grey") + NoLegend() + labs(title=s))
	dev.off()
}

# *********************
# Pseudotime simulation
# *********************

so$dayInt <- dplyr::recode(so$day, d8 = 20, d38 = 15, d60 = 5)
so$dayInt<-as.numeric(so$dayInt)

png(paste0("featurePlotD8-D60.png"), width = 10, height = 8, units = "in", res = 300)
print(FeaturePlot(so, "dayInt"))
dev.off()

so$mockInt <- dplyr::recode(so$mock_day, Mock = 5, d8 = 20, d38 = 15, d60 = 10)
so$mockInt<-as.numeric(so$mockInt)

png(paste0("featurePlotD0-D60.png"), width = 10, height = 8, units = "in", res = 300)
print(FeaturePlot(so, "mockInt"))
dev.off()

# ***************************
# Cell proportion in clusters
# ***************************

DefaultAssay(so) <- "SCT"
so$SCT_clusters_0.6 <- factor(so$SCT_clusters_0.6, levels = c(0:16))

png("cellProp_condition.png", width = 10, height = 8, units = "in", res = 300)
print(dittoBarPlot(so, "condition", group.by="SCT_clusters_0.6",
		x.reorder = c(1,2,10,11,12,13,14,15,16,17,3,4,5,6,7,8,9)))
dev.off()
png("cellProp_conditionSplitNoOrder.png", width = 12, height = 8, units = "in", res = 300)
print(dittoBarPlot(so, "condition", split.by = "day", group.by = "SCT_clusters_0.6")),
		#x.reorder = c(1,2,10,11,12,13,14,15,16,17,3,4,5,6,7,8,9)))
dev.off()
png("cellProp_condition_day.png", width = 10, height = 8, units = "in", res = 300)
print(dittoBarPlot(so, "condition_day", group.by = "SCT_clusters_0.6",
		x.reorder = c(1,2,10,11,12,13,14,15,16,17,3,4,5,6,7,8,9)))
dev.off()
png("cellProp_day.png", width = 10, height = 8, units = "in", res = 300)
print(dittoBarPlot(so, "day", group.by = "SCT_clusters_0.6",
		x.reorder = c(1,2,10,11,12,13,14,15,16,17,3,4,5,6,7,8,9))) 
dev.off()

md <- so@meta.data %>% as.data.table
md[, .N, by = c("condition_day", "SCT_clusters_0.6")]
write.csv2(as.data.frame(md[, .N, by = c("condition_day", "SCT_clusters_0.6")] %>% dcast(., condition_day ~ SCT_clusters_0.6, value.var = "N")), file="clusterTable.csv")

png("cellProp_ncellClusters.png")
ggplot(so@meta.data, aes(SCT_clusters_0.6, fill=condition_day)) + geom_bar(stat = "count")
dev.off()

# # FindAllMarkers 

clusGeneMarkers <- FindAllMarkers(so, group.by = "SCT_clusters_0.6")
clusGeneMarkersTable <- as_tibble(clusGeneMarkers) %>% dplyr::select(gene, cluster, everything())
write_csv(clusGeneMarkersTable, "FindAllGeneMarkersOfClusters_SCT_0.6.csv")

# Gene curation
# clusGeneMarkersTableC <- clusGeneMarkersTable  %>% filter(!grepl('^Rp', gene))
# clusGeneMarkersTableC <- clusGeneMarkersTableC %>% filter(!grepl('^mt*', gene))
# clusGeneMarkersTableC <- clusGeneMarkersTableC %>% filter(!grepl('^Tub', gene))
# clusGeneMarkersTableC <- clusGeneMarkersTableC %>% filter(!grepl('^Gm', gene))
# clusGeneMarkersTableC <- clusGeneMarkersTableC %>% filter(!grepl('*Rik', gene))

# Select top upregulated markers of each cluster
topGeneMarkers <- clusGeneMarkersTable %>%
  dplyr::filter(avg_log2FC > 0, pct.1 - pct.2 > 0.1, p_val_adj < 0.05) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(5, avg_log2FC)

# Make a dot plot
geneVector <- unique(topGeneMarkers$gene)

png("DotPlotClusMarkers.png", width = 18, height = 6, units = "in", res = 300)
print(DotPlot(so, features = geneVector, group.by = "SCT_clusters_0.6") + 
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
		labs(x="Genes", y="Clusters", title="Gene expression per clusters"))
dev.off()

png("DotPlotClusteredClusMarkers.png", width = 8, height = 10, units = "in", res = 300)
print(Clustered_DotPlot(so, features = geneVector, k=17, group.by = "SCT_clusters_0.6"))
dev.off()

# Violin plot Arthur annotation

geneVector <- c("Gngt2","Ear2", # Ly6c low => 0
				"Vcan","Ccr2","Ly6c1","Ly6c2", # Ly6c high => 2 VCAN- + 10 VCAN+
				"Itgax","Zbtb46","Ciita","H2-ia-ie-AbSeq", # DC => 11, 14 , 15 CIITA+ & 16
				"Cd274","Ly6a","CD274-AbSeq", # Activated Mo => 9
				"Pf4","Cavin2","Mef2c","Pbx1","Rab27b", # MPP2 plateled oriented => 7
				"Rhd","Kel","Ermap",#"Mpp2",# MPP2 erythrocyte oriented => 12
				"Esam","Fgd5","Mmrn1","Pdzk1ip1","Plscr4","Procr","Tinagl1", # LT-HSC => 8
				"Cd34","Hlf","Hoxa9","Flt3", # HSC => ?
				"Elane","Mpo","Ms4a3",#"Mpp3", # MPP3 => 5
				"Dntt",#"Mpp4", # MMP4 cycling cells => 6
				"Cdk1" # cycling cells => 3 ?
				)

png("vlnAPOannotation.png", width = 20, height = 30, units = "in", res = 300)
print(VlnPlot(so, assay = "RNA", group.by = "SCT_clusters_0.6", features = geneVector, stack = TRUE, flip = TRUE))
dev.off()

png("DotPlotClusMarkersAPOonClus.png", width = 8, height = 10, units = "in", res = 300)
print(Clustered_DotPlot(so, features = geneVector, k=17, group.by = "SCT_clusters_0.6"))
dev.off()

png("DotPlotClusMarkersAPOannotation.png", width = 8, height = 10, units = "in", res = 300)
print(Clustered_DotPlot(so, features = geneVector, k=13, group.by = "APO_annotation"))
dev.off()

# ADT
DefaultAssay(so) <- "ADT"

png("adtCD274.png", width = 10, height = 8, units = "in", res = 300)
FeaturePlot(so, features="CD274-AbSeq")
dev.off()

png("adtH2iaie.png", width = 10, height = 8, units = "in", res = 300)
FeaturePlot(so, features="H2-ia-ie-AbSeq")
dev.off()

png("adtVlnH2iaie.png", width = 10, height = 10, units = "in", res = 300)
print(VlnPlot(so, assay = "ADT", group.by = "SCT_clusters_0.6", features = "H2-ia-ie-AbSeq", pt.size = 0.1))
dev.off()

png("adtVlnCD274.png", width = 10, height = 10, units = "in", res = 300)
print(VlnPlot(so, assay = "ADT", group.by = "SCT_clusters_0.6", features = "CD274-AbSeq", pt.size = 0.1))
dev.off()

png("adtVlnH2iaieSplit1.png", width = 30, height = 10, units = "in", res = 300)
print(VlnPlot(so, assay = "ADT", group.by = "condition_day", split.by = "SCT_clusters_0.6", features = "H2-ia-ie-AbSeq", pt.size = 0.1))
dev.off()

png("adtVlnH2iaieSplit2.png", width = 30, height = 10, units = "in", res = 300)
print(VlnPlot(so, assay = "ADT", group.by = "SCT_clusters_0.6", split.by = "condition_day", features = "H2-ia-ie-AbSeq", 
			 layer = "data", pt.size = 0.1))
dev.off()

png("adtVlncD274Split1.png", width = 30, height = 10, units = "in", res = 300)
print(VlnPlot(so, assay = "ADT", group.by = "condition_day", split.by = "SCT_clusters_0.6", features = "CD274-AbSeq", 
			 layer = "data", pt.size = 0.1))
dev.off()

png("adtVlncD274Split2.png", width = 30, height = 10, units = "in", res = 300)
print(VlnPlot(so, assay = "ADT", group.by = "SCT_clusters_0.6", split.by = "condition_day", features = "CD274-AbSeq", 
			 layer = "data", pt.size = 0.1))
dev.off()

for (c in 0:16){
png(paste0("adtVlnH2iaie_Clus",c,".png"), width = 10, height = 10, units = "in", res = 300)
print(VlnPlot(so, assay = "ADT", group.by = "SCT_clusters_0.6", split.by = "condition_day", features = "H2-ia-ie-AbSeq", 
				idents = c, layer = "data", pt.size = 0.1))
dev.off()
}

for (c in 0:16){
png(paste0("adtVlnCD274_Clus",c,".png"), width = 10, height = 10, units = "in", res = 300)
print(VlnPlot(so, assay = "ADT", group.by = "SCT_clusters_0.6", split.by = "condition_day", features = "CD274-AbSeq", 
				idents = c, layer = "data", pt.size = 0.1))
dev.off()
}

# FindMarkers in each clusters between condition
DefaultAssay(so) <- "ADT"

geneVector <- c("CD274-AbSeq", "H2-ia-ie-AbSeq")

sampleIDVector <- c("Mock_d8","Mock_d38","Mock_d120",
					"MuHV4_d8","MuHV4_d38","MuHV4_d120",
					"del73_d8","del73_d38","del73_d120" )

markers <- list() 
for (clus in 0:16){
Idents(so) <- "SCT_clusters_0.6"
subClus <- subset(so, ident = clus)
Idents(subClus) <- "condition_day"
markers[[as.character(clus)]] <- FindAllMarkers(subClus, only.pos = FALSE)

png(paste0("DotPlotADTmarkers_cluster_",clus,".png"), width = 4, height = 6, units = "in", res = 300)
print(DotPlot(subClus, features = geneVector, group.by = "condition_day") + 
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
		labs(x="Genes", y="treatment", title=paste0("AbSeq expression in cluster ",clus))+
		theme(plot.title = element_text(size=10)))
dev.off()
					
exprsMatrix <- as(subClus[["ADT"]]$data[geneVector,], "matrix")
exprsTable <- tibble(sampleID = character(0), Gene = character(0), Total = numeric(0), Positive = numeric(0), AvgExpr = numeric(0), MaxExpr = numeric(0))
for(gene1 in geneVector){
	for(sampleID1 in sampleIDVector){
		sampleBc <- colnames(subClus)[subClus$condition_day == sampleID1]
		nTotal <- length(sampleBc)
		nPositive <- sum(exprsMatrix[gene1, sampleBc] > 0)
		exprSum <- sum(exprsMatrix[gene1, sampleBc], na.rm = TRUE)
		if(nPositive == 0){
          avgExpr1 <- 0
        } else{
          avgExpr1 <- exprSum / nPositive
        }
		maxExpr1 <- max(exprsMatrix[gene1, sampleBc], na.rm = TRUE)
		exprsTable <- exprsTable %>% add_row(sampleID = sampleID1, Gene = gene1, Total = nTotal, Positive = nPositive, AvgExpr = avgExpr1, MaxExpr = maxExpr1)
		}
	}

dfSel <- exprsTable %>% dplyr::select(sampleID, Gene, AvgExpr)
dfWide <- tidyr::pivot_wider(dfSel, names_from = "sampleID", values_from = "AvgExpr")
exprsData <- data.matrix(dfWide[, 2:ncol(dfWide)])
rownames(exprsData) <- dfWide$Gene
exprsData[is.na(exprsData)] <- 0
hmDataScale <- t(scale(t(exprsData)))
# hmFilename <- gsub("diffExp", "heatmapTop", comparisubClusnList[[c1]]$deFilename)
# hmFilename <- gsub("csv", "png", hmFilename)
widthEst <- 2 + length(sampleIDVector) * 0.5
heightEst <- 3 + length(geneVector) * 0.1
colTitle <- paste(length(geneVector), "top DEGs")
png(paste0("heatmapADTcluster_",clus,".png"), width = widthEst, height = heightEst, units = "in", res = 300)
# print(Heatmap(hmDataScale, name = "Z-score", column_names_rot = 45, column_title = colTitle, row_names_gp = gpar(fontsize = 8)))
print(Heatmap(hmDataScale,cluster_columns=FALSE, name = "Z-score", show_column_dend = FALSE, column_names_rot = 45, column_title = paste0("ADT in cluster ",clus), row_names_gp = gpar(fontsize = 8)))
dev.off()
# row_km = 2, , column_km = 4
rm(subClus)
gc()
}

for (i in 1:17){
markers[[i]] <- as_tibble(markers[[i]]) %>% dplyr::select(gene, cluster, everything())
}
save(markers, file = "adtPerClusters.Rdata")

# **************************
# BoxPlot in cluster 1,5 & 8
# **************************
# https://enblacar.github.io/SCpubr-book/functions/BoxPlots.html
library(SCpubr)
DefaultAssay(so) <- "RNA"

sampleIDVector <- c("Mock_d8","Mock_d38","Mock_d120",
					"MuHV4_d8","MuHV4_d38","MuHV4_d120",
					"del73_d8","del73_d38","del73_d120" )
					
so$condition_day <- factor(so$condition_day, levels = c("Mock_d8","MuHV4_d8","del73_d8",
														"Mock_d38", "MuHV4_d38","del73_d38",
														"Mock_d120","MuHV4_d120","del73_d120"))
geneVector <- c("CD274-AbSeq", "H2-ia-ie-AbSeq")

for (clus in 0:6){
	Idents(so) <- "SCT_clusters_0.6"
	subClus <- subset(so, ident = clus)
		for (marker in geneVector){
		png(paste0("bloxPlotADT_",marker,"_cluster_",clus,".png"), width = 8, height = 8, units = "in", res = 300)
		print(do_BoxPlot(subClus, feature = marker, group.by="condition_day", test = "wilcox.test", use_test = TRUE, assay = "ADT",
												comparisons = list(c("Mock_d8","MuHV4_d8"),
																	c("Mock_d8","del73_d8"),
																	c("Mock_d38","MuHV4_d38"),
																	c("Mock_d38","del73_d38"),
																	c("Mock_d120","MuHV4_d120"),
																	c("Mock_d120","del73_d120")))
																	)
		dev.off()
	}
}

# Compare ADT to mRNA expression


features <- c("CD274-AbSeq", "H2-ia-ie-AbSeq", "Cd274", "H2.Aa") #,"H2.Ab1", "H2.Eb1", "H2.Eb2")
p<-VlnPlot(so, features = features, group.by = "SCT_clusters_0.6", assay="RNA", ncol=2) #, split.by="condition_day"
# SCpubr::do_ViolinPlot(so, features = features, group.by = "SCT_clusters_0.6")
for (alpha in 1:4){p[[alpha]]$layers[[2]]$aes_params$alpha <- 0.05}
png(paste0("VlnPlotADTvsmRNAclus4.png"), width = 12, height = 12, units = "in", res = 300)
print(p)
dev.off()


features <- c("CD274-AbSeq", "H2-ia-ie-AbSeq", "Cd274", "H2.Aa")
plots <- list()
plots[[1]] <- do_BoxPlot(so, assay = "ADT", feature = "CD274-AbSeq", group.by = "SCT_clusters_0.6", legend.position = "none")
plots[[2]] <- do_BoxPlot(so, assay = "ADT", feature = "H2-ia-ie-AbSeq", group.by = "SCT_clusters_0.6", legend.position = "none")  
plots[[3]] <- do_BoxPlot(so, assay = "RNA", feature = "Cd274", group.by = "SCT_clusters_0.6", legend.position = "none")  
plots[[4]] <- do_BoxPlot(so, assay = "RNA", feature = "H2.Aa", group.by = "SCT_clusters_0.6", legend.position = "none")   
png(paste0("BoxPlotADTvsmRNA.png"), width = 12, height = 12, units = "in", res = 300)
print(wrap_plots(plots, ncol = 2, byrow = TRUE)) #+
		# plot_annotation(paste0("cluster ",clus, " in experiment integration"), theme=theme(plot.title=element_text(hjust=0.5)))
dev.off()

png(paste0("umapADTvsmRNA_CD274.png"), width = 20, height = 6, units = "in", res = 300)
FeaturePlot(so, features = c("CD274-AbSeq","Cd274"), blend = TRUE, blend.threshold=0.1, reduction = "umap_SCT")
dev.off()
png(paste0("umapADTvsmRNA_H2aa.png"), width = 20, height = 6, units = "in", res = 300)
FeaturePlot(so, features = c("H2-ia-ie-AbSeq","H2.Aa"), blend = TRUE, blend.threshold=0.1, reduction = "umap_SCT")
dev.off()

# Focus on cluster 4
subSet <- subset(so, idents="4", group.by = "SCT_clusters_0.6")

ncol(subset(so, subset = `H2-ia-ie-AbSeq` >0, idents = "4"))
# 4749
ncol(subset(so, subset = `CD274-AbSeq` >0, idents = "4"))
# 4257
ncol(subset(so, idents="4"))
# 4767

FindMarkers(so, ident.1 = "del73_d8", ident.2 = "Mock_d8",group.by = "condition_day", subset.ident = "4", features=c("Cd274","H2.Aa"))
FindMarkers(so, ident.1 = "MuHV4_d8", ident.2 = "Mock_d8",group.by = "condition_day", subset.ident = "4", features=c("Cd274","H2.Aa"))


png(paste0("VlnPlotADTvsmRNAH2clus4.png"), width = 10, height = 8, units = "in", res = 300)
print(do_ViolinPlot(subClus, features = "H2-ia-ie-AbSeq", group.by = "condition_day", assay="ADT", ncol = 2))
dev.off()

png(paste0("VlnPlotADTvsmRNAH2Aaclus4.png"), width = 10, height = 8, units = "in", res = 300)
print(do_ViolinPlot(subClus, features = "H2.Aa", group.by = "condition_day", assay="RNA", ncol = 2))
dev.off()

png(paste0("VlnPlotADTvsmRNA_CD274adtclus4.png"), width = 10, height = 8, units = "in", res = 300)
print(do_ViolinPlot(subClus, features = "CD274-AbSeq", group.by = "condition_day", assay="ADT"))
dev.off()

png(paste0("VlnPlotADTvsmRNA_Cd274rnaclus4.png"), width = 10, height = 8, units = "in", res = 300)
print(do_ViolinPlot(subClus, features = "Cd274", group.by = "condition_day", assay="RNA"))
dev.off()


# ***************************************************
# Heatmap on condition_day inside cluster 0,1,2,4,5,8
# ***************************************************

DefaultAssay(so) <- "RNA"
lfc = 1
for (clus in c(0,1,2,4,5,8)){
	Idents(so) <- "SCT_clusters_0.6"
	subClus <- subset(so, ident = clus)
	Idents(subClus) <- "condition_day"
	sample.markers <- FindAllMarkers(subClus, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, slot = "data")
	sampleGeneMarkersTable <- as_tibble(sample.markers) %>% dplyr::select(gene, cluster, everything())
	colnames(sample.markers)[6] <- "condition_day"
	colnames(sampleGeneMarkersTable)[2] <- "condition_day"
	write_csv2(sampleGeneMarkersTable, paste0("FindAllGeneMarkersOfcluster_",clus,".csv"))
	clusGeneMarkersTable <- read.csv(paste0("FindAllGeneMarkersOfcluster_",clus,".csv"), sep= ";", dec = ",")

	topMarkers <- clusGeneMarkersTable %>%
	  dplyr::filter(abs(avg_log2FC) > lfc, pct.1 - pct.2 > 0.1, p_val_adj < 0.05) %>%
	  dplyr::group_by("condition_day") # %>%
	  # dplyr::top_n(50, avg_log2FC)			

	geneVector <- topMarkers$gene %>% unique

	sampleIDVector <- c("Mock_d8","Mock_d38","Mock_d120",
						"MuHV4_d8","MuHV4_d38","MuHV4_d120",
						"del73_d8","del73_d38","del73_d120" )
						
	# exprsMatrix <- as(is.na(subClus[["RNA"]]@data[geneVector,], "matrix"))
	# convert seuObj v5 into v3 version
	# seuObjSubV3 <- Convert_Assay(seurat_object = seuObjSub, convert_to = "V3")
	# subClus[["RNAv3"]] <- as(object = subClus[["RNA"]], Class = "Assay")

	exprsMatrix <- as(subClus[["RNA"]]$data[geneVector,], "matrix")
	exprsTable <- tibble(sampleID = character(0), Gene = character(0), Total = numeric(0), Positive = numeric(0), AvgExpr = numeric(0), MaxExpr = numeric(0))
	for(gene1 in geneVector){
		for(sampleID1 in sampleIDVector){
			sampleBc <- colnames(subClus)[subClus$condition_day == sampleID1]
			nTotal <- length(sampleBc)
			nPositive <- sum(exprsMatrix[gene1, sampleBc] > 0)
			exprSum <- sum(exprsMatrix[gene1, sampleBc], na.rm = TRUE)
			if(nPositive == 0){
			  avgExpr1 <- 0
			} else{
			  avgExpr1 <- exprSum / nPositive
			}
			maxExpr1 <- max(exprsMatrix[gene1, sampleBc], na.rm = TRUE)
			exprsTable <- exprsTable %>% add_row(sampleID = sampleID1, Gene = gene1, Total = nTotal, Positive = nPositive, AvgExpr = avgExpr1, MaxExpr = maxExpr1)
			}
		}

	dfSel <- exprsTable %>% dplyr::select(sampleID, Gene, AvgExpr)
	dfWide <- tidyr::pivot_wider(dfSel, names_from = "sampleID", values_from = "AvgExpr")
	exprsData <- data.matrix(dfWide[, 2:ncol(dfWide)])
	rownames(exprsData) <- dfWide$Gene
	exprsData[is.na(exprsData)] <- 0
	hmDataScale <- t(scale(t(exprsData)))
	# hmFilename <- gsub("diffExp", "heatmapTop", comparisonList[[c1]]$deFilename)
	# hmFilename <- gsub("csv", "png", hmFilename)
	widthEst <- 2 + length(sampleIDVector) * 0.5
	heightEst <- 3 + length(geneVector) * 0.1
	colTitle <- paste0(length(geneVector)," top DEGs, log2FC>",lfc)
	png(paste0("heatmapAbsCluster",clus,"log2FC_",lfc,".2.png"), width = widthEst, height = heightEst, units = "in", res = 300)
	# print(Heatmap(hmDataScale, name = "Z-score", column_names_rot = 45, column_title = colTitle, row_names_gp = gpar(fontsize = 8)))
	print(Heatmap(hmDataScale, name = "Z-score", #row_km = 3, column_km = 3, cluster_columns = TRUE, show_column_dend = FALSE, 
	column_names_rot = 45, column_title = colTitle, row_names_gp = gpar(fontsize = 8), col = mr))
	dev.off()
	# , column_km = 4, col = mr

	# sink("heatmapAbsSCT_clusters_0.4Markers.txt")
	# print(geneVector)
	# sink()
}
# colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
brewer_palette <- "RdBu"
ramp <- colorRampPalette( RColorBrewer::brewer.pal(11, brewer_palette))
mr <- ramp(256)[256:1]

# ************************************************************
# Heatmap in each clusters for each time point with replicates
# ************************************************************

DefaultAssay(so) <- "RNA"
lfc = 0
for (clus in 0:16){
	Idents(so) <- "SCT_clusters_0.6"
	subClus <- subset(so, ident = clus)
	Idents(subClus) <- "condition_day"
	sample.markers <- FindAllMarkers(subClus, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, slot = "data")
	sampleGeneMarkersTable <- as_tibble(sample.markers) %>% dplyr::select(gene, cluster, everything())
	colnames(sample.markers)[6] <- "condition_day"
	colnames(sampleGeneMarkersTable)[2] <- "condition_day"
	write_csv2(sampleGeneMarkersTable, paste0("FindAllGeneMarkersOfcluster_",clus,".csv"))
	clusGeneMarkersTable <- read.csv(paste0("FindAllGeneMarkersOfcluster_",clus,".csv"), sep= ";", dec = ",")

	topMarkers <- clusGeneMarkersTable %>%
	  dplyr::filter(abs(avg_log2FC) > lfc, pct.1 - pct.2 > 0.1, p_val_adj < 0.05) %>%
	  dplyr::group_by("condition_day") # %>%
	  # dplyr::top_n(50, avg_log2FC)			

	geneVector <- topMarkers$gene %>% unique

	sampleIDVector <- c("Mock_d8","Mock_d38","Mock_d120",
						"MuHV4_d8","MuHV4_d38","MuHV4_d120",
						"del73_d8","del73_d38","del73_d120" )
						
	# exprsMatrix <- as(is.na(subClus[["RNA"]]@data[geneVector,], "matrix"))
	# convert seuObj v5 into v3 version
	# seuObjSubV3 <- Convert_Assay(seurat_object = seuObjSub, convert_to = "V3")
	# subClus[["RNAv3"]] <- as(object = subClus[["RNA"]], Class = "Assay")

	exprsMatrix <- as(subClus[["RNA"]]$data[geneVector,], "matrix")
	exprsTable <- tibble(sampleID = character(0), Gene = character(0), Total = numeric(0), Positive = numeric(0), AvgExpr = numeric(0), MaxExpr = numeric(0))
	for(gene1 in geneVector){
		for(sampleID1 in sampleIDVector){
			sampleBc <- colnames(subClus)[subClus$condition_day == sampleID1]
			nTotal <- length(sampleBc)
			nPositive <- sum(exprsMatrix[gene1, sampleBc] > 0)
			exprSum <- sum(exprsMatrix[gene1, sampleBc], na.rm = TRUE)
			if(nPositive == 0){
			  avgExpr1 <- 0
			} else{
			  avgExpr1 <- exprSum / nPositive
			}
			maxExpr1 <- max(exprsMatrix[gene1, sampleBc], na.rm = TRUE)
			exprsTable <- exprsTable %>% add_row(sampleID = sampleID1, Gene = gene1, Total = nTotal, Positive = nPositive, AvgExpr = avgExpr1, MaxExpr = maxExpr1)
			}
		}

	dfSel <- exprsTable %>% dplyr::select(sampleID, Gene, AvgExpr)
	dfWide <- tidyr::pivot_wider(dfSel, names_from = "sampleID", values_from = "AvgExpr")
	exprsData <- data.matrix(dfWide[, 2:ncol(dfWide)])
	rownames(exprsData) <- dfWide$Gene
	exprsData[is.na(exprsData)] <- 0
	hmDataScale <- t(scale(t(exprsData)))
	# hmFilename <- gsub("diffExp", "heatmapTop", comparisonList[[c1]]$deFilename)
	# hmFilename <- gsub("csv", "png", hmFilename)
	widthEst <- 2 + length(sampleIDVector) * 0.5
	heightEst <- 3 + length(geneVector) * 0.1
	colTitle <- paste0(length(geneVector)," top DEGs, log2FC>",lfc)
	png(paste0("heatmapAbsCluster",clus,"log2FC_",lfc,".2.png"), width = widthEst, height = heightEst, units = "in", res = 300)
	# print(Heatmap(hmDataScale, name = "Z-score", column_names_rot = 45, column_title = colTitle, row_names_gp = gpar(fontsize = 8)))
	print(Heatmap(hmDataScale, name = "Z-score", #row_km = 3, column_km = 3, cluster_columns = TRUE, show_column_dend = FALSE, 
	column_names_rot = 45, column_title = colTitle, row_names_gp = gpar(fontsize = 8), col = mr))
	dev.off()
	# , column_km = 4, col = mr

	# sink("heatmapAbsSCT_clusters_0.4Markers.txt")
	# print(geneVector)
	# sink()
}
# colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
brewer_palette <- "RdBu"
ramp <- colorRampPalette( RColorBrewer::brewer.pal(11, brewer_palette))
mr <- ramp(256)[256:1]


# ^^^^^^^^^
# Save data

write_rds(so, file="../../RawData/seuObjRmCycle.rds")
gc()

library(loupeR)

DefaultAssay(so) <- "RNA"
create_loupe_from_seurat(so, output_name = "loupeRmCycle", output_dir = "../../LoupeData", force = T)

# DefaultAssay(so) <- "ADT"
# create_loupe_from_seurat(so, output_name = "loupeAbSeq", output_dir = "../LoupeData", force = T)

keep(so, sure = T)
gc()
# ^^^^^^^^^


gc()

# load(file="~/CellDataBase/cellCycle_g2m_genes_tinyAtlas_mm.Rdata")
# load(file="~/CellDataBase/cellCycle_s_genes_tinyAtlas_mm.Rdata")

# test <- RunPCA(so, features = c(s_genes, g2m_genes))
# DimPlot(test)