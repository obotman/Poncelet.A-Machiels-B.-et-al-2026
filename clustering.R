rm(list=ls())
library(Seurat)
library(SeuratData)
library(tidyseurat)
library(clustree)
library(dplyr)
library(tidyverse)

so <- readRDS(file="~/seurat_object.rds")

# **********************************
# Normalization, UMAP and clustering 
# **********************************
so <- ScaleData(so, assay = "RNA", vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")) 
so <- SCTransform(so,method = "glmGamPoi",return.only.var.genes = FALSE, vst.flavor = "v2", vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
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
so <- FindClusters(so, resolution = r1, cluster.name = paste0("SCT_clusters_",r1))


