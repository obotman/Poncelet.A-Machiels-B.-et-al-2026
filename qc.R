rm(list = ls())
library(ggplot2)
library(Seurat)
library(readr)
library(dplyr)
library(tidyverse)
library(tidyseurat)
library(dittoSeq)
library(gdata)

plan("multicore", workers = 10)
options(future.globals.maxSize= 20 * 1024 * 1024 * 1024)

setwd("~/Projects/Arthur/scRNAseq_Analysis2/QC")
so <- readRDS("../RawData/so.rds") #soList

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
  png(paste0("QcNfeature_",name,".png"), width = 6, height = 5, units = "in", res = 300)
  print(VlnPlot(so[[name]], features = "nFeature_RNA", group.by = "Sample_Name", pt.size=0.00000000001, raster = FALSE) + 
          labs(x="SampleID", title="Nb of genes detected per cell")+
          geom_hline(yintercept = 6500, color = "red", linetype= "dashed", size=1) + # &
		  geom_hline(yintercept = 1000, color = "red", linetype= "dashed", size=1))
  dev.off()
  
  png(paste0("QcNcount_",name,".png"), width = 6, height = 5, units = "in", res = 300)
  print(VlnPlot(so[[name]], features = "nCount_RNA", group.by = "Sample_Name", pt.size=0.00000000001, raster = FALSE) + 
          labs(x="SampleID", title="No. of UMIs per cell"))
  dev.off()
  
  png(paste0("QcPercentMt_",name,".png"), width = 6, height = 5, units = "in", res = 300)
  print(VlnPlot(so[[name]], features = "percent.mt", group.by = "Sample_Name", pt.size=0.00000000001, raster = FALSE) + NoLegend() + 
          labs(x="SampleID", title="% of mitochondrial UMIs per cell") +
          geom_hline(yintercept = 10, color = "red", linetype= "dashed", size=1) &
          geom_text(aes(2, 14, label = "10%"), color="red", size=5))
  dev.off()
  
  png(paste0("QcPercentMtVsNcount_",name,".png"), width = 6, height = 5, units = "in", res = 300)
  print(FeatureScatter(so[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt", raster = FALSE)  + 
          labs(x="No. of UMIs per cell", y="% of mitochondrial UMIs per cell", title=NULL) + NoLegend() +
          geom_hline(yintercept = 10, color = "red", linetype= "dashed", size=1) &
          geom_text(aes(150000, 14, label = "10%"), color="red", size=5))
  dev.off()
  
  png(paste0("QcPercentMtVsNfeature_",name,".png"), width = 6, height = 5, units = "in", res = 300)
  print(FeatureScatter(so[[name]], feature1 = "nFeature_RNA", feature2 = "percent.mt", raster = FALSE) + 
          labs(x="No. of genes detected per cell", y="% of mitochondrial UMIs per cell", title=NULL) + NoLegend())
  dev.off()
  
  png(paste0("QcNfeatureVsNcount_",name,".png"), width = 6, height = 5, units = "in", res = 300)
  print(FeatureScatter(so[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE) + 
          labs(x="No. of UMIs per cell", y="No. of genes detected per cell", title=NULL) + NoLegend())
  dev.off()
  
  png(paste0("QcNfeatureHist_",name,".png"), width = 6, height = 5, units = "in", res = 300)
  print(hist(so[[name]]@meta.data$`nFeature_RNA`, 50, xlab = "No. of features", ylab = "Frequency", main = "Histogram of features/cell"))
  dev.off()
  
  png(paste0("QcNmtHist_",name,".png"), width = 6, height = 5, units = "in", res = 300)
  print(hist(so[[name]]@meta.data$`percent.mt`,50, xlab = "Percent of mitochondrial genes", ylab = "Frequency", main = "Histogram of mitochondrial genes/cell"))
  dev.off()
  
  # https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

  png(paste0("QcComplexicity_",name,".png"), width = 6, height = 5, units = "in", res = 300)
  d <- so[[name]]@meta.data %>% ggplot(aes(x=log10GenesPerUMI, color = Sample_Name, fill=Sample_Name)) +
  geom_density(alpha = 0.2) + theme_classic() + geom_vline(xintercept = 0.8)
  print(d)
  dev.off()
	
  png(paste0("QcBoxPlotNrna_",name,".png"), width = 6, height = 5, units = "in", res = 300)
  b <- so[[name]]@meta.data %>% ggplot(aes(x=Sample_Name, y=log10(nFeature_RNA), fill=Sample_Name)) + 
  geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) + ggtitle("NCells vs NGenes") %>% print()
  print(b)
  dev.off()
  rm(b,d)
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

lapply(so,ncol)
# $d8
# [1] 18477
# $d38
# [1] 22668
# $d120
# [1] 22397
# total = 63542 cells after filtering



