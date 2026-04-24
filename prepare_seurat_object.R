library(readr)
library(Seurat)
library(dplyr)
library(Matrix)
library(future)
library(stringr)
library(tidyverse)

# get RNA & Ab count
expName <- c("d8","d38","d120")
so <- list()
for (name in expName){
  counts <- read.csv(paste0("Combined_DBEC_MolsPerCell_",name,".csv"), 
                     sep = ",", header = T, row.names=1, skip = 7)
  Ab <- counts[grep(names(counts),pattern = "pAbO")] 
  # colnames(Ab) <- colnames(Ab) %>% gsub("[.]","_",.)
  RNA <- counts[grep(names(counts),pattern = c("pAbO|Cell"), invert = T)] 
  #%>% clean_names()
  
  # get metadata
  tag <- read.table(paste0("Sample_Tag_Calls_",name,".csv"),
                    sep = ",", header = T, skip = 0)
  metadata <- column_to_rownames(tag, var="Cell_Index")
  barcodeAll <- rownames(counts)
  barcodeToKeep <- tag %>% dplyr::filter(!(Sample_Tag == "Multiplet" | Sample_Tag == "Undetermined")) %>% .$Cell_Index
  barcodeSelection <- barcodeAll[barcodeAll %in% barcodeToKeep]
  
  AbToKeep <- select(Ab, c( #"Cell_Index",
    #"Siglec.F.Siglecf.AMM2013.pAbO",
    "I.A_I.E.H2.Ab_Ad_Aq_Ed_Ek.AMM2019.pAbO",
    "CD274.Cd274.AMM2038.pAbO",
    #"CD11c.HL3.Itgax.AMM2008.pAbO",
    #"Ly.6G.Ly6g.AMM2009.pAbO",
    #"Ly.6A_Ly.6E.Ly6a_Ly6e.AMM2026.pAbO",
  ))
  colnames(AbToKeep) <- c("H2_ia_ie_AbSeq","Cd274_AbSeq") # "Siglecf_AbSeq","Cd11c","Ly_6g_AbSeq","Ly_6a_AbSeq"
  barcodeAbFiltered <- barcodeSelection[barcodeSelection %in% rownames(AbToKeep)]
  
  # create a seurat object 
  AbToKeep <- as.data.frame(as.matrix(t(AbToKeep)))
  RNA <- as.data.frame(as.matrix(t(RNA)))
  seuObj <- CreateSeuratObject(counts = RNA, meta.data = metadata, project = name)
  seuObj[["ADT"]] <- CreateAssayObject(AbToKeep[, colnames(seuObj)])
  seuObj <- seuObj[, barcodeAbFiltered] # barcodeSelection to switch with barcodeAbFiltered if we filter cells on Ab
  
  rm(Ab, AbToKeep, counts, metadata, RNA, tag, barcodeAbFiltered, barcodeAll, barcodeSelection, barcodeToKeep)
  gc()
  
  # create additional metadata
  seuObj$Sample_Name <- str_replace(seuObj$Sample_Name,pattern="[-]", replacement = "_")
  seuObj$time <- name #"d8","d38","d120"
  seuObj$condition_replicate <- seuObj$Sample_Name
  seuObj$condition <- gsub("_.*","",seuObj$Sample_Name)
  seuObj$condition_time <- paste0(seuObj$condition,"_",name)
  seuObj$Sample_Name <- paste0(seuObj$Sample_Name,"_",name)
  so[[name]] <- seuObj
  rm(seuObj)
  gc()
}

