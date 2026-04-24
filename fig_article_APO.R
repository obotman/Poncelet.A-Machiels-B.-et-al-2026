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

setwd("~/Projects/Arthur/scRNAseq_Analysis3/Article/Figures")
# so <- readRDS(file="~/Projects/Arthur/scRNAseq_Analysis3/RawData/seuObjClustersFiltered.rds")
so <- readRDS(file="~/Projects/Arthur/scRNAseq_Analysis3/RawData/seuObjRmCycle.rds")
# so <- readRDS(file="../so.rds")
DefaultAssay(so) <- "RNA"
Idents(so) <- "SCT_clusters_0.6"
so$SCT_clusters_0.6 <- droplevels(so$SCT_clusters_0.6)
DimPlot(so, group.by = "SCT_clusters_0.6", reduction = "umap_SCT")
table(so$SCT_clusters_0.6)

so <- subset(so, idents= c("14","15","16"), invert = TRUE)
table(so$SCT_clusters_0.6)

so$cell_group <- dplyr::recode(so$SCT_clusters_0.6, "0"="MO", "2"="MO","4"="MO","9"="MO","10"="MO","11"="MO",
								"1"="HSC","2"="HSC","3"="HSC","5"="HSC","6"="HSC","7"="HSC","8"="HSC","12"="HSC","13"="HSC")

table(so$cell_group)

   MO   HSC 
27990 31901

# ********
# Figure 2
# ********

# - fig A - 
# umapSCT
p <- dittoDimPlot(so, var = "SCT_clusters_0.6", reduction = "umap_SCT", 
					do.label = TRUE, labels.highlight = FALSE, labels.size = 10)
# {png(file = paste0("A_umapSCT_clusters.png"), width = 10, height= 8, units = "in", res = 300);plot(p);dev.off()}
{pdf(file = paste0("A_umapSCT_clusters.pdf"), width = 12, height= 8);plot(p);dev.off()}

so$cluster_annotated <- dplyr::recode(so$SCT_clusters_0.6, "0"="0 (Patrolling MOs)", "1"="1 (HSCs)", 
    "2"="2 (CCR2+ MOs)","3"="3 (Cycling cells)","4"="4 (DC-like MOs)","5"="5 (MPP3)","6"="6 (MPP4)",
	"7"="7 (MPP2-platelets)","8"="8 (LT-HSC)","9"="9 (activated patrolling MOs)","10"="10 (activated CCR2+ MOs)",
    "11"="11 (DCs)","12"="12 (MPP2-erythrocytes)","13"="13 (MPP2-basophiles)")

p <- dittoDimPlot(so, var = "cluster_annotated", reduction = "umap_SCT")
pb <- ggplot_build(p)
color_data <- unique(pb$data[[1]][, c("group", "colour")])
original_levels <- levels(so$cluster_annotated)
colors_named <- setNames(color_data$colour, original_levels[color_data$group])

desired_order <- c("0 (Patrolling MOs)","2 (CCR2+ MOs)","4 (DC-like MOs)","9 (activated patrolling MOs)",
    "10 (activated CCR2+ MOs)","11 (DCs)","1 (HSCs)","3 (Cycling cells)","5 (MPP3)","6 (MPP4)","7 (MPP2-platelets)",
    "8 (LT-HSC)","12 (MPP2-erythrocytes)","13 (MPP2-basophiles)")
so$cluster_annotated <- factor(so$cluster_annotated, levels = desired_order)

p <- dittoDimPlot(so, var = "cluster_annotated", reduction = "umap_SCT", color.panel = colors_named)
# {png(file = paste0("A_umapSCT_annotated.png"), width = 10, height= 8, units = "in", res = 300);plot(p);dev.off()}
{pdf(file = paste0("A_umapSCT_annotated.pdf"), width = 12, height= 8);plot(p);dev.off()}

# - fig B -
# canonical markers of clusters
Idents(so) <- "SCT_clusters_0.6"
MO <- subset(so, idents = c(0,2,4,9,10,11))
MO$SCT_clusters_0.6 <- factor(MO$SCT_clusters_0.6, levels = c(0,9,4,2,10,11))
p <- DotPlot(MO, features = unique(c("Csf1r", "Gngt2", "Ear2", "Ly6c2", 
								"Ccr2", "Cd74", "H2.Aa",  
								"H2.Ab1", "H2.Eb1", "Slamf7", "Ly6a", 
								"Cd274", "Fcgr4", "Chil3", "Ifi213", 
								"Mmp8", "Cd177", "Zbtb46", "Flt3", "Itagx", "Cd209a", "Itgam")),
		group.by = "SCT_clusters_0.6", dot.scale = 5) +
		# coord_flip() +
		scale_color_gradientn(colors = c("lightblue", "white", "red")) +
		labs(x="Genes", y="Clusters", title="Markers of MO clusters") + 
		# subtitle = "avg_log2FC > 0.5, pct.1 - pct.2 > 0.5, p_val_adj < 0.05") +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(p)		
# {png("B_dotPlotHSCcanonicalMarker.png", width = 10, height = 4, units = "in", res = 300);plot(p);dev.off()}
{pdf("B_dotPlotMOcanonicalMarker.pdf", width = 10, height = 4);plot(p);dev.off()}

# - fig C -
# final_fig_APO.R
# countDEGs_pMO_cMO
MO <- subset(so, idents = c(0,2,4,9,10))
umap_coords <- Embeddings(MO, reduction = "umap_SCT")
cells_to_keep <- rownames(umap_coords)[umap_coords[, "umapSCT_1"] < 0 & umap_coords[, "umapSCT_2"] < 5]
MO <- subset(MO, cells = cells_to_keep)

MO$MO_cp <- dplyr::recode(MO$SCT_clusters_0.6, 
							"0"="pMO", "9"="pMO", 
							"2"="cMO", "4"="cMO", "10"="cMO")

MO$condition_day <- paste0(MO$condition,"_",MO$day)
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
write.csv2(allDEGs, file = "C_cMO_pMO_allDEGs.csv")

degs_filtered <- list()
degs_filtered <- lapply(degs, function(x) {x %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.75)})
# degs_filtered <- lapply(degs, function(x) {x %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5, pct.1 > 0.1 & pct.2 > 0.1)})
all_degs <- unique(unlist(lapply(degs_filtered, rownames)))

# Numbers of DEGs in cluster cMO | pMO per condition_day
summary_data <- data.frame()

for (name in names(degs_filtered)) {
  parts <- strsplit(name, "_")[[1]]
  clus <- parts[1]
  condition <- paste(parts[2], collapse = "_")
  day <- parts[3]
 
  deg_table <- degs_filtered[[name]]
  
  pos <- sum(deg_table$avg_log2FC > 0, na.rm = TRUE)
  neg <- sum(deg_table$avg_log2FC < 0, na.rm = TRUE)
  
  summary_data <- rbind(summary_data, data.frame(
    Cluster = clus,
    Condition = condition,
    Day = day,
    Upregulated = pos,
    Downregulated = neg
  ))
}
write.csv2(summary_data, file = "C_cMO_pMO_countDEGs_summary.csv")

# Melt the data for ggplot
summary_long <- pivot_longer(summary_data, 
                            cols = c(Upregulated, Downregulated),
                            names_to = "Direction",
                            values_to = "Count")

summary_long$Genes <- mapply(
  function(clus, cond, day, dir) {
    get_deg_names(clus, cond, day, dir)
  },
  summary_long$Cluster,
  summary_long$Condition,
  summary_long$Day,
  summary_long$Direction
)

summary_long$Count_signed <- with(summary_long,
  ifelse(Direction == "Upregulated", Count, -Count)
)
summary_long <- summary_long %>%
    mutate(Condition = ifelse(Condition == "MuHV4", "WT", Condition))

summary_long$Day <- factor(summary_long$Day, levels = c("d8","d38","d120"))
# summary_long$Condition <- gsub("^[0-9]_","",summary_long$Condition)
# summary_long$Cluster <- gsub("clus","",summary_long$Cluster)
summary_long$Cluster <- factor(summary_long$Cluster, levels = c("pMO","cMO"))
write_csv2(summary_long, file = "C_cMO_pMO_countDEGs_summaryLong.csv")

# Plot
p <- ggplot(summary_long, aes(x = Cluster, y = Count_signed)) +
    # First, the negative bars (downregulated)
    geom_bar(data = filter(summary_long, Count_signed < 0),
             aes(fill = "Downregulated"),
             stat = "identity",
             position = "identity",
             width = 0.7) +
    # Then, the positive bars (upregulated) on top
    geom_bar(data = filter(summary_long, Count_signed > 0),
             aes(fill = "Upregulated"),
             stat = "identity",
             position = "identity",
             width = 0.7) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    facet_grid(Condition ~ Day, scales = "free_x") +
    scale_y_continuous(
        breaks = scales::pretty_breaks(n = 8),
        expand = expansion(mult = c(0.1, 0.1))
    ) +
    scale_fill_manual(
        values = c("Upregulated" = "firebrick2", "Downregulated" = "steelblue3"),
        name = "Direction"
    ) +
    labs(
        title = "Differentially Expressed Genes per day and condition",
		subtitle = "padj < 0.05, abs(avg_log2FC) > 0.75",
        x = "Cell Cluster",
        y = "#DEGs (to mock)"
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.position = "top")
    ) +     
	coord_flip()
print(p)

# {png("C_pMO_cMO_countDEGs.1.png", width = 12, height = 10, units = "in", res = 300);plot(p);dev.off()}
# {pdf("C_pMO_cMO_countDEGs.1.pdf", width = 8, height = 10);plot(p);dev.off()}
{pdf("C_pMO_cMO_countDEGs.2.pdf", width = 12, height = 10);plot(p);dev.off()}

# - fig D - 
# volcanoHMpathway_cMO_2+4+10_MuHV4_del73_d8
# 10-2025_fig_APO.R
Idents(so) <- "SCT_clusters_0.6"
soSub <- subset(so, ident = c("2","4","10"))
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

write.csv2(rest_act, file = "D_HMpathway_results_cMO_2+4+10_MuHV4_del73_d8.csv")

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
  labs(title = "cMO (2+4+10): MuHV4 & del73 vs mock at d8", 
  subtitle = "Strong enrichment: FC > 5 & padj < 0.05 \nMild enrichment: −5 < FC < 5 & padj < 0.05 \nNot significant: padj > 0.05")
# ggsave(paste0("D_volcanoPathway_cMO_2+4+10_MuHV4_del73_d8.png"), width = 7, height = 7, units = "in", dpi = 300)
ggsave(paste0("D_volcanoHMpathway_cMO_2+4+10_MuHV4_del73_d8.pdf"), width = 7, height = 7, units = "in", dpi = 300)

# - fig E - 
# VennMuHV4_d8d38d120-Up_Log2FC_0.5
# deOutList_APO.R
# Clusters_2_4_10_Heatmap_APO
# Up & Dn-regulated genes shared genes through time and virus
# cMO = c(2,4,10)

# # To complete with deOut function !!!! from deOutList_APO.R # #

load("~/Projects/Arthur/scRNAseq_Analysis3/Clusters_2_4_10_agrr_Seurat/deOut_aggr_cluster_2_4_10.Rdata")
names(deOut) <- c("MuHV4_d8","del73_d8","MuHV4_del73_d8","MuHV4_d38","del73_d38","MuHV4_del73_d38","MuHV4_d120","del73_d120","MuHV4_del73_d120")
combined_df <- bind_rows(
  deOut[[7]] %>% mutate(comparison = "MuHV4_vs_Mock_clus2+4+10_d120"),
  deOut[[8]] %>% mutate(comparison = "del73_vs_Mock_clus2+4+10_d120"),
  deOut[[9]] %>% mutate(comparison = "MuHV4_vs_del73_clus2+4+10_d120")
)
write.csv2(combined_df, file = "E_DEGs_vs_Mock_clus2+4+10_d120.csv")

combined_df <- bind_rows(
	deOut[[1]] %>% mutate(comparison = "MuHV4_vs_Mock_clus2+4+10_d8"),
	deOut[[2]] %>% mutate(comparison = "del73_vs_Mock_clus2+4+10_d8"),
	deOut[[3]] %>% mutate(comparison = "MuHV4_vs_del73_clus2+4+10_d8"),
	deOut[[4]] %>% mutate(comparison = "MuHV4_vs_Mock_clus2+4+10_d38"),
	deOut[[5]] %>% mutate(comparison = "del73_vs_Mock_clus2+4+10_d38"),
	deOut[[6]] %>% mutate(comparison = "MuHV4_vs_del73_clus2+4+10_d38"),
	deOut[[7]] %>% mutate(comparison = "MuHV4_vs_Mock_clus2+4+10_d120"),
	deOut[[8]] %>% mutate(comparison = "del73_vs_Mock_clus2+4+10_d120"),
	deOut[[9]] %>% mutate(comparison = "MuHV4_vs_del73_clus2+4+10_d120")
)
write.csv2(combined_df, file = "E_DEGs_vs_Mock_clus2+4+10_d8d38d120.csv")

sampleID <- c("MuHV4_", "del73_")
virusDayUp <- list()
# virusDayDn <- list()

lfc = 0.5
for(virus in sampleID){
	for(d in c("d8","d38","d120")){
		virusDayUp[[paste0(virus,d)]] <- deOut[[paste0(virus,d)]] %>% dplyr::filter(p_val_adj < 0.05, avg_log2FC > lfc) %>%.$geneSymbol
		# virusDayDn[[paste0(virus,d)]] <- deOut[[paste0(virus,d)]] %>% dplyr::filter(p_val_adj < 0.05, avg_log2FC < -lfc) %>%.$geneSymbol
	}
	names(virusDayUp) <- c("d8","d38","d120")
	# names(virusDayDn) <- c("d8","d38","d120")
	totalGene <- calculate.overlap(virusDayUp) %>% sapply(., length) %>% sum
	# png(file=paste0("E_Venn_",virus,"d8d38d1200-Up_Log2FC_",lfc,".png"), width=7, height=7, units="in", res=300)
	pdf(file=paste0("E_Venn_",virus,"d8d38d120-Up_Log2FC_",lfc,".pdf"), width=7, height=7)
	print(ggVennDiagram(virusDayUp, label_alpha = 0, set_color = c("d8" = "blue","d38" ="red",'d120' = 'forestgreen'), 
						label_percent_digit=1,label_color="black") + 
						scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
						labs(title= paste0(virus,"up-genes Log2FC>", lfc,", total genes=", totalGene)) + 
						theme(plot.title = element_text(margin = margin(t = 0, b = 20), hjust = 0.5)) + NoLegend())
	dev.off()
	printOverlapGenes4(virusDayUp, file=paste0("E_",virus,"d8d38d120-Up_Log2FC_",lfc,".csv"))
	virusDayUp[1:4] <- NULL
	# totalGene <- calculate.overlap(virusDayDn) %>% sapply(., length) %>% sum
	# png(file=paste0("Venn",virus,"d8d38d120-Dn_Log2FC_",lfc,".png"), width=7, height=7, units="in", res=300)
	# print(ggVennDiagram(virusDayDn, label_alpha = 0, set_color = c("d8" = "blue","d38" ="red",'d120' = 'forestgreen'), 
						# label_percent_digit=1,label_color="black") + 
						# scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
						# labs(title= paste0(virus,"dn-genes Log2FC>", lfc,", total genes=", totalGene)) + 
						# theme(plot.title = element_text(margin = margin(t = 0, b = 20), hjust = 0.5)) + NoLegend())
	# dev.off()
	# printOverlapGenes4(virusDayDn, file=paste0(virus,"d8d38d120-Dn_Log2FC_",lfc,".csv"))
	# virusDayDn[1:4] <- NULL
}

# - fig F -
# VolcanoPlotcMo_MuHV4vsMock_D120
# cMO_2+4+10
Idents(so) <- "SCT_clusters_0.6"
soSub <- subset(so, ident = c("2","4","10"))
Idents(soSub) <- "day"
soSub <- subset(soSub, ident = "d120")

Idents(soSub) <- "condition"
markers <- FindMarkers(soSub, ident.1 = "MuHV4", ident.2 = "Mock")
markers$gene <- rownames(markers)
write.csv2(markers, file = "F_degs_cMo_MuHV4vsMock_D120.csv")

highlight_genes <- c("Il27","Ifi205","Ciita","Trem2","Cd177","Nlrc5","Ly6a","Gbp5","Gbp3","H2.Q6","Ifi47","Tap1","Irgm2",
			"Lair1","Cd274","Cd36","Cd300e","Cd74","H2.Eb1","H2.Ab1","H2.Aa","H2.Q7","AW112010","Ly6i","Cd300e","Irf1",
			"Fcgr4","Stat1","Gbp2") #,"Cxcl9","Cxcl10","Acod1"

# markers <- markers%>%filter(!str_detect(gene, 'MT-'))
markers$label_gene <- ifelse(markers$gene %in% highlight_genes, markers$gene, NA)
markers$is_highlight <- ifelse(markers$gene %in% highlight_genes, TRUE, FALSE)
markers$signif <- with(markers, ifelse(p_val_adj < 0.05, "Significant", "NS")) # & abs(avg_log2FC) > 0.5
markers$plot_pval <- ifelse(markers$p_val_adj == 0, 1e-300, markers$p_val_adj)

p <- volcanoCustom(markers, main = "cMO (2+4+10)", submain = "MuHV4 vs Mock at D120")
print(p)
# {png(paste0("F_volcanoCustom_cMo_MuHV4vsMock_D120.png"), width=10, height=10,units="in", res=300);print(p);dev.off()}
{pdf(paste0("F_volcanoCustom_cMo_MuHV4vsMock_D120.pdf"), width=10, height=10);print(p);dev.off()}

# - fig G -
# vln_SCENIC_Ird1_Stat2_Stat2_MO
# SCENIC_integrated_APO.R
MO <- subset(so, idents = c(0,2,4,9,10))
Idents(so) <- "SCT_clusters_0.6" 

optionMO <- c("MOd8_SCENIC","MOd38_SCENIC","MOd120_SCENIC")
scenicOptionsListMO <- list()
for(opt in optionMO){
scenicOptionsListMO[[opt]] <- readRDS(file = paste0("~/Projects/Arthur/scRNAseq_Analysis3/Scenic/",opt,"/int/scenicOptions.Rds"))
}

cellsAUClistMO <- list()
for(opt in optionMO){
cellsAUClistMO[[opt]] <- readRDS(file = paste0("~/Projects/Arthur/scRNAseq_Analysis3/Scenic/",opt,"/int/3.4_regulonAUC.rds"))
}

cellsAUC_matricesMO <- list()
for(opt in optionMO) {
  cellsAUC_matricesMO[[opt]] <- getAUC(cellsAUClistMO[[opt]])
}

all_regulonsMO <- unique(unlist(lapply(cellsAUC_matricesMO, rownames)))

# FUN
# Create a function to fill missing regulons with 0
fill_missing_regulons <- function(auc_matrix, all_regulons) {
  missing_regulons <- setdiff(all_regulons, rownames(auc_matrix))
  if (length(missing_regulons) > 0) {
    # Create matrix of zeros for missing regulons
    missing_matrix <- matrix(0, nrow = length(missing_regulons), ncol = ncol(auc_matrix),
                            dimnames = list(missing_regulons, colnames(auc_matrix)))
    # Combine with existing matrix
    auc_matrix <- rbind(auc_matrix, missing_matrix)
  }
  # Reorder to match all_regulons
  auc_matrix[all_regulons, ]
}
# FUN

cellsAUC_filledMO <- lapply(cellsAUC_matricesMO, fill_missing_regulons, all_regulons = all_regulonsMO)
all_regulonAUCMO <- do.call(cbind, cellsAUC_filledMO)

# FUN
clean_regulon_names <- function(regulon_names) {
  gsub("\\s*\\([0-9]+g\\)$", "", regulon_names)
}
# FUN

# FUN
# Aggregate duplicates by taking the maximum AUC value for each regulon
aggregate_regulons <- function(auc_matrix) {
  clean_names <- clean_regulon_names(rownames(auc_matrix))
  unique_regulons <- unique(clean_names)
  
  # Create aggregated matrix
  aggregated_matrix <- matrix(0, nrow = length(unique_regulons), ncol = ncol(auc_matrix),
                             dimnames = list(unique_regulons, colnames(auc_matrix)))
  
  # For each unique regulon, take max across all versions
  for(i in seq_along(unique_regulons)) {
    regulon <- unique_regulons[i]
    regulon_rows <- which(clean_names == regulon)
    
    if(length(regulon_rows) > 1) {
      # Take maximum AUC across all versions of this regulon
      aggregated_matrix[i, ] <- apply(auc_matrix[regulon_rows, , drop = FALSE], 2, max)
    } else {
      aggregated_matrix[i, ] <- auc_matrix[regulon_rows, ]
    }
  }
  
  return(aggregated_matrix)
}
# FUN

all_regulonAUC_cleanMO <- aggregate_regulons(all_regulonAUCMO)
all_regulonAUC_cleanMO <- all_regulonAUC_cleanMO[, colnames(MO)]
MO[["SCENIC"]] <- CreateAssayObject(data = all_regulonAUC_cleanMO)

auc_metaMO <- t(all_regulonAUC_cleanMO)
MO <- AddMetaData(MO, auc_metaMO)

MO$day <- factor(MO$day, levels = c("d8","d38","d120"))
MO$condition <- factor(MO$condition, levels = c("Mock", "MuHV4", "del73"))
MO <- AddMetaData(MO, auc_metaMO)
auc_metaMO <- t(all_regulonAUC_cleanMO)

DefaultAssay(MO) <- "SCENIC"
cellsMO <- colnames(MO)[MO$condition %in% c("Mock", "MuHV4", "del73")]

DefaultAssay(MO) <- "SCENIC"
p <- VlnPlot2(MO, assay = "SCENIC", pt.alpha = 0.1,
        features = c("Stat1","Stat2-extended","Irf1"), 
        group.by = "day",
		split.by = "condition",		  
        cells = cellsMO,
        stat.method = "wilcox.test", #"t.test",
              p.adjust.method = "holm", #"BH",			
        comparisons = list(c(1,2), c(1,3), c(2,3)), 
        hide.ns = FALSE, label ="p.format") + labs(title = "MO scenic score", subtitle = "clusters 0+2+4+9+10")
print(p)
# {png("G_vln_SCENIC_Ird1_Stat2_Stat2_MO.png", width = 10, height = 10, units = "in", res = 300);plot(p);dev.off()}
{pdf("G_vln_SCENIC_Ird1_Stat2_Stat2_MO.pdf", width = 10, height = 10);plot(p);dev.off()}


condition <- c("MuHV4","del73")
MO$condition_day <- gsub("^1_|^2_|^3_","",MO$condition_day)

Idents(MO) <- "condition_day"
degs <- list()
for(con in condition){
	for(d in c("d8","d38","d120")){
		Idents(MO) <- "condition_day"
		deg_key <- paste0(con, "_vs_Mock_", d)			
		deg_result <- FindMarkers(MO, ident.1 = paste0(con, "_", d), ident.2 = paste0("Mock_", d), only.pos = FALSE,
						features = c("Stat1","Stat2-extended","Irf1"), logfc.threshold = 0)
		deg_result$gene <- rownames(deg_result)
		deg_result$group <- deg_key
		degs[[deg_key]] <- deg_result
		deg_result <- FindMarkers(MO, ident.1 = paste0("MuHV4_", d), ident.2 = paste0("del73_", d), only.pos = FALSE,
						features = c("Stat1","Stat2-extended","Irf1"), logfc.threshold = 0)
		deg_result$gene <- rownames(deg_result)
		deg_result$group <- paste0("MuHV4_vs_del73_",d)
		degs[[paste0("MuHV4_vs_del73_",d)]] <- deg_result
		deg_result <- FindMarkers(MO, ident.1 = paste0("del73_", d), ident.2 = paste0("MuHV4_", d), only.pos = FALSE,
						features = c("Stat1","Stat2-extended","Irf1"), logfc.threshold = 0)
		deg_result$gene <- rownames(deg_result)
		deg_result$group <- paste0("del73_vs_MuHV43_",d)
		degs[[paste0("del73_vs_MuHV43_",d)]] <- deg_result
	}
}

combined_degs <- bind_rows(degs)
write.csv2(combined_degs, file = "Table/G_vln_SCENIC_Ird1_Stat2_Stat2_MO.csv")

# - fig H - 
# cellPropDitto_MO_
# Misc_APO.R
Idents(so) <- "SCT_clusters_0.6"
MO <- subset(so, idents = c(0,2,4,9,10))

all_colors <- dittoColors()
cluster_colors <- all_colors[1:14]
names(cluster_colors) <- 0:13 
cluster_colors
 
paint <- cluster_colors[c(1,3,5,10,11)]

for(d in c("d8","d38","d120")){
	Idents(MO) <- "day"
	soSub <- subset(MO, idents = d)
	soSub$condition <- dplyr::recode(soSub$condition, "MuHV4"="WT", "del73"="Del73")
	soSub$condition <- factor(soSub$condition, levels = c("Mock","WT","Del73"))
	p <- dittoBarPlot(soSub, "SCT_clusters_0.6", group.by = "condition", retain.factor.levels = TRUE,
						main = paste("Proportion of MOs at",d), legend.title = "MO clusters (0+2+4+9+10)", 
						color.panel = paint) + # color.panel = paint,
						theme(plot.title = element_text(size=10))
	# {png(paste0("cellPropDitto_MO_",d,".png"), width = 3, height = 5, units = "in", res = 300);plot(p);dev.off()}
	print(p)
	{pdf(paste0("H_cellPropDitto_MO_",d,".pdf"), width = 3, height = 5);plot(p);dev.off()}
}

# boxPlot_clus0+2+4+9+10
sce <- as.SingleCellExperiment(MO)
sce$cluster_col <- MO@meta.data[["SCT_clusters_0.6"]]
	
prop <- propeller(sce, clusters = sce$cluster_col, sample = sce$Sample_Name, group = sce$condition_day)
write.csv2(prop, file = paste0("H_propellerPropReplicate_MO.csv"))

propT <- getTransformedProps(clusters = sce$cluster_col, sample = sce$Sample_Name, transform="logit")
write.csv2(propT[["Proportions"]], file = paste0("H_transformedPropReplicate_MO.csv"))

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

dataProp$day <- factor(dataProp$day, levels = c("d8", "d38", "d120"))
dataProp$condition <- factor(dataProp$condition, levels = c("Mock", "MuHV4", "del73"))

day <- unique(sce$day)
clusters <- unique(sce$SCT_clusters_0.6)
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
	write.csv2(stats, file = paste0("H_propMO_stats_", d, ".csv"), row.names = FALSE)
	write.csv2(summary_data, file = paste0("H_propMO_summary_", d, ".csv"), row.names = FALSE)

	print(p)
	# {png(paste0("H_barPlot_propMO_",d,".png"), width = 10, height = 6, units = "in", res = 300);plot(p);dev.off()}
	{pdf(paste0("H_barPlot_propMO_",d,".pdf"), width = 6, height = 10);plot(p);dev.off()}
}

# - fig I -
# dotPlotcMOtopMarker
soSub <- subset(so, idents = c(2,4,10))
p <- DotPlot(soSub, features = unique(c("Cxcl10","Ifit1bl1","Chil3","Ifit1","Ifit2","Ifit3b","Ifi208","Ifit3","Gbp5","Gbp7",
										"Vcan","Oasl1","Gbp6","Fn1","Itga1","Nlrp3","Gas7","Racgap1","Cd177","Slpi","S100a4",
										"Lyz2","Mmp8","Cd300lg","Fos","Rbks","Sgms2","Rfx2","H2.Eb1","H2.Aa","H2.Ab1","Cd74",
										"Ciita","H2.DMb1","Slamf7","TremI4","Itgax","Fcgr4","Slc12a2","Eno3","Lair1","Ly6i",
										"Cd300e","Cd36","Ace","Tcf7I2","Pparg","Tgm2")),
		group.by = "SCT_clusters_0.6", dot.scale = 5) +
		# coord_flip() +
		scale_color_gradientn(colors = c("lightblue", "white", "red")) +
		labs(x="Genes", y="Clusters", title="Markers of cMO clusters") + 
		# subtitle = "avg_log2FC > 0.5, pct.1 - pct.2 > 0.5, p_val_adj < 0.05") +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
		coord_flip()
print(p)		
# {png("I_dotPlotcMOtopMarker.png", width = 10, height = 4, units = "in", res = 300);plot(p);dev.off()}
{pdf("I_dotPlotcMOtopMarker.pdf", width = 4, height = 10);plot(p);dev.off()}

# - fig J -
# heatmapAbsNaive_clus4_d120_log2FC_1
# DEGs_per_cluster_APO.R & final_fig_APO.R
load("~/Projects/Arthur/scRNAseq_Analysis3/DEGs_clusters_Seurat/deOut_aggr_cluster_4.Rdata")
names(deOut) <- c("MuHV4_d8","del73_d8","MuHV4_del73_d8","MuHV4_d38","del73_d38","MuHV4_del73_d38","MuHV4_d120","del73_d120","MuHV4_del73_d120")

clus = 4
lfc = 1

temp.1 <- deOut[[7]] %>% dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > lfc) %>%.$geneSymbol # , pct.1 - pct.2 > 0.1
temp.2 <- deOut[[8]] %>% dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > lfc) %>%.$geneSymbol # dplyr::filter(pct.1 > 0.1 & pct.2 > 0.1)
temp.3 <- deOut[[9]] %>% dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > lfc) %>%.$geneSymbol
		
geneVectorAbs <- temp.1 %>% union(temp.2) %>% union(temp.3) #%>% union(temp.4) %>% union(temp.5) %>% union(temp.6) 
gene <- c("Sqle","Cd209a","Olfm1","Fabp5","Cadps2","Clec4a2")
geneVectorAbs <- setdiff(geneVectorAbs, gene)

# sampleIDVector <- c(paste0("Mock_",d),paste0("MuHV4_",d),paste0("del73_",d))
sampleIDVector <- c("Mock","MuHV4","del73")

Idents(so) <- "SCT_clusters_0.6"
subClus <- subset(so, idents = clus)
subClus$condition_day <- gsub("^[0-9]_","", subClus$condition_day)
Idents(subClus) <- "day"
subClus <- subset(subClus, idents = "d120")

exprsMatrix <- as(subClus[["RNA"]]$data[geneVectorAbs,], "matrix")
exprsTable <- tibble(sampleID = character(0), Gene = character(0), Total = numeric(0), Positive = numeric(0), AvgExpr = numeric(0), MaxExpr = numeric(0))
for(gene1 in geneVectorAbs){
	for(sampleID1 in sampleIDVector){
		sampleBc <- colnames(subClus)[subClus$condition == sampleID1]
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
hmDataScale <- na.omit(hmDataScale)
# hmFilename <- gsub("diffExp", "heatmapTop", comparisonList[[c1]]$deFilename)
# hmFilename <- gsub("csv", "png", hmFilename)
widthEst <- 2 + length(sampleIDVector) * 0.5
heightEst <- 3 + length(geneVectorAbs) * 0.1
colTitle <- paste0(length(geneVectorAbs)," top DEGs, log2FC>",lfc)

p <- Heatmap(hmDataScale, name = "Z-score", row_km = 4, cluster_columns = FALSE, cluster_rows = TRUE, show_column_dend = FALSE, 
column_names_rot = 45, column_title = colTitle, row_names_gp = gpar(fontsize = 8), col = mr)
# , column_km = 4, col = mr
print(p)

# {png("J_heatmapAbsNaive_clus4_log2FC_1.png", width = widthEst, height = heightEst, units = "in", res = 300);plot(p);dev.off()}
{pdf("J_heatmapAbsNaive_clus4_log2FC_1.pdf", width = widthEst, height = heightEst);print(p);dev.off()}

# - fig K -
DefaultAssay(so) <- "RNA"
Idents(so) <- "SCT_clusters_0.6"
MO <- subset(so, ident = c(0,2,4,9,10))
umap_coords <- Embeddings(MO, reduction = "umap_SCT")
cells_to_keep <- rownames(umap_coords)[umap_coords[, "umapSCT_1"] < 0 & umap_coords[, "umapSCT_2"] < 5]
MO <- subset(MO, cells = cells_to_keep)

# AUC/NES scoring and dotplot for GMP & MDP
MO <- AddModuleScore(MO, features = list(c("Cd177", "Cebpb", "Chil3", "Fpr2", "Hopx", "Lcn2", "Sell", "Slpi", "Padi4")), 
						name = "NES_GMP", nbin = 50, ctrl = 100) # bin = 24, ctrl = 50
{MO$NES_GMP <- MO$NES_GMP1;MO$NES_GMP1 <- NULL}

p <- FeaturePlot(MO, features = "NES_GMP", reduction = "umap_SCT", cols = c("grey","blue")) + 
						labs(title = "GMP NES score", subtitle = "Cd177, Cebpb, Chil3, Fpr2, Hopx, Lcn2, Sell, Slpi, Padi4")
print(p)
# {png("K_featurePlot_NES_GMP.png", width = 10, height = 8, units = "in", res = 300);print(p);dev.off()}
{pdf("K_featurePlot_NES_GMP.pdf", width = 8, height = 8);print(p);dev.off()}

MO <- AddModuleScore(MO, features = list(c("Batf3", "Btla", "C1qb", "Cd209a", "Cd40", "H2.Ab1", "Hpgd", "Id3", "L1cam", "Slamf7", "Tmem176a")), 
						name = "NES_MDP", nbin = 50, ctrl = 100) # bin = 24, ctrl = 50
{MO$NES_MDP <- MO$NES_MDP1;MO$NES_MDP1 <- NULL}

p <- FeaturePlot(MO, features = "NES_MDP", reduction = "umap_SCT", cols = c("grey","blue")) + 
		labs(title = "MDP NES score", subtitle = "Batf3, Btla, C1qb, Cd209a, Cd40, H2.Ab1, Hpgd, Id3, L1cam, Slamf7, Tmem176a")
print(p)
# {png("K_featurePlot_NES_DP.png", width = 10, height = 8, units = "in", res = 300);print(p);dev.off()}
{pdf("K_featurePlot_NES_MDP.pdf", width = 8, height = 8);print(p);dev.off()}

# GMP_cells <- WhichCells(object = MO, expression = "NES_GMP" > 0)
GMP_cells <- colnames(MO)[MO$NES_GMP > 0.25]
cells <- list("GMP" = GMP_cells)
p <- Cell_Highlight_Plot(MO, cells_highlight = cells, reduction = "umap_SCT", highlight_color = c("red")) + NoLegend() +
						labs(title = "GMP NES score > 0.25", subtitle = "Cd177, Cebpb, Chil3, Fpr2, Hopx, Lcn2, Sell, Slpi, Padi4")
print(p)
# {png("K_featurePlot_NES_GMP.png", width = 10, height = 8, units = "in", res = 300);print(p);dev.off()}
{pdf("K_featurePlot_NES_GMP_notScaled.pdf", width = 8, height = 8);print(p);dev.off()}

# MDP_cells <- WhichCells(object = MO, expression = "NES_MDP" > 0)
MDP_cells <- colnames(MO)[MO$NES_MDP > 0.25]
cells <- list("MDP" = MDP_cells)
p <- Cell_Highlight_Plot(MO, cells_highlight = cells, reduction = "umap_SCT", highlight_color = c("red")) + NoLegend() +
		labs(title = "MDP NES score > 0.25", subtitle = "Batf3, Btla, C1qb, Cd209a, Cd40, H2.Ab1, Hpgd, Id3, L1cam, Slamf7, Tmem176a")
print(p)
# {png("K_featurePlot_NES_MDP.png", width = 10, height = 8, units = "in", res = 300);print(p);dev.off()}
{pdf("K_featurePlot_NES_MDP_notScaled.pdf", width = 8, height = 8);print(p);dev.off()}

# FeaturePlot(MO, features = c("NES_MDP","H2-ia-ie-AbSeq"), blend = TRUE, reduction = "umap_SCT")
# {png("featurePlot_MDP_H2-ia-ie-AbSeq.png", width = 16, height = 4, units = "in", res = 300);plot(p);dev.off()}
# {pdf("featurePlot_MDP_H2-ia-ie-AbSeq.pdf", width = 16, height = 4);plot(p);dev.off()}

# - fig L -
Idents(so) <- "SCT_clusters_0.6"
MO <- subset(so, ident = c(0,2,4,9,10))
umap_coords <- Embeddings(MO, reduction = "umap_SCT")
cells_to_keep <- rownames(umap_coords)[umap_coords[, "umapSCT_1"] < 0 & umap_coords[, "umapSCT_2"] < 5]
MO <- subset(MO, cells = cells_to_keep)

cells <- WhichCells(object = MO, expression = H2.Ab1 > 2)
cells <- list("H2.Ab1" = cells)
p <- Cell_Highlight_Plot(MO, cells_highlight = cells, reduction = "umap_SCT", highlight_color = "blue") + NoLegend() + 
					labs(title =  "H2.Ab1 gene expression > 2")
print(p)
# {png("L_featurePlot_MO_H2.Ab1.png", width = 16, height = 4, units = "in", res = 300);plot(p);dev.off()}
{pdf("L_featurePlot_MO_H2.Ab1.pdf", width = 10, height = 10);plot(p);dev.off()}

DefaultAssay(MO) <- "ADT"
cells <- WhichCells(object = MO, expression = `H2-ia-ie-AbSeq` > 1.5)
cells <- list("H2.Ab1" = cells)
p <- Cell_Highlight_Plot(MO, cells_highlight = cells, reduction = "umap_SCT", highlight_color = "red") + NoLegend() + 
					labs(title =  "H2-ia-ie-AbSeq abSeq expression > 1.5")
print(p)
# {png("L_featurePlot_MO_H2-ia-ie-AbSeq.png", width = 16, height = 4, units = "in", res = 300);plot(p);dev.off()}
{pdf("L_featurePlot_MO_H2-ia-ie-AbSeq.pdf", width = 10, height = 10);plot(p);dev.off()}

# - fig M - 
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
adt_expression <- GetAssayData(soSub, slot = "data")["H2-ia-ie-AbSeq", ]
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

# # Scatter plot of MDP and MHCI score in cluster 4
DefaultAssay(soSub) <- "RNA"
# soSub <- AddModuleScore(soSub, features = list(c("Cd177", "Cebpb", "Chil3", "Fpr2", "Hopx", "Lcn2", "Sell", "Slpi", "Padi4")), 
						# name = "NES_GMP", nbin = 50, ctrl = 100) # bin = 24, ctrl = 50
# {soSub$NES_GMP <- soSub$NES_GMP1;soSub$NES_GMP1 <- NULL}
soSub <- AddModuleScore(soSub, features = list(c("Batf3", "Btla", "C1qb", "Cd209a", "Cd40", "H2.Ab1", "Hpgd", "Id3", "L1cam", 
													"Slamf7", "Tmem176a",
													"Ccl5", "C1qb", "C1qc","Mgl2")), 
						name = "NES_MDP", nbin = 50, ctrl = 100) # bin = 24, ctrl = 50
{soSub$NES_MDP <- soSub$NES_MDP1;soSub$NES_MDP1 <- NULL}

p <- FeatureScatter(soSub, feature1 = "MHCII_RNA", feature2 = "NES_MDP", cells = NULL, group.by = "SCT_clusters_0.6", 
					split.by = NULL, log = FALSE) + labs(subtitle = "MO MHCII & MDP signature correlation") +
					geom_smooth(method = "lm", color = "red", se = FALSE)
print(p)
{pdf(file = paste0("M_scatterPlot_MO_MHCII_MDP_sig.pdf"), width = 8, height= 8);plot(p);dev.off()}

p <- FeatureScatter(soSub, feature1 = "MHCII_RNA", feature2 = "NES_MDP", cells = NULL, group.by = "SCT_clusters_0.6", 
					split.by = NULL, log = TRUE) + labs(subtitle = "MO MHCII & MDP signature correlation") #+
					# geom_smooth(method = "lm", color = "red", se = FALSE)
print(p)
{pdf(file = paste0("M_scatterPlot_MO_MHCII_MDP_log_sig.pdf"), width = 8, height= 8);plot(p);dev.off()}

# - fig N -
# volcanoCustom_MOclus0+2+4+9+10_MHCII_High_vs_low
Idents(soSub) <- "MHCII_AUC_status"
markers <- FindMarkers(soSub, ident.1 = "MHCII high", ident.2 = "MHCII low")
markers$gene <- rownames(markers)
write.csv2(markers, file = "N_DEGs_MO_0+2+4+9+10_MHCII_highVsLow.csv")
# markers <- read.csv2(file = "DEGs_MO_0+2+4+9+10_MHCII_highVsLow.csv")

highlight_genes <- c("Trem1","Cd177","F5","Sell","Fn1","Chil3","Slpi","Fpr1","Mmp8","Vcan","S100a8","Lcn2","H2.Eb1","C1qc","C1qb",
"H2.Aa","H2.Ab1","C1qa","Car4","Cd74","Ciita","Mertk","Itga9","Ddr1","Slamf7","Adgrg5","H2.DMb1","Ccl5",
"Cxcl16","Axl","Tmem176a","Mgl2","Tmem176b","H2.Eb2","Id3","Igf1","Mmp14","Creb5","Dok2","C3ar1","Ms4a7",
"Aif1","Ccr5","Traf4","H2.DMa","Cxcl9","Cxcl10","Acod1")

# markers <- markers%>%filter(!str_detect(gene, 'MT-'))
markers$label_gene <- ifelse(markers$gene %in% highlight_genes, markers$gene, NA)
markers$is_highlight <- ifelse(markers$gene %in% highlight_genes, TRUE, FALSE)
markers$signif <- with(markers, ifelse(p_val_adj < 0.05, "Significant", "NS")) # & abs(avg_log2FC) > 0.5
markers$plot_pval <- ifelse(markers$p_val_adj == 0, 1e-300, markers$p_val_adj)

p <- volcanoCustom(markers, "scRNAseq: DEG between MO MHCII High vs low cells", "clusters 0+2+4+9+10")
print(p)
# {png(paste0("N_volcanoCustom_MOclus0+2+4+9+10_MHCII_High_vs_low.png"), width=10, height=10,units="in", res=300);print(p);dev.off()}
{pdf(paste0("N_volcanoCustom_MOclus0+2+4+9+10_MHCII_High_vs_low.pdf"), width=10, height=10);print(p);dev.off()}

# - fig O -
# canonical markers of clusters
Idents(so) <- "SCT_clusters_0.6"
HSC <- subset(so, idents = c(1,3,5,6,7,8,12,13))
 
HSC <- subset(so, idents = c(1,3,5,6,8))
p <- DotPlot(HSC, features = unique(c("Hlf", "Mecom", "Mllt3", "Fgd5", # (LT-HSC) 8
								"Procr", "Mycn", "Cd34", # (ST-HSC) 1?
								"Mki67", "Pcna", "Mcm3", "Rrm2", # (cycling) 3
								"Spi1", "Irf8", "Fcgr3", "Ly6c2", #"Slp1", # (MMP3) 5
								"Il7r", "Rag2", "Dntt", "Flt3")), #(MPP4) 6
		group.by = "SCT_clusters_0.6", dot.scale = 5) +
		# coord_flip() +
		scale_color_gradientn(colors = c("lightblue", "white", "red")) +
		labs(x="Genes", y="Clusters", title="Markers of HSC clusters") + 
		# subtitle = "avg_log2FC > 0.5, pct.1 - pct.2 > 0.5, p_val_adj < 0.05") +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(p)		
# {png("O_dotPlotHSCcanonicalMarker.png", width = 10, height = 4, units = "in", res = 300);plot(p);dev.off()}
{pdf("O_dotPlotHSCcanonicalMarker.pdf", width = 10, height = 4);plot(p);dev.off()}

# - fig P -
# HSC 1+3+5+8 d8 d38 d120 virus vs mock
HSC <- subset(so, idents = c(1,3,5,8)) #6?
# umap_coords <- Embeddings(HSC, reduction = "umap_SCT")
# cells_to_keep <- rownames(umap_coords)[umap_coords[, "umapSCT_1"] > 0 & umap_coords[, "umapSCT_2"] < 6]
# HSC <- subset(HSC, cells = cells_to_keep)

HSC$condition_day <- paste0(HSC$condition,"_",HSC$day)
degs <- list()
condition <- c("MuHV4","del73")

for(con in condition){
	for(d in c("d8","d38","d120")){
		Idents(HSC) <- "condition_day"
		deg_key <- paste0(con, "_", d)			
		deg_result <- FindMarkers(HSC, ident.1 = paste0(con, "_", d), ident.2 = paste0("Mock_", d), only.pos = FALSE)
		deg_result$gene <- rownames(deg_result)
		deg_result$group <- deg_key
		degs[[deg_key]] <- deg_result
	}
}

allDEGs <- bind_rows(degs)
write.csv2(allDEGs, file = "P_HSC_clus1+3+5+8_virus_vs_mock_allDEGs.csv")

degs_filtered <- list()
degs_filtered <- lapply(degs, function(x) {x %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)})
# degs_filtered <- lapply(degs, function(x) {x %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5, pct.1 > 0.1 & pct.2 > 0.1)})
all_degs <- unique(unlist(lapply(degs_filtered, rownames)))

# Numbers of DEGs in HSCs per condition_day
summary_data <- data.frame()

for (name in names(degs_filtered)) {
  parts <- strsplit(name, "_")[[1]]
  # clus <- parts[1]
  condition <- paste(parts[1], collapse = "_")
  day <- parts[2]
 
  deg_table <- degs_filtered[[name]]
  
  pos <- sum(deg_table$avg_log2FC > 0, na.rm = TRUE)
  neg <- sum(deg_table$avg_log2FC < 0, na.rm = TRUE)
  
  summary_data <- rbind(summary_data, data.frame(
    # Cluster = clus,
    Condition = condition,
    Day = day,
    Upregulated = pos,
    Downregulated = neg
  ))
}
# write.csv2(summary_data, file = "P_HSC_clus1+3+5+8_virus_vs_mock_countsDEGs_summary.csv")

# Melt the data for ggplot
summary_long <- pivot_longer(summary_data, 
                            cols = c(Upregulated, Downregulated),
                            names_to = "Direction",
                            values_to = "Count")

summary_long$Genes <- mapply(
  function(cond, day, dir) { # clus,
    get_names_conditions_day(cond, day, dir) # clus,
  },
  # summary_long$Cluster,
  summary_long$Condition,
  summary_long$Day,
  summary_long$Direction
)

summary_long$Count_signed <- with(summary_long,
  ifelse(Direction == "Upregulated", Count, -Count)
)
summary_long <- summary_long %>%
    mutate(Condition = ifelse(Condition == "MuHV4", "WT", Condition))

summary_long$Day <- factor(summary_long$Day, levels = c("d8","d38","d120"))
# summary_long$Condition <- gsub("^[0-9]_","",summary_long$Condition)
# summary_long$Cluster <- gsub("clus","",summary_long$Cluster)
# summary_long$Cluster <- factor(summary_long$Cluster, levels = c("pHSC","cHSC"))
# write_csv2(summary_long, file = "P_HSC_clus1+3+5+8_virus_vs_mock_countsDEGs_summaryLong.csv")

# Plot
p <- ggplot(summary_long, aes(x = Day, y = Count_signed)) + # Cluster
    # First, the negative bars (downregulated)
    geom_bar(data = filter(summary_long, Count_signed < 0),
             aes(fill = "Downregulated"),
             stat = "identity",
             position = "identity",
             width = 0.7) +
    # Then, the positive bars (upregulated) on top
    geom_bar(data = filter(summary_long, Count_signed > 0),
             aes(fill = "Upregulated"),
             stat = "identity",
             position = "identity",
             width = 0.7) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    facet_grid(Condition ~ Day, scales = "free_x") +
    scale_y_continuous(
        breaks = scales::pretty_breaks(n = 8),
        expand = expansion(mult = c(0.1, 0.1))
    ) +
    scale_fill_manual(
        values = c("Upregulated" = "firebrick2", "Downregulated" = "steelblue3"),
        name = "Direction"
    ) +
    labs(
        title = "Differentially Expressed Genes per day and condition\nHSC clusters (1+3+5+8)",
		subtitle = "padj < 0.05, abs(avg_log2FC) > 0.75",
        x = "Day",
        y = "#DEGs (to mock)"
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.position = "top")
    # ) +     
	# coord_flip()
print(p)

# {png("P_HSC_clus1+3+5+8_virus_vs_mock_countsDEGs.1.png", width = 12, height = 10, units = "in", res = 300);plot(p);dev.off()}
{pdf("P_HSC_clus1+3+5+8_virus_vs_mock_countsDEGs.1.pdf", width = 6, height = 10);plot(p);dev.off()}
# {pdf("P_HSC_clus1+3+5+8_virus_vs_mock_countsDEGs.2.pdf", width = 12, height = 10);plot(p);dev.off()}

# - fig Q -
# HSC 1+3+5+8 d8 virus vs mock 
HSC <- subset(so, idents = c(1,3,5,8))
# umap_coords <- Embeddings(HSC, reduction = "umap_SCT")
# cells_to_keep <- rownames(umap_coords)[umap_coords[, "umapSCT_1"] > 0 & umap_coords[, "umapSCT_2"] < 6]
# HSC <- subset(HSC, cells = cells_to_keep)

HSC$condition_day <- paste0(HSC$condition,"_",HSC$day)
degs <- list()
condition <- c("MuHV4","del73")

Idents(HSC) <- "day"
soSub <- subset(HSC, ident = "d8")

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
save(resList, file = "Q_resListHM_HSC_clus1+3+5+8_d8.Rdata")
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

# filter(grepl(pattern = "Hallmark_interferon", ignore.case = T, x = Pathway))

# Create a separate column to identify the aa_path points
# aa_path <- rest_act %>% filter(FC > 0.5 & adjPval < 0.05) %>% arrange(desc(FC)) %>% head(,n=5)
aa_path <- rest_act %>% filter(Pathway%in% c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
								"HALLMARK_ALLOGRAPH_REJECTION","HALLMARK_MTORC1_SIGNALING","HALLMARK_GLYCOLYSIS"))
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

write.csv2(rest_act, file = "Table/Q_HMpathway_results_HSC_clus1+3+5+8_d8.csv")

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
  labs(title = "HSC (1+3+5+8): MuHV4 & del73 vs mock at d8", 
  subtitle = "Strong enrichment: FC > 5 & padj < 0.05 \nMild enrichment: −5 < FC < 5 & padj < 0.05 \nNot significant: padj > 0.05")
# ggsave(paste0("Q_volcanoPathway_HSC_clus1+3+5+8_d8.png"), width = 7, height = 7, units = "in", dpi = 300)
ggsave(paste0("Q_volcanoHMpathway_HSC_clus1+3+5+8_d8.pdf"), width = 7, height = 7, units = "in", dpi = 300)

# - fig R -
# HSC 1+3+5+8 d120 MuHV4 vs mock 
HSC <- subset(so, idents = c(1,3,5,8))
# umap_coords <- Embeddings(HSC, reduction = "umap_SCT")
# cells_to_keep <- rownames(umap_coords)[umap_coords[, "umapSCT_1"] > 0 & umap_coords[, "umapSCT_2"] < 6]
# HSC <- subset(HSC, cells = cells_to_keep)

HSC$condition_day <- paste0(HSC$condition,"_",HSC$day)
degs <- list()
condition <- c("MuHV4","del73")

Idents(HSC) <- "day"
soSub <- subset(HSC, ident = "d120")

Idents(soSub) <- "condition"
markers <- FindMarkers(soSub, ident.1 = "MuHV4", ident.2 = "Mock")
markers$gene <- rownames(markers)
write.csv2(markers, file = "R_DEGs_HSC_clus1+3+5+8_MuHV4vsMock_d120.csv")
# markers <- read.csv2(file = "R_DEGs_HSC_clus1+3+5+8_MuHV4vsMock_d120.csv")

highlight_genes <- c("Nlrc5","Stat1","Igtp","H2.T22","Fgl2","Junb","Fosb","Fos","H2.Q7","H2.Eb1","H2.Ab1","H2.Aa")

# markers <- markers%>%filter(!str_detect(gene, 'MT-'))
markers$label_gene <- ifelse(markers$gene %in% highlight_genes, markers$gene, NA)
markers$is_highlight <- ifelse(markers$gene %in% highlight_genes, TRUE, FALSE)
markers$signif <- with(markers, ifelse(p_val_adj < 0.05, "Significant", "NS")) # & abs(avg_log2FC) > 0.5
markers$plot_pval <- ifelse(markers$p_val_adj == 0, 1e-300, markers$p_val_adj)

p <- volcanoCustom(markers, "scRNAseq: DEG between HSC MuHV4 vs Mock at D120", "clusters 1+3+5+8")
print(p)
# {png(paste0("R_volcanoCustom_HSC_clus1+3+5+8_MuHV4vsMock_d120.png"), width=10, height=10,units="in", res=300);print(p);dev.off()}
{pdf(paste0("R_volcanoCustom_HSC_clus1+3+5+8_MuHV4vsMock_d120.pdf"), width=10, height=10);print(p);dev.off()}

# - fig S -
# HSC cluster proportion (1,3,5,6,7,8,12,13)
# cellPropDitto_HSC_
# Misc_APO.R
Idents(so) <- "SCT_clusters_0.6"
HSC <- subset(so, idents = c(1,3,5,6,7,8,12,13))
all_colors <- dittoColors()
cluster_colors <- all_colors[1:14]
names(cluster_colors) <- 0:13 
cluster_colors
paint <- cluster_colors[c(2,4,6,7,8,9,13,14)]

# # Cell proportions
for(d in c("d8","d38","d120")){
	Idents(HSC) <- "day"
	soSub <- subset(HSC, idents = d)
	
	p <- dittoBarPlot(soSub, "SCT_clusters_0.6", group.by = "condition_day", retain.factor.levels = TRUE,
						main = paste("Proportion of HSCs at",d), legend.title = "Mo clusters",
						color.panel = paint) + # , color.panel = paint
						theme(plot.title = element_text(size=10))
	print(p)
	# {png(file =paste0("S_cellPropDitto_HSC_",d,".png"), width = 3, height = 5, units = "in", res = 300);plot(p);dev.off()}
	{pdf(file = paste0("S_cellPropDitto_HSC_",d,".pdf"), width = 3, height = 5);plot(p);dev.off()}
}

# boxPlot_clus1+3+5+5+6+7+8+12+13
sce <- as.SingleCellExperiment(HSC)
sce$cluster_col <- HSC@meta.data[["SCT_clusters_0.6"]]
	
prop <- propeller(sce, clusters = sce$cluster_col, sample = sce$Sample_Name, group = sce$condition_day)
write.csv2(prop, file = paste0("S_propellerPropReplicate_HSC.csv"))

propT <- getTransformedProps(clusters = sce$cluster_col, sample = sce$Sample_Name, transform="logit")
write.csv2(propT[["Proportions"]], file = paste0("S_transformedPropReplicate_HSC.csv"))

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

dataProp$day <- factor(dataProp$day, levels = c("d8", "d38", "d120"))
dataProp$condition <- factor(dataProp$condition, levels = c("Mock", "MuHV4", "del73"))

day <- unique(sce$day)
clusters <- unique(sce$SCT_clusters_0.6)
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
															p.adjust.method = "none", #"bonferroni", 
															ref.group = "Mock",
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
				palette = condition_colors, position = position_dodge(0.8)) + labs(title = paste0("HSC at ",d), 
																				subtitle = "clusters 1+3+4+5+6+7+8+12+13")
	
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
	write.csv2(stats, file = paste0("S_propHSC_stats_", d, ".csv"), row.names = FALSE)
	write.csv2(summary_data, file = paste0("S_propHSC_summary_", d, ".csv"), row.names = FALSE)

	print(p)
	# {png(paste0("S_barPlot_propHSC_",d,".png"), width = 10, height = 6, units = "in", res = 300);plot(p);dev.off()}
	{pdf(paste0("S_barPlot_propHSC_",d,".pdf"), width = 6, height = 10);plot(p);dev.off()}
}

# - fig T -
# waterfallCustom_MuHV4vsMock_1_HSCsub_d8
# Scenic_HSCsub_d8_APO.R 
Idents(so) <- "SCT_clusters_0.6"
HSC <- subset(so, idents = c(1,3,5,8))

ID = "HSCsub_d8"
day = "d8"
Idents(HSC) <- "day"
soSub <- subset(HSC, ident = day)

setwd(paste0("~/Projects/Arthur/scRNAseq_Analysis3/Scenic/",ID,"_SCENIC/"))
scenicOptions <- readRDS(file = "~/Projects/Arthur/scRNAseq_Analysis3/Scenic/HSCsub_d8_SCENIC/int/scenicOptions.Rds")
regulonAUC <- readRDS("~/Projects/Arthur/scRNAseq_Analysis3/Scenic/HSCsub_d8_SCENIC/int/3.4_regulonAUC.rds")
aucMatrix <- as.matrix(getAUC(regulonAUC))
common_cells <- intersect(colnames(soSub), colnames(aucMatrix))
aucMatrix <- aucMatrix[, common_cells]
soSub[[paste0(ID,"_SCENIC")]] <- CreateAssayObject(counts = aucMatrix)
auc_meta <- t(aucMatrix)
soSub <- AddMetaData(soSub, auc_meta)

DefaultAssay(soSub) <- paste0(ID,"_SCENIC")
soSub <- subset(soSub, cells = WhichCells(soSub, expression = umapSCT_1 >= 0))
soSub <- subset(soSub, cells = WhichCells(soSub, expression = umapSCT_2 < 5))

regulons <- loadInt(scenicOptions, "aucell_regulons")

# regulons[grepl("Stat2",names(regulons))]

setwd("~/Projects/Arthur/scRNAseq_Analysis3/Article/Figures")
soSub@meta.data <- soSub@meta.data %>% dplyr::rename("Irf3 (20g)"="Irf3_extended (20g)","Stat2 (31g)"="Stat2_extended (31g)",
									"Ikzf2 (37g)"="Ikzf2_extended (37g)","Myb (98g)"="Myb_extended (98g)")

regulons <- c("Irf1 (27g)","Cebpd (33g)","Stat1 (61g)","Stat2 (31g)","Irf9 (31g)","Cebpa (69g)","Maf (16g)","Myb (98g)","Myc (80g)",
				"Sp1 (88g)","Klf10 (55g)","Fosb (20g)","Junb (17g)")   

Idents(soSub) <- "condition"
p <- WaterfallPlot(soSub,features = regulons, ident.1 = "MuHV4", ident.2 = "Mock", exp.transform = FALSE, 
					log.base = 2, angle = 45, title = "HSC: MuHV4 vs Mock at d8\nclusters 1+3+5+8") 
					# color_theme = "D", # center_color = FALSE, # top.n = 20
print(p)
# {png(paste0("T_waterfallCustom_MuHV4vsMock_",ID,".png"), width = 10, height = 6, units = "in", res = 300);plot(p);dev.off()}
{pdf(paste0("T_waterfallCustom_MuHV4vsMock_",ID,".pdf"), width = 10, height = 6);plot(p);dev.off()}

ID = "HSCsub_d120"
day = "d120"
Idents(HSC) <- "day"
soSub <- subset(HSC, ident = day)

setwd(paste0("~/Projects/Arthur/scRNAseq_Analysis3/Scenic/",ID,"_SCENIC/"))
scenicOptions <- readRDS(file = "~/Projects/Arthur/scRNAseq_Analysis3/Scenic/HSCsub_d120_SCENIC/int/scenicOptions.Rds")
regulonAUC <- readRDS("~/Projects/Arthur/scRNAseq_Analysis3/Scenic/HSCsub_d120_SCENIC/int/3.4_regulonAUC.rds")
aucMatrix <- as.matrix(getAUC(regulonAUC))
common_cells <- intersect(colnames(soSub), colnames(aucMatrix))
aucMatrix <- aucMatrix[, common_cells]
soSub[[paste0(ID,"_SCENIC")]] <- CreateAssayObject(counts = aucMatrix)
auc_meta <- t(aucMatrix)
soSub <- AddMetaData(soSub, auc_meta)

DefaultAssay(soSub) <- paste0(ID,"_SCENIC")
soSub <- subset(soSub, cells = WhichCells(soSub, expression = umapSCT_1 >= 0))

regulons <- loadInt(scenicOptions, "aucell_regulons")
regulons[grepl("Stat1",names(regulons))]

setwd("~/Projects/Arthur/scRNAseq_Analysis3/Article/")
soSub@meta.data <- soSub@meta.data %>% dplyr::rename("Irf3 (16g)"="Irf3_extended (16g)", "Maf (43g)"="Maf_extended (43g)",
							"Stat2 (36g)"="Stat2_extended (36g)", "Ikzf2 (52g)"="Ikzf2_extended (52g)",
							"Fosb (77g)"="Fosb_extended (77g)", "Klf10 (15g)"="Klf10_extended (15g)") 

regulons <- c("Irf1 (27g)","Stat1 (33g)","Maf (43g)","Hoxb2 (15g)","Irf9 (25g)","Irf3 (16g)","Stat2 (36g)","Myb (17g)",
	"Cebpa (10g)", "Myc (73g)","Sp1 (91g)","Ikzf2 (52g)","Cebpd (121g)","Rel (22g)","Fosb (77g)","Klf10 (15g)","Junb (19g)")  

Idents(soSub) <- "condition"
p <- WaterfallPlot(soSub, features = regulons, ident.1 = "MuHV4", ident.2 = "Mock", exp.transform = FALSE, 
					log.base = 2, angle = 45, title = "HSC: MuHV4 vs Mock at d120\nclusters 1+3+5+8") 
					# color_theme = "D", # center_color = FALSE, # top.n = 20
print(p)
# {png(paste0("T_waterfallCustom_MuHV4vsMock_",ID,".png"), width = 10, height = 6, units = "in", res = 300);plot(p);dev.off()}
{pdf(paste0("T_waterfallCustom_MuHV4vsMock_",ID,".pdf"), width = 10, height = 6);plot(p);dev.off()}

# - fig U -
# dotPlot_MHCII_sig
# 10-2025_fig_APO.R
Idents(so) <- "SCT_clusters_0.6"
soSub <- subset(so, ident = c("1","3","5","8"))

cellList <- list()
cellList[["MHCII_scores"]] <- c("H2-ia-ie-AbSeq", "H2.Aa", "H2.Eb1", "H2.Ab1","Ciita","H2.DMb1", "AW112010") #"H2.DMb2",

DefaultAssay(so) <- "SCT"
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

# png(paste0("U_HSC_clus1+3+5+8_MHCIIsig_cellAssignAuto.png"), width = 6, height = 6, units = "in", res = 300)
pdf(paste0("U_HSC_clus1+3+5+8_MHCIIsig_cellAssignAuto.pdf"), width = 12, height = 4)   
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
dev.off()
cells_assignment[["MHCII_scores"]][["aucThr"]]
# $selected
# Global_k1 
# 0.1079959 

# $thresholds
                 # threshold nCells
# tenPercentOfMax 0.08749857   5359
# Global_k1       0.10799587   3960

# $comment
# [1] ""

# png(paste0("U_HSC_clus1+3+5+8_MHCIIsig_cellAssignHistCutOff.png"), width = 6, height = 6, units = "in", res = 300)
pdf(paste0("U_HSC_clus1+3+5+8_MHCIIsig_cellAssignHistCutOff.pdf"), width = 12, height = 4)   
AUCell_plotHist(cells_AUC["MHCII_scores",], aucThr=c(0.01))
abline(v=c(0.01))
dev.off()

threshold = 0.01
MHCII_class <- ifelse(MHCII_scores > threshold, "MHCII high", "MHCII low")
table(MHCII_class)
# MHCII_class
# MHCII_high  MHCII_low 
     # 10538      11975 

soSub$MHCII_AUC_status <- MHCII_class
soSub$MHCII_AUC_score <- as.numeric(MHCII_scores[colnames(soSub)])

soTemp <- subset(soSub, cells = WhichCells(soSub, expression = umapSCT_1 >= 0))
soTemp <- subset(soTemp, cells = WhichCells(soTemp, expression = umapSCT_2 < 5))
p <- FeaturePlot(soTemp, reduction = "umap_SCT", features  = cellList[["MHCII_scores"]], label = TRUE, ncol = 3)
print(p)
# {png(paste("U_featuresPlot_HSC_clus1+3+5+8_MHCII.png"), width=20, height=10, units="in", res=300);plot(p);dev.off()}
{pdf(paste("U_featuresPlot_HSC_clus1+3+5+8_MHCII.pdf"), width=20, height=10);plot(p);dev.off()}

p <- DimPlot(soTemp, group.by = "MHCII_AUC_status", label = FALSE, reduction = "umap_SCT") +
	scale_color_manual(values = c("MHCII_high" = "red","MHCII_low" = "gray")) +
	 labs(title = "HSC MHCII signature without AbSeq integration")
print(p)
# {png(file = paste0("M_umap_HSC_clus1+3+5+8_MHCII.png"), width = 6, height= 3, units = "in", res = 300);plot(p);dev.off()}
{pdf(file = paste0("M_umap_HSC_clus1+3+5+8_MHCII.pdf"), width = 8, height= 8);plot(p);dev.off()}
rm(soTemp)

p <- DotPlot(soSub, assay = "RNA", features = cellList[["MHCII_scores"]], 
        group.by = "MHCII_AUC_status", scale = TRUE) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        legend.box = "horizontal",           # Horizontal arrangement
        legend.box.just = "center",          # Center alignment
        legend.spacing.x = unit(0.5, 'cm')) + # Spacing between legends
    guides(color = guide_colorbar(title = "Avg Expression", order = 1),
        size = guide_legend(title = "Pct Exp", order = 2)) +
	labs(title = "HSC MHCII cell signature without AbSeq", y = "MHCII cells")	+
	theme(plot.title = element_text(size = 10))	 
p[["data"]][["features.plot"]][is.na(p[["data"]][["features.plot"]])] <- "H2-ia-ie-AbSeq"
print(p)
# {png(file = paste0("U_dotPlot_HSC_clus1+3+5+8_MHCII_sig.png"), width = 6, height= 3, units = "in", res = 300);plot(p);dev.off()}
{pdf(file = paste0("U_dotPlot_HSC_clus1+3+5+8_MHCII_sig.pdf"), width = 6, height= 3);plot(p);dev.off()}

# # Integrate ADT with RNA scores
soSub$MHCII_RNA <- soSub$MHCII_AUC_score

DefaultAssay(soSub) <- "ADT"
adt_expression <- GetAssayData(soSub, slot = "data")["H2-ia-ie-AbSeq", ]
hist(adt_expression)
# quantile(adt_expression, 0.75) !!!!!!!!!!!!need to improve the selection
     # 75% 
# 0.805331 
soSub$MHCII_ADT_binary <- adt_expression > 0.4 # median(adt_expression)

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

soTemp <- subset(soSub, cells = WhichCells(soSub, expression = umapSCT_1 >= 0))
soTemp <- subset(soTemp, cells = WhichCells(soTemp, expression = umapSCT_2 < 5))

p <- DimPlot(soTemp, group.by = "MHCII_AbSeq_integration", label = FALSE, reduction = "umap_SCT") +
	scale_color_manual(values = c("MHCII high" = "red", "MHCII low" = "gray")) +
	 labs(title = "HSC MHCII signature with AbSeq integration")
print(p)
# {png(file = paste0("M_umap_HSC_MHCII_sig2.png"), width = 6, height= 3, units = "in", res = 300);plot(p);dev.off()}
{pdf(file = paste0("M_umap_HSC_MHCII_sig2.pdf"), width = 8, height= 8);plot(p);dev.off()}
rm(soTemp)

# dotPlot_MO_MHCII_sig2
p <- DotPlot(soSub, assay = "RNA", features = cellList[["MHCII_scores"]], group.by = "MHCII_AbSeq_integration", scale = TRUE) +
	 theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        legend.box = "horizontal",           # Horizontal arrangement
        legend.box.just = "center",          # Center alignment
        legend.spacing.x = unit(0.5, 'cm')) + # Spacing between legends
    guides(color = guide_colorbar(title = "Avg Expression", order = 1),
        size = guide_legend(title = "Pct Exp", order = 2)) +
	labs(title = "HSC MHCII cell signature with AbSeq", y = "MHCII cells")	+
	theme(plot.title = element_text(size = 10))	 
p[["data"]][["features.plot"]][is.na(p[["data"]][["features.plot"]])] <- "H2-ia-ie-AbSeq"
print(p)
# {png(file = paste0("M_dotPlot_HSC_MHCII_sig2.png"), width = 6, height= 3, units = "in", res = 300);plot(p);dev.off()}
{pdf(file = paste0("M_dotPlot_HSC_MHCII_sig2.pdf"), width = 6, height= 4);plot(p);dev.off()}

# - fig V -
# boxPlot_MHCII_low_clus1+3+5+8
# 10-2025_fig_APO.R
soSub$day <- factor(soSub$day, levels = c("d8","d38","d120"))
p <- dittoBarPlot(soSub, "MHCII_AUC_status", group.by = "condition", split.by = "day", sub = "cluster 1+3+5+8")
# {png(paste0("V_cellProp_MHCII_clus1+3+5+8.png"), width=6, height=4, units = "in", res = 300);plot(p);dev.off()}
{pdf(paste0("V_cellProp_MHCII_clus1+3+5+8.pdf"), width=6, height=4);plot(p);dev.off()}

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
    write.csv2(stats, file=paste0("test_V_boxPlot_",cell,"_clus1+3+5+8_stats.csv"))

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
    pdf(paste0("V_boxPlot_",cell,"_clus1+3+5+8.pdf"), width = 5, height = 7)
    print(b)
    dev.off()
}

# - fig W -
# volcanoCustom_HSCclus1+3+5+8_MHCII_High_vs_low
# 10-2025_fig_APO.R
Idents(soSub) <- "MHCII_AUC_status"
markers <- FindMarkers(soSub, ident.1 = "MHCII high", ident.2 = "MHCII low")
markers$gene <- rownames(markers)
write.csv2(markers, file = "W_DEGs_HSCs1+3+5+8_MHCII_highVsLow.csv")
# markers <- read.csv2(file = "W_DEGs_HSCs1+3+5+8_MHCII_highVsLow.csv")

# significant_genes <- markers %>% filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.05)
# top_up <- significant_genes %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 20)
# top_down <- significant_genes %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC) %>% slice_head(n = 20)
# highlight_genes <- bind_rows(top_up, top_down) %>% pull(gene)

highlight_genes <- c("AW112010","H2.Eb1","H2.Aa","H2.Ab1","Iigp1","Cd74","Ciita","Gbp4","Gbp6","Gbp8","Ifit1","Ifit3","Ifit1bL1","Ifit3b",
"Clec14a","Tm4sf1","Cyp26b1","Ifi208","Slamf7","Mpo","Slpi","Elane","Ms4a3","Clec2a","Olr1","Tmem119","Tmem178")          

# markers <- markers%>%filter(!str_detect(gene, 'MT-'))
markers$label_gene <- ifelse(markers$gene %in% highlight_genes, markers$gene, NA)
markers$is_highlight <- ifelse(markers$gene %in% highlight_genes, TRUE, FALSE)
markers$signif <- with(markers, ifelse(p_val_adj < 0.05, "Significant", "NS")) # & abs(avg_log2FC) > 0.5
markers$plot_pval <- ifelse(markers$p_val_adj == 0, 1e-300, markers$p_val_adj)

p <- volcanoCustom(markers, "scRNAseq: DEG between HSC MHCII high vs low", "clusters 1+3+5+8")
print(p)
# {png(paste0("W_volcanoCustom_HSC_clus1+3+5+8_MHCII_highVsLow.png"), width=10, height=10,units="in", res=300);print(p);dev.off()}
{pdf(paste0("W_volcanoCustom_HSC_clus1+3+5+8_MHCII_highVsLow.pdf"), width=10, height=10);print(p);dev.off()}

# ^^^^^^^^^
# Save data

write_rds(so, file="so.rds")
gc()

# library(loupeR)

# DefaultAssay(so) <- "RNA"
# create_loupe_from_seurat(so, output_name = "loupeRNA") #, output_dir = "../LoupeData", force = T)

# DefaultAssay(so) <- "ADT"
# create_loupe_from_seurat(so, output_name = "loupeADT") #, output_dir = "../LoupeData", force = T)

# keep(so, sure = T)
gc()
# ^^^^^^^^^

# ********
# Figure 3 
# ********
# bulk RNAseq; RNAseq_mono_EPV_vs_Mock

setwd("~/Projects/Arthur/scRNAseq_Analysis3/Article/Figures")
# - fig S3C -
# # ciider_MA_ plot MuHV-4 versus Mock D8
tf <- read_excel("../Table/3C_ciiider_MuHV4_D8_vs_Mock_D8_MostSigDeficit.xlsx")

tf$`Transcription Factor Name` <- gsub("\\(var,2\\)","",tf$`Transcription Factor Name`)
tf$`Transcription Factor Name` <- gsub("\\(var,3\\)","",tf$`Transcription Factor Name`)
TF<- c("IRF1","IRF7","IRF2","STAT1::STAT2","Stat2","IRF8","IRF3","CEBPB","CEBPE","CEBPA","CEBPD","CEBPG",
		"IRF4","JUN::JUNB","FOSL1::JUN","FOSL1::JUND","FOS::JUNB","FOS::JUND","FOSL1::JUNB","FOS","FOSL2::JUND",
		"FOSB::JUNB","IRF9","FOSL2","HOXB3","IRF5","ATF4","NFIL3","BACH2","STAT1","REL","NFKB2","KLF2","KLF11","KLF10",
		"EGR1","Plagl1","TFAP2C","TFAP2B")

colnames(tf) <- gsub(" ", ".", colnames(tf))
tf_highlight <- tf[tf$Transcription.Factor.Name %in% TF, ]
# tf_highlight$Significance <- ifelse(tf_highlight$Site.P.Value < 0.05, "Significant", "Not significant")
tf_highlight$Significance <- ifelse(tf_highlight$Significance.Score > 1.30103 | tf_highlight$Significance.Score < -1.30103, "Significant", "Not significant")

# tf$Significance <- ifelse(tf$Site.P.Value < 0.05, "Significant", "Not significant")
tf$Significance <- ifelse(tf$Significance.Score > 1.30103 | tf$Significance.Score < -1.30103, "Significant", "Not significant")
tf$Deficit_Level <- ifelse(tf$Log2.Enrichment < 0, "Deficit", "Enriched")

# Add highlight flag
tf$Highlight <- tf$Transcription.Factor.Name %in% TF

# Create the plot with size mapped to significance score
p <- ggplot(tf, aes(x = Average.Log2.Proportion.Bound, 
               y = Log2.Enrichment)) +
    # Main points with size based on absolute significance score
    geom_point(aes(color = Significance, 
                   alpha = abs(Significance.Score), #Deficit_Level,
                   size = abs(Significance.Score)),  # Use absolute value
               shape = 16) +
    # Highlighted TFs as hollow circles with size based on their significance
    geom_point(data = tf_highlight, #subset(tf, Highlight),
               aes(size = abs(Significance.Score),
               color = Significance)) + #,"#E98415", #"#E63946",
               # shape = 1,  # Hollow circles
               # stroke = 1.2) +  # Thicker border for hollow circles
    # Labels for highlighted TFs
    geom_text_repel(data = subset(tf, Highlight),
                    aes(label = Transcription.Factor.Name),
                    size = 2.8,
                    max.overlaps = 20,
                    box.padding = 0.5,
                    segment.color = "gray50",
                    segment.alpha = 0.5) +
    # Horizontal line at y=0
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    # Color scale
    scale_color_manual(values = c("Not significant" = "#457B9D", 
                                  "Significant" = "#E63946")) +
    # Alpha scale
    # scale_alpha_manual(values = c("Deficit" = 0.9, "Enriched" = 0.5)) +
	scale_alpha_continuous(name = "|Log2 Enrichment|",
                         range = c(0.3, 1),  # Min and max transparency
                         breaks = seq(0, max(abs(tf$Log2.Enrichment)), 
                                      length.out = 5)) +
    # Size scale - adjust range as needed
    scale_size_continuous(name = "Significance Score\n(absolute value)",
                          range = c(0.5, 4),  # Min and max point sizes
                          breaks = seq(0, max(abs(tf$Significance.Score)), 
                                       length.out = 5)) +
    # Labels
    labs(x = "Average.Log2.Proportion.Bound",
         y = "Log2.Enrichment",
         title = "Transcription Factor Enrichment Plot",
         subtitle = "MuHV4 vs Mock at d8") +
    # Theme
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right",
          legend.box = "vertical")
print(p)
# {png("3C_ciider_MA_plot_MuHV4_vs_Mock_D8.png", width = 10, height = 8, units = "in", res = 300);print(p);dev.off()}
{pdf("3C_ciider_MA_plot_MuHV4_vs_Mock_D8.pdf", width = 10, height = 8);print(p);dev.off()}

# - fig S3D -
tf <- read.csv("../Table/3D_ciiider_MuHV4_D120_vs_Mock_D120_MostSigDeficit.csv")

TF <- c("RORA","Rela","Ebf2","NFIB","TCF21","Bach1::Mafk","BACH2","JUN::JUNB","XBP1","Sp1","Irf7","Irf3","Cred5","Irf5","Irf1",
	"Nrf1","Crebpb","Cebpe","Irf8","Atf4","Sox21","Nfil3","Rarb","Irf9")

TF <- c("RORA","RELA","Ebf2","RORC","NFIB","TCF21","Bach1::Mafk","BACH2","JUN::JUNB","XBP1","SP1","IRF7","IRF3","Creb5","IRF5","Irf1",
"Nrf1","IRF4","CEBPB","CEBPE","IRF8","ATF4","SOX21","NFIL3","RARB","IRF9")

tf_highlight <- tf[tf$Transcription.Factor.Name %in% TF, ]
# tf_highlight$Significance <- ifelse(tf_highlight$Site.P.Value < 0.05, "Significant", "Not significant")
tf_highlight$Significance <- ifelse(tf_highlight$Significance.Score > 1.30103 | tf_highlight$Significance.Score < -1.30103, "Significant", "Not significant")

# tf$Significance <- ifelse(tf$Site.P.Value < 0.05, "Significant", "Not significant")
tf$Significance <- ifelse(tf$Significance.Score > 1.30103 | tf$Significance.Score < -1.30103, "Significant", "Not significant")
tf$Deficit_Level <- ifelse(tf$Log2.Enrichment < 0, "Deficit", "Enriched")

# Add highlight flag
tf$Highlight <- tf$Transcription.Factor.Name %in% TF

# Create the plot with size mapped to significance score
p <- ggplot(tf, aes(x = Average.Log2.Proportion.Bound, 
               y = Log2.Enrichment)) +
    # Main points with size based on absolute significance score
    geom_point(aes(color = Significance, 
                   alpha = abs(Significance.Score), #Deficit_Level,
                   size = abs(Significance.Score)),  # Use absolute value
               shape = 16) +
    # Highlighted TFs as hollow circles with size based on their significance
    geom_point(data = tf_highlight, #subset(tf, Highlight),
               aes(size = abs(Significance.Score),
               color = Significance)) + #,"#E98415", #"#E63946",
               # shape = 1,  # Hollow circles
               # stroke = 1.2) +  # Thicker border for hollow circles
    # Labels for highlighted TFs
    geom_text_repel(data = subset(tf, Highlight),
                    aes(label = Transcription.Factor.Name),
                    size = 2.8,
                    max.overlaps = 20,
                    box.padding = 0.5,
                    segment.color = "gray50",
                    segment.alpha = 0.5) +
    # Horizontal line at y=0
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    # Color scale
    scale_color_manual(values = c("Not significant" = "#457B9D", 
                                  "Significant" = "#E63946")) +
    # Alpha scale
    # scale_alpha_manual(values = c("Deficit" = 0.9, "Enriched" = 0.5)) +
	scale_alpha_continuous(name = "|Log2 Enrichment|",
                         range = c(0.3, 1),  # Min and max transparency
                         breaks = seq(0, max(abs(tf$Log2.Enrichment)), 
                                      length.out = 5)) +
    # Size scale - adjust range as needed
    scale_size_continuous(name = "Significance Score\n(absolute value)",
                          range = c(0.5, 4),  # Min and max point sizes
                          breaks = seq(0, max(abs(tf$Significance.Score)), 
                                       length.out = 5)) +
    # Labels
    labs(x = "Average.Log2.Proportion.Bound",
         y = "Log2.Enrichment",
         title = "Transcription Factor Enrichment Plot",
         subtitle = "MuHV4 vs Mock at d120") +
    # Theme
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right",
          legend.box = "vertical")
print(p)
# {png("3D_ciider_MA_plot_MuHV4_vs_Mock_D120.png", width = 10, height = 8, units = "in", res = 300);print(p);dev.off()}
{pdf("3D_ciider_MA_plot_MuHV4_vs_Mock_D120.pdf", width = 10, height = 8);print(p);dev.off()}

# - fig S3E -
# # ciider_MA_ plot MuHV-4 versus del73 D120
tf <- read_excel("../Table/3E_ciiider_MuHV4_D120_vs_Del73_D120_MostSigDeficit.xlsx")
# tf <- read.csv("ciiider_MuHV4_D120_vs_Del73_D120_MostSigDeficit.csv")

tf$`Transcription Factor Name` <- gsub("\\(var,2\\)","",tf$`Transcription Factor Name`)
tf$`Transcription Factor Name` <- gsub("\\(var,3\\)","",tf$`Transcription Factor Name`)
TF <- c("CEBPD","CEBPACEBPE","PLAGL2","STAT1::STAT2","CEBPB","Stat2","IRF2","IRF1","HIF1A","IRF7","STAT3","FOSL2::JUNB",
"FOSL1","JUNB","POU3F2","SP1","FOSL2::JUN")

colnames(tf) <- gsub(" ", ".", colnames(tf))
tf_highlight <- tf[tf$Transcription.Factor.Name %in% TF, ]
# tf_highlight$Significance <- ifelse(tf_highlight$Site.P.Value < 0.05, "Significant", "Not significant")
tf_highlight$Significance <- ifelse(tf_highlight$Significance.Score > 1.30103 | tf_highlight$Significance.Score < -1.30103, "Significant", "Not significant")

# tf$Significance <- ifelse(tf$Site.P.Value < 0.05, "Significant", "Not significant")
tf$Significance <- ifelse(tf$Significance.Score > 1.30103 | tf$Significance.Score < -1.30103, "Significant", "Not significant")
tf$Deficit_Level <- ifelse(tf$Log2.Enrichment < 0, "Deficit", "Enriched")

# Add highlight flag
tf$Highlight <- tf$Transcription.Factor.Name %in% TF

# Create the plot with size mapped to significance score
p <- ggplot(tf, aes(x = Average.Log2.Proportion.Bound, 
               y = Log2.Enrichment)) +
    # Main points with size based on absolute significance score
    geom_point(aes(color = Significance, 
                   alpha = abs(Significance.Score), #Deficit_Level,
                   size = abs(Significance.Score)),  # Use absolute value
               shape = 16) +
    # Highlighted TFs as hollow circles with size based on their significance
    geom_point(data = tf_highlight, #subset(tf, Highlight),
               aes(size = abs(Significance.Score),
               color = Significance)) + #,"#E98415", #"#E63946",
               # shape = 1,  # Hollow circles
               # stroke = 1.2) +  # Thicker border for hollow circles
    # Labels for highlighted TFs
    geom_text_repel(data = subset(tf, Highlight),
                    aes(label = Transcription.Factor.Name),
                    size = 2.8,
                    max.overlaps = 20,
                    box.padding = 0.5,
                    segment.color = "gray50",
                    segment.alpha = 0.5) +
    # Horizontal line at y=0
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    # Color scale
    scale_color_manual(values = c("Not significant" = "#457B9D", 
                                  "Significant" = "#E63946")) +
    # Alpha scale
    # scale_alpha_manual(values = c("Deficit" = 0.9, "Enriched" = 0.5)) +
	scale_alpha_continuous(name = "|Log2 Enrichment|",
                         range = c(0.3, 1),  # Min and max transparency
                         breaks = seq(0, max(abs(tf$Log2.Enrichment)), 
                                      length.out = 5)) +
    # Size scale - adjust range as needed
    scale_size_continuous(name = "Significance Score\n(absolute value)",
                          range = c(0.5, 4),  # Min and max point sizes
                          breaks = seq(0, max(abs(tf$Significance.Score)), 
                                       length.out = 5)) +
    # Labels
    labs(x = "Average.Log2.Proportion.Bound",
         y = "Log2.Enrichment",
         title = "Transcription Factor Enrichment Plot",
         subtitle = "MuHV4 vs Del73 at d120") +
    # Theme
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right",
          legend.box = "vertical")
print(p)
# {png("3E_ciider_MA_plot_MuHV4_vs_Del73_D120.png", width = 10, height = 8, units = "in", res = 300);print(p);dev.off()}
{pdf("3E_ciider_MA_plot_MuHV4_vs_Del73_D120.pdf", width = 10, height = 8);print(p);dev.off()}

