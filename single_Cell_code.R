library(devtools)
library(Seurat)
library(tidyverse)
library(patchwork)
library(readxl)
library(hdf5r)
library(presto)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(msigdbr)
library(ggpubr)
library(patchwork)
library(SCP)
library(openxlsx)
source("D:\\sn-rna seq\\fp_vln_creator.R")
merge_scRNA<-readRDS("D:\\Nature2022\\human_ adipocytes.rds")
merge_scRNA$condition <- dplyr::case_when(merge_scRNA$bmi < 25 ~ "Lean",
                                          merge_scRNA$bmi >= 25 & merge_scRNA$bmi < 30 ~ "Overweight",
                                          merge_scRNA$bmi >= 30 ~ "Obese", TRUE ~ NA_character_)
merge_scRNA_obese <- subset(merge_scRNA, subset = condition == "Obese")
merge_scRNA_Lean <- subset(merge_scRNA, subset = condition == "Lean")
merge_scRNA$condition <- factor(merge_scRNA$condition, levels = c("Lean", "Obese"))
conditions <- unique(merge_scRNA$condition)
SCP::FeatureStatPlot(srt = merge_scRNA, group.by = "condition", bg.by = "condition",
                     stat.by = c("TSC22D3"),add_box = TRUE,comparisons = list(c("Lean",
                                                                                "Obese")),palcolor = c("#47A535", "#EF851B"),bg_palcolor = c("#47A535",
                                                                                                                                             "#EF851B"),alpha = 0.7)
中的改变
create_fp(merge_scRNA_obese,feature_vector = "NR3C1",reduction = "umap",cutoff =
            TRUE,pt.size = 3,ncolumn = 10) + theme(aspect.ratio = 1)

custom_genes_up <- read.csv("D:\\Nature 2022\\common_up_genes.csv",
                            stringsAsFactors = FALSE)
head(custom_genes_up)
gc_pathway_genes_up <- custom_genes_up$up_genes 
gc_pathway_genes_up <- unique(na.omit(gc_pathway_genes_up))
genes_present_up<-intersect(gc_pathway_genes_up,rownames(merge_scRNA))
cat(paste0("通路共有", length(gc_pathway_genes_up), " 个基因，在数据中找到",
           length(genes_present_up), " 个\n"))
merge_scRNA <- AddModuleScore(object = merge_scRNA,features =
                                list(genes_present_up),name = "GC_AddModuleScore_up", assay = "RNA",slot =
                                "data",nbin = 24,ctrl = 100 )

colnames(merge_scRNA@meta.data)
merge_scRNA$AddModule_score_up <- merge_scRNA$GC_AddModule Score_up1
addmodule_data_up<-
  FetchData(merge_scRNA,vars=c("AddModule_score_up","condition")) %>%
  filter(condition %in% c("Lean", "Obese"))
addmodule_stats_up <- addmodule_data_up %>% group_by(condition) %>%
  summarise(mean_score = mean(AddModule_score_up, na.rm = TRUE),sem_score =
              sd(AddModule_score_up, na.rm = TRUE) / sqrt(n()),n = n()) %>%mutate(method =
                                                                                    "AddModuleScore_up")

cells_of_interest <- c(GC_activated_cells, GC_nonactivated_cells)
GC_subset <- subset(merge_scRNA, cells = cells_of_interest)
GC_subset$GC_group <- droplevels(GC_subset$GC_group)
table(GC_subset$GC_group)
DefaultAssay(GC_subset) <- "RNA"
GC_markers<-FindMarkers(GC_subset,ident.1 = "GC_activated", ident.2 =
                          "GC_nonactivated", group.by = "GC_group",assay = "RNA",slot =
                          "data",logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox")
GC_markers$gene <- rownames(GC_markers)
write.csv(GC_markers,"D:\\GC_activated_vs_nonactivated_DEGs.csv", row.names =
            FALSE)