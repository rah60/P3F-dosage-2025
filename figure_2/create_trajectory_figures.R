rm(list = ls())
graphics.off()

#plots trajectory analysis figures added in revision, Supplemental Figure 5

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(harmony)
library(fgsea)
library(msigdbr)
library(clustree)
library(Cairo)
library(DoubletFinder)
library(stringr)
library(readr) 
library(viridis)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(qs)
library(ggVennDiagram)
library(slingshot)
library(SingleCellExperiment)
library(igraph)

tol_muted <- c('#332288','#88CCEE','#44AA99','#117733','#999933','#DDCC77','#CC6677','#882255','#AA4499')
source("~/dosage_manuscript/figure_2/figure_2_functions.R")

########## p3f withdrawal

dbt_3_trajectory <- qread("/home/gdstantonlab/lab/sharing/Ambuj/PAX3_FOXO1/2026/PAX3_FOXO1_Rachel_March13/Phase2/slingshot_obj_phase2.qs")
dbt_3 <- qread("/home/gdstantonlab/lab/sharing/Ambuj/PAX3_FOXO1/2026/PAX3_FOXO1_Rachel_March13/Phase2/seurat_phase2_annotated.qs")

metadata <- dbt_3$cell_label
metadata2 <- sapply(str_replace(metadata, "PAX3-FOXO1","PAX3::FOXO1"), "[",1)
names(metadata2) <- names(metadata)
metadata2 <- factor( metadata2, levels= c("Cycling","Progenitor-like","PAX3::FOXO1-specific","Differentiated"))
dbt_3$subclusters_3 <- metadata2

cluster_colors <- c(
  "Progenitor-like" = tol_muted[3],
  "Cycling" = tol_muted[1],
  "PAX3::FOXO1-specific" = tol_muted[8],
  "Differentiated" = tol_muted[7]
)

cell_plot <- plot_slingshot_trajectory_trimmed(
  dbt_3,
  dbt_3_trajectory,
  annot_column = "subclusters_3",
  cluster_colors = cluster_colors,
  point_size = 1.0,
  point_alpha = 0.4
) + guides(colour = guide_legend(override.aes = list(size=4))) + theme(legend.title=element_blank())

png("~/dosage_manuscript/figure_2/p3f_removal_trajectory.png", width = 8, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(cell_plot)
dev.off()

######### Long timepoint

dbt_2_trajectory <- qread("/home/gdstantonlab/lab/sharing/Ambuj/PAX3_FOXO1/2026/PAX3_FOXO1_Rachel_March13/Phase1/slingshot_obj_phase1.qs")
dbt_2 <- qread("/home/gdstantonlab/lab/sharing/Ambuj/PAX3_FOXO1/2026/PAX3_FOXO1_Rachel_March13/Phase1/seurat_phase1_annotated.qs")

metadata <- dbt_2$cell_label
metadata2 <- sapply(str_replace(metadata, "PAX3-FOXO1","PAX3::FOXO1"), "[",1)
metadata2 <- sapply(str_replace(metadata2, "Ground-state","Ground state"), "[",1)
names(metadata2) <- names(metadata)
metadata2 <- factor( metadata2, levels= c("Progenitor-like","Cycling","Ground state","PAX3::FOXO1-specific","Differentiated"))
dbt_2$subclusters_3 <- metadata2

Idents(dbt_2) <- "subclusters_3"

cluster_colors2 <- c(
  "Progenitor-like" = tol_muted[3],
  "Cycling" = tol_muted[1],
  "Ground state" = tol_muted[4],
  "PAX3::FOXO1-specific" = tol_muted[8],
  "Differentiated" = tol_muted[7]
)

cell_plot2 <- plot_slingshot_trajectory_trimmed(
  dbt_2,
  dbt_2_trajectory,
  annot_column = "subclusters_3",
  cluster_colors = cluster_colors2,
  point_size = 1.0,
  point_alpha = 0.4
) + guides(colour = guide_legend(override.aes = list(size=4))) + theme(legend.title=element_blank())

png("~/dosage_manuscript/figure_2/long_timept_trajectory.png", width = 8, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(cell_plot2)
dev.off()
