rm(list = ls())
graphics.off()

#plots all figures in S2 & S3, panels A-D in S4, and panel A in S6

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
library(ggsignif)
library(rstatix)

source("~/dosage_manuscript/figure_2/figure_2_functions.R")

tol_muted <- c('#332288','#88CCEE','#44AA99','#117733','#999933','#DDCC77','#CC6677','#882255','#AA4499')

#distribution of cells across 2 wk clusters  colors: colors from figure 3
col.v1 <- viridis(8)
col.v <- col.v1[1:7] # 0, 75, 150, 250, 500, 1000, common  so 0, 75, 500 is 1, 2, 5
col_3 <- c(col.v[1],col.v[2],col.v[5])
col_4 <- c(tol_muted[2],tol_muted[9])
 
library(enrichR)
websiteLive <- getOption("enrichR.live")
options(enrichR.base.address="https://amp.pharm.mssm.edu/Enrichr/")
dbs_to_use <- c("MSigDB_Hallmark_2020", "ChEA_2022", "GO_Biological_Process_2025")

##################################### short time points

dbt_1 <- qread("/home/gdstantonlab/lab/sharing/Ambuj/PAX3_FOXO1/2026/PAX3_FOXO1_Rachel_March13/Phase_8_24_hrs/seurat_induction_8_24_hrs_annotated.qs")
dbt_1$sample_names <- factor(dbt_1$sample_names, levels= c("0 ng/mL, 0 hrs","500 ng/mL, 8 hrs" ,"500 ng/mL, 24 hrs"))

#Extended Data Figure 2A

#cell cycle 

p2 <- DimPlot(dbt_1, group.by="Phase", pt.size=0.8) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4))) + labs(title="")

png("~/dosage_manuscript/figure_2/short_timept_S2_dimplot_revision.png", width = 9, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()

metadata <- dbt_1$cell_label
metadata2 <- sapply(str_replace(metadata, "Ground-state","Ground state"), "[",1)
names(metadata2) <- names(metadata)
metadata2 <- factor(metadata2, levels= c("Cycling","Ground state"))
dbt_1$subclusters_3 <- metadata2

Idents(dbt_1) <- "subclusters_3"

p2_1 <- DimPlot(dbt_1, cols=c(tol_muted[9],tol_muted[2]), split.by="sample_names", pt.size=0.8) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4))) + labs(title="")

png("~/dosage_manuscript/figure_2/short_timept_reviewer_fig_dimplot.png", width = 13, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2_1)
dev.off()

# Extended Data Figure SA, B

markers <- FindAllMarkers(dbt_1, logfc.threshold = 0.25, min.pct=0.1,  min.diff.pct=0.1)

#saveRDS(markers, "~/dosage_manuscript/figure_2/markers_short_timepoint_S2_revision.rds" )
markers <- readRDS("~/dosage_manuscript/figure_2/markers_short_timepoint_S2_revision.rds")

markers <- markers[which(markers$avg_log2FC > 0),] 

markers_gs <- enrichr(markers[which(markers$cluster  == "Ground state"),7], dbs_to_use)
#plot_gs <- plot_gene_sets2(markers_gs[["MSigDB_Hallmark_2020"]],  plot_title = "Ground state", num_sets=5, fill_col=tol_muted[2]) #NULL

#labeled_gs <- plot_gs + geom_text(aes(label=markers_gs[["MSigDB_Hallmark_2020"]]$Overlap), y=0.75 , x=1, color="white",size=6) + theme_classic(base_size=20)

markers_cyc <- enrichr(markers[which(markers$cluster  == "Cycling"),7], dbs_to_use)
plot_cyc <- plot_gene_sets2(markers_cyc[["MSigDB_Hallmark_2020"]],  plot_title = "Cycling", num_sets=5, fill_col=tol_muted[2])

labeled_cyc <- plot_cyc + geom_text(aes(label=markers_cyc[["MSigDB_Hallmark_2020"]]$Overlap[1:5]), y=c(20,20,45,45,45) , x=c(5,4,3,2,1), color=c("white","white","black","black","black"),size=6) + theme_classic(base_size=20)

p4 <- labeled_cyc

png("~/dosage_manuscript/figure_2/short_timept_gsea_plot_revision.png", width =6, height = 3, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p4)
dev.off()

######### again with ChEA

 plot_cyc <- plot_gene_sets2(markers_cyc[["ChEA_2022"]],  plot_title = "Cycling", num_sets=5, fill_col=tol_muted[2])
 labeled_cyc <- plot_cyc + geom_text(aes(label=markers_cyc[["ChEA_2022"]]$Overlap[1:5]), y=c(15,15,15,15,15) , x=c(5,4,3,2,1), color=c("white"),size=4)

 #plot_gs <- plot_gene_sets2(markers_gs[["ChEA_2022"]],  plot_title = "Ground state", num_sets=5, fill_col=tol_muted[9]) #NULL
 p4 <- labeled_cyc

png("~/dosage_manuscript/figure_2/short_timept_chea_plot_revision.png", width =9, height = 3, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p4)
dev.off()

#Extended Data Figure 2E, top

top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

#marker dotplot
p5 <- DotPlot(dbt_1, features = unique(top10$gene))+
        labs(x="Marker genes",y="Cluster")+
        theme(axis.text.x = element_text(angle=45, hjust=1))

png("~/dosage_manuscript/figure_2/short_timept_markers_plot_revision.png", width = 6, height = 3, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p5)
dev.off()
        

#danielli scores by cluster

#Extended Data Figure 2C, left

dbt_1_scores <- FetchData(dbt_1, vars=c("subclusters_3","sample_names","Progenitor_Score1","Proliferative_Score2","Differentiated_Score3")) 

aov_res <- aov(Progenitor_Score1 ~ subclusters_3, data = dbt_1_scores) %>% tukey_hsd()

aov_res2 <- aov(Proliferative_Score2 ~ subclusters_3, data = dbt_1_scores) %>% tukey_hsd()

aov_res3 <- aov(Differentiated_Score3 ~ subclusters_3, data = dbt_1_scores) %>% tukey_hsd()

my_comparisons <- list(c("Cycling","Ground state"))
sig_levels <- c("****")

p6 <- ggplot(dbt_1_scores,aes(y=Progenitor_Score1, x=subclusters_3, fill=subclusters_3))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 25)+
    scale_fill_manual(values=col_4)+
    labs(y="Progenitor\nScore")+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.title.y=element_text(size=25))+
    scale_y_continuous(limits=c(min(dbt_1_scores$Progenitor_Score1),0.7))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=sig_levels, y_position=c(0.6), show.legend=F, color = "black",textsize=6)

p7 <- ggplot(dbt_1_scores,aes(y=Proliferative_Score2, x=subclusters_3, fill=subclusters_3))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 25)+
    scale_fill_manual(values=col_4)+
    labs(y="Proliferative\nScore")+
    theme(axis.title.x = element_blank(),axis.text.x =element_blank(),axis.title.y=element_text(size=25))+
    scale_y_continuous(limits=c(min(dbt_1_scores$Proliferative_Score2),1.2))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=sig_levels, y_position=c(0.9), show.legend=F, color = "black",,textsize=6)

p8 <- ggplot(dbt_1_scores,aes(y=Differentiated_Score3, x=subclusters_3, fill=subclusters_3))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 25)+
    scale_fill_manual(values=col_4)+
    labs(y="Differentiated\nScore")+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45, vjust=1, hjust=1,size=25),axis.title.y=element_text(size=25))+
    scale_y_continuous(limits=c(min(dbt_1_scores$Differentiated_Score3),0.8))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=sig_levels, y_position=c(0.65), show.legend=F, color = "black",textsize=6)

png("~/dosage_manuscript/figure_2/short_timept_score_boxplots_revision.png", width = 6, height = 10, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p6 / p7 / p8)
dev.off()


#Extended Figure 2D
###########
n_cells <- FetchData(dbt_1, #look at cells per sample in each cluster
                     vars=c("subclusters_3","sample_names")) %>%
          count(subclusters_3, sample_names, name="value") 

n_totals <- FetchData(dbt_1, #look at cells per sample
                     vars=c("sample_names")) %>%
          count(sample_names) 

n_cells$cluster_total <- NA
for(i in 1:nrow(n_cells)){

  row_num <- which(n_totals[,1] %in% n_cells[i,2])
  n_cells[i,4] <- n_totals[row_num,2]

}

n_cells$percent <- n_cells$value/n_cells$cluster_total

p3 <- ggplot(n_cells, aes( x=sample_names, fill=sample_names, y=percent ))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(x="Doxycycline (ng/mL)",y="Fraction of\nsample",fill="Sample")+
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),legend.position = "none")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = 0 ) )+ #don't expand y-axis on bottom, expand y-axis on top by 10%
  facet_rep_wrap(vars(subclusters_3), nrow=1,scales="free" )+ 
  scale_fill_manual(values = col_3)

png("~/dosage_manuscript/figure_2/short_timept_cluster_distribution_revision.png", width = 7, height = 4, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p3)
dev.off()

##################################### long time points

dbt_2 <- qread("/home/gdstantonlab/lab/sharing/Ambuj/PAX3_FOXO1/2026/PAX3_FOXO1_Rachel_March13/Phase1/seurat_phase1_annotated.qs")

col.v2 <- c(tol_muted[3], tol_muted[1],tol_muted[4],tol_muted[8],tol_muted[7])

metadata <- dbt_2$cell_label
metadata2 <- sapply(str_replace(metadata, "PAX3-FOXO1","PAX3::FOXO1"), "[",1)
metadata2 <- sapply(str_replace(metadata2, "Ground-state","Ground state"), "[",1)
names(metadata2) <- names(metadata)
metadata2 <- factor( metadata2, levels= c("Progenitor-like","Cycling","Ground state","PAX3::FOXO1-specific","Differentiated"))
dbt_2$subclusters_3 <- metadata2

Idents(dbt_2) <- "subclusters_3"

#Extended Data Figure 2B

#cell cycle 

p2 <- DimPlot(dbt_2, group.by="Phase", pt.size=0.8) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4))) + labs(title="")

png("~/dosage_manuscript/figure_2/long_timept_S2_dimplot_revision.png", width = 9, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()

dbt_2$sample_names <- factor(dbt_2$sample_names, levels=c("0 ng/mL, 2 wks","75 ng/mL, 2 wks","500 ng/mL, 2 wks"))
cols_dbt2 <-  c(tol_muted[3], tol_muted[1],tol_muted[4],tol_muted[8],tol_muted[7])
p2_2 <- DimPlot(dbt_2,cols=cols_dbt2, split.by="sample_names", pt.size=0.8) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4))) + labs(title="")

png("~/dosage_manuscript/figure_2/long_timept_reviewer_figure_dimplot.png", width = 13, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2_2)
dev.off()


#Extended Data Figure 3C, D

markers <- FindAllMarkers(dbt_2, logfc.threshold = 0.25, min.pct=0.1,  min.diff.pct=0.1)

#saveRDS(markers, "~/dosage_manuscript/figure_2/markers_long_timepoint_S2_revised.rds" )

markers <- readRDS("~/dosage_manuscript/figure_2/markers_long_timepoint_S2_revised.rds" )

markers <- markers[which(markers$avg_log2FC > 0),] 

markers_cyc <- enrichr(markers[which(markers$cluster  == "Cycling"),7], dbs_to_use)
plot_cyc <- plot_gene_sets2(markers_cyc[["MSigDB_Hallmark_2020"]],  plot_title = "Cycling", num_sets=5, fill_col=tol_muted[1])

labeled_cyc <- plot_cyc + geom_text(aes(label=markers_cyc[["MSigDB_Hallmark_2020"]]$Overlap[1:5]), y=c(25,25,50,50,50) , x=c(5,4,3,2,1), color=c("white","white","black","black","black"),size=6) + theme_classic(base_size=20)


markers_pro <- enrichr(markers[which(markers$cluster  == "Progenitor-like"),7], dbs_to_use)
markers_pro[["MSigDB_Hallmark_2020"]][which(markers_pro[["MSigDB_Hallmark_2020"]] == "Epithelial Mesenchymal Transition"),1] <- "EMT"
plot_pro <- plot_gene_sets2(markers_pro[["MSigDB_Hallmark_2020"]],  plot_title = "Progenitor-like", num_sets=5, fill_col=tol_muted[3])

labeled_pro <- plot_pro + geom_text(aes(label=markers_pro[["MSigDB_Hallmark_2020"]]$Overlap[1:5]), y=c(8,8,8,15,15) , x=c(5,4,3,2,1), color=c("white","white","white","black","black"),size=6) + theme_classic(base_size=20)


markers_gs <- enrichr(markers[which(markers$cluster  == "Ground state"),7], dbs_to_use)
plot_gs <- plot_gene_sets2(markers_gs[["MSigDB_Hallmark_2020"]],  plot_title = "Ground state", num_sets=5, fill_col=tol_muted[4])

labeled_gs <- plot_gs + geom_text(aes(label=markers_gs[["MSigDB_Hallmark_2020"]]$Overlap[1:3]), y=c(5,5,5) , x=c(3,2,1), color=c("white","black","black"),size=6) + theme_classic(base_size=20)


markers_diff <- enrichr(markers[which(markers$cluster  == "Differentiated"),7], dbs_to_use)
plot_diff <- plot_gene_sets2(markers_diff[["MSigDB_Hallmark_2020"]],  plot_title = "Differentiated", num_sets=5, fill_col=tol_muted[7])

labeled_diff <- plot_diff + geom_text(aes(label=markers_diff[["MSigDB_Hallmark_2020"]]$Overlap[1:5]), y=c(12,12,12,12,12) , x=c(5,4,3,2,1), color=c("white","black","black","black","black"),size=6) + theme_classic(base_size=20)


markers_p3f <- enrichr(markers[which(markers$cluster  == "PAX3::FOXO1-specific"),7], dbs_to_use)
plot_p3f <- plot_gene_sets2(markers_p3f[["MSigDB_Hallmark_2020"]],  plot_title = "P3F-specific", num_sets=5, fill_col=tol_muted[8])

labeled_p3f <- plot_p3f + geom_text(aes(label=markers_p3f[["MSigDB_Hallmark_2020"]]$Overlap[1:2]), y=c(1,1) , x=c(2,1), color="white",size=6) + theme_classic(base_size=20)

p4 <- plot_grid(labeled_cyc ,labeled_pro, labeled_gs, labeled_p3f, labeled_diff, ncol=1)

png("~/dosage_manuscript/figure_2/long_timept_gsea_plot_revised.png", width =6, height = 15, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p4)
dev.off()

#### again with chea

#Extended Data Figure 2G

plot_cyc <- plot_gene_sets2(markers_cyc[["ChEA_2022"]],  plot_title = "Cycling", num_sets=5, fill_col=tol_muted[1])
labeled_cyc <- plot_cyc + geom_text(aes(label=markers_cyc[["ChEA_2022"]]$Overlap[1:5]), y=c(30,30,30,30,72) , x=c(5,4,3,2,1), color=c("white","white","white","white","black"),size=6) + theme_classic(base_size=20)

plot_pro <- plot_gene_sets2(markers_pro[["ChEA_2022"]],  plot_title = "Progenitor-like", num_sets=5, fill_col=tol_muted[3])
labeled_pro <- plot_pro + geom_text(aes(label=markers_pro[["ChEA_2022"]]$Overlap[1:5]), y=c(7,7,7,7,7) , x=c(5,4,3,2,1), color="white",size=6) + theme_classic(base_size=20)

plot_gs <- plot_gene_sets2(markers_gs[["ChEA_2022"]],  plot_title = "Ground state", num_sets=5, fill_col=tol_muted[4])
labeled_gs <- plot_gs + geom_text(aes(label=markers_gs[["ChEA_2022"]]$Overlap[1:5]), y=c(5,5,5,5,5) , x=c(5,4,3,2,1), color="white",size=6) + theme_classic(base_size=20)

plot_diff <- plot_gene_sets2(markers_diff[["ChEA_2022"]],  plot_title = "Differentiated", num_sets=5, fill_col=tol_muted[7])
labeled_diff <- plot_diff + geom_text(aes(label=markers_diff[["ChEA_2022"]]$Overlap[1:5]), y=c(4,4,4,4,4) , x=c(5,4,3,2,1), color="white",size=6) + theme_classic(base_size=20)

plot_p3f <- plot_gene_sets2(markers_p3f[["ChEA_2022"]],  plot_title = "PAX3::FOXO1-specific", num_sets=5, fill_col=tol_muted[8])
labeled_p3f <- plot_p3f + geom_text(aes(label=markers_p3f[["ChEA_2022"]]$Overlap[1:5]), y=c(9,9,9,9,9) , x=c(5,4,3,2,1), color=c("white","black","black","black","black"),size=6) + theme_classic(base_size=20)

p6 <- plot_grid(labeled_cyc, labeled_pro, labeled_diff, labeled_p3f, labeled_gs, ncol=1) 

png("~/dosage_manuscript/figure_2/long_timept_gsea_plot_chea_revised.png", width = 12, height = 15, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p6)
dev.off()


#Extended Data Figure 2E, bottom

top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top10_2 <-  markers %>% group_by(cluster) %>% top_n(n = 20, wt = pct.1)
top10_2[which(top10_2$gene %in% top10$gene),]

cycling <- c("TOP2A", "PCLAF", "CKS1B","UBE2S", "HMGB2")

progenitor <- c("THBS1", "TIMP3", "TAGLN","COL1A1", "TPM1")

ground_state <- c("ATP5F1D","DDT","MZT2A","SNRPD3","DRAP1")

differentiated <- c("MYLPF", "TNNC1","ACTC1","TNNI1", "MYL4")

p3f_specific <- c("MYOD1","FGF8","XRN1","GAS1","FGFR4")

other_markers <- unique( c(differentiated,p3f_specific,ground_state,cycling, progenitor))

#marker dotplot
p5 <- DotPlot(dbt_2, features = other_markers, scale=T)+
        labs(x="Marker genes",y="Cluster")+
        theme(axis.text.x = element_text(angle=45, hjust=1))

png("~/dosage_manuscript/figure_2/long_timept_markers_plot_2_revisions.png", width = 12, height = 3, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p5)
dev.off()


#Extended Data Figure 2C, right

dbt_2_scores <- FetchData(dbt_2, vars=c("subclusters_3","sample_names","Progenitor_Score1","Proliferative_Score2","Differentiated_Score3")) 

aov_res <- aov(Progenitor_Score1 ~ subclusters_3, data = dbt_2_scores) %>% tukey_hsd()
my_comparisons <- list(c("Cycling","Ground state"), c("Cycling","PAX3::FOXO1-specific"),c("Cycling","Differentiated"),
    c("Progenitor-like","Cycling"),c("Progenitor-like","Ground state"),c("Progenitor-like","PAX3::FOXO1-specific"),c("Progenitor-like","Differentiated"))
sig_levels <- c("****","****","****","****","****","****","****")

aov_res2 <- aov(Proliferative_Score2 ~ subclusters_3, data = dbt_2_scores) %>% tukey_hsd()
my_comparisons2 <- list(c("Progenitor-like","Cycling"),c("Progenitor-like","Ground state"),c("Progenitor-like","PAX3::FOXO1-specific"),
c("Cycling","Ground state"), c("Cycling","PAX3::FOXO1-specific"),c("Cycling","Differentiated"),
    c("Ground state","PAX3::FOXO1-specific"),c("Ground state","Differentiated"),c("PAX3::FOXO1-specific","Differentiated"))
sig_levels2 <- c("****","****","****","****","****","****","****","****","****")

aov_res3 <- aov(Differentiated_Score3 ~ subclusters_3, data = dbt_2_scores) %>% tukey_hsd()
my_comparisons3 <- list(c("Progenitor-like","Cycling"),c("Progenitor-like","Ground state"),c("Progenitor-like","PAX3::FOXO1-specific"),c("Progenitor-like","Differentiated"),
c("Cycling","Ground state"), c("Cycling","PAX3::FOXO1-specific"),c("Cycling","Differentiated"),
    c("Ground state","PAX3::FOXO1-specific"),c("Ground state","Differentiated"),c("PAX3::FOXO1-specific","Differentiated"))
sig_levels3 <- c("****","****","****","****","****","****","****","****","****","****")

p6 <- ggplot(dbt_2_scores,aes(y=Progenitor_Score1, x=subclusters_3, fill=subclusters_3))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 25)+
    scale_fill_manual(values=col.v2)+
    labs(y="Progenitor\nScore")+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.title.y=element_text(size=25))+
    scale_y_continuous(limits=c(min(dbt_2_scores$Progenitor_Score1),1.3))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=sig_levels, y_position=c(0.6,0.7,0.8,0.9,1,1.1,1.2), show.legend=F, color = "black",textsize = 3)


p7 <- ggplot(dbt_2_scores,aes(y=Proliferative_Score2, x=subclusters_3, fill=subclusters_3))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 25)+
    scale_fill_manual(values=col.v2)+
    labs(y="Proliferative\nScore")+
    theme(axis.title.x = element_blank(),axis.text.x =element_blank(),axis.title.y=element_text(size=25))+
    scale_y_continuous(limits=c(min(dbt_2_scores$Proliferative_Score2),1.8))+
    geom_signif(comparisons = my_comparisons2, map_signif_level = TRUE, annotation=sig_levels2, y_position=c(0.95, 1.1,1.25,1.4,1.55,1.7,0.65,0.8,1.1), show.legend=F, color = "black",textsize = 3)

p8 <- ggplot(dbt_2_scores,aes(y=Differentiated_Score3, x=subclusters_3, fill=subclusters_3))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 25)+
    scale_fill_manual(values=col.v2)+
    labs(y="Differentiated\nScore")+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45, vjust=1, hjust=1,size=25),axis.title.y=element_text(size=25))+
    scale_y_continuous(limits=c(min(dbt_2_scores$Differentiated_Score3),1.6))+
    geom_signif(comparisons = my_comparisons3, map_signif_level = TRUE, annotation=sig_levels3, y_position=c(0.8, 0.9, 1,1.1, 0.6, 0.7,1.2, 0.5, 1.3, 1.45), show.legend=F, color = "black",textsize = 3)


png("~/dosage_manuscript/figure_2/long_timept_score_boxplots_revision.png", width = 8, height = 12, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p6 / p7 / p8)
dev.off()

########## checking neural signature, Figure S4

#import neural signature markers
n_markers <- read_csv("~/neural_markers.csv")

#create neural signature module score, plot across clusters, S4A

dbt_2 <- AddModuleScore(object=dbt_2, features=list(n_markers$gene), name=c("Neuronal_Score"), assay="RNA", search=T)

dbt_2_scores2 <- FetchData(dbt_2, vars=c("subclusters_3","sample_names","Neuronal_Score1")) 

aov_res4 <- aov(Neuronal_Score1 ~ subclusters_3, data = dbt_2_scores2) %>% tukey_hsd() #all signifcantly different
my_comparisons4 <- list(c("Progenitor-like","Cycling"),c("Progenitor-like","Ground state"),c("Progenitor-like","PAX3::FOXO1-specific"),c("Progenitor-like","Differentiated"),
c("Cycling","Ground state"), c("Cycling","PAX3::FOXO1-specific"),c("Cycling","Differentiated"),
    c("Ground state","PAX3::FOXO1-specific"),c("Ground state","Differentiated"),c("PAX3::FOXO1-specific","Differentiated"))
sig_levels4 <- c("****","****","****","****","****","****","****","****","****","****")

p9 <- ggplot(dbt_2_scores2,aes(y=Neuronal_Score1, x=subclusters_3, fill=subclusters_3))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 30)+
    scale_fill_manual(values=col.v2)+
    labs(y="Neuronal Score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45, vjust=1, hjust=1),axis.title.y=element_text(size=30))+
    scale_y_continuous(limits=c(min(dbt_2_scores2$Neuronal_Score1),0.8))+
    geom_signif(comparisons = my_comparisons4, map_signif_level = TRUE, annotation=sig_levels4, y_position=c(0.4, 0.45, 0.5,0.55, 0.35, 0.6,0.65, 0.3,0.7, 0.4), show.legend=F, color = "black",textsize = 3)

png("~/dosage_manuscript/figure_2/neuronal_score_vln_plot.png", width = 6, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p9)
dev.off()

#plot across samples, S4B
aov_res5 <- aov(Neuronal_Score1 ~ sample_names, data = dbt_2_scores2) %>% tukey_hsd() #all signifcantly different
my_comparisons5 <- list(c("0 ng/mL, 2 wks","75 ng/mL, 2 wks"),c("500 ng/mL, 2 wks","75 ng/mL, 2 wks"),c("0 ng/mL, 2 wks","500 ng/mL, 2 wks"))
sig_levels5 <- c("****","****","****")
dbt_2_scores2$sample_names <- factor(dbt_2_scores2$sample_names, levels=c("0 ng/mL, 2 wks","75 ng/mL, 2 wks","500 ng/mL, 2 wks"))

p13 <- ggplot(dbt_2_scores2,aes(y=Neuronal_Score1, x=sample_names, fill=sample_names))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 30)+
    scale_fill_manual(values=col_3)+
    labs(y="Neuronal Score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45, vjust=1, hjust=1),axis.title.y=element_text(size=30))+
    scale_y_continuous(limits=c(min(dbt_2_scores$Differentiated_Score3),0.55))+
    geom_signif(comparisons = my_comparisons5, map_signif_level = TRUE, annotation=sig_levels5, y_position=c(0.4, 0.45, 0.5), show.legend=F, color = "black",textsize = 3)


png("~/dosage_manuscript/figure_2/neuronal_score_vln_plot_bysample.png", width = 6, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p13)
dev.off()

### compare neuronal score to differentiated score for the differentiated cluster, S4D

danielli_markers <- read_csv("~/Dbt_scRNA/danielli_2023_markers.csv")
differentiated_markers <- as.vector(danielli_markers$Differentiated_markers) 
differentiated_markers  <- differentiated_markers[!is.na(differentiated_markers)]

length(which(n_markers$gene %in% differentiated_markers))/length(n_markers$gene)
#### 16% of neuronal markers in differentiated markers

length(which(differentiated_markers %in% n_markers$gene))/length(differentiated_markers)
#### 21% of the differentiated markers in neuronal markers

scores.df <- data.frame(matrix(ncol=2, nrow=503))
colnames(scores.df) <- c("gene","cluster") 

#set neuronal & differentiated cluster markers for gsea
cluster_names <- c(rep("Neuronal_cluster",times=(nrow(n_markers))), rep("Differentiated_cluster",times=(length(differentiated_markers))))
gene_names <- c(n_markers$gene,differentiated_markers )

scores.df[,1] <- gene_names
scores.df[,2] <- cluster_names

gene_sets <- split(x=scores.df$gene, f=scores.df$cluster)

cluster_names <- unique(markers$cluster)
perm_num <- 2000
fgsea.l <- vector("list") 

for(k in 1:length(unique(markers$cluster))){ #go through clusters & run fgsea on each, 
      print(k)
      clust_num <- k - 1

      markers.v <- markers[which(markers$cluster == cluster_names[k]),2] #avglogFC

      names(markers.v) <- rownames(markers[which(markers$cluster == cluster_names[k]),])
      
      if(min(markers.v) > 0){ #determine if there are only positive logFC markers
      fgsea.l[[k]] <- fgsea(pathways = gene_sets, stats=markers.v , maxSize=500, minSize=5, nPerm=perm_num, scoreType="pos" ) %>%
                        arrange(-NES) 
                        print("pos") 
                        next }
      if(max(markers.v) < 0){ #negative markers only
      fgsea.l[[k]] <- fgsea(pathways = gene_sets, stats=markers.v , maxSize=500, minSize=5, nPerm=perm_num , scoreType="neg") %>%
                        arrange(-NES)
                        print("neg")
                        next }
      fgsea.l[[k]] <- fgsea(pathways = gene_sets, stats=markers.v , maxSize=500, minSize=5, nPerm=perm_num ) %>% #or both & run correct FGSEA settings
                        arrange(-NES)
                        print("std")
      }

##only differentiated cluster scores as significant for differentiated & neuronal markers
#plot S4D

p10 <- ggplot(fgsea.l[[5]]) +
        geom_col(aes(x=factor(pathway, levels=c("Neuronal_cluster","Differentiated_cluster"),labels=c("Neuronal\nmarkers","Differentiated\nmarkers")), y=-log10(padj), fill="blank"))+
        coord_flip()+
        theme_classic(base_size=28)+
        scale_fill_manual(values=tol_muted[7])+
        scale_y_continuous(expand = c(0,0))+
        labs(x= "Gene sets" , y="-Log10 adjusted p-value", title="Differentiated")+
        guides(fill="none")

p11 <- p10 + geom_text(aes(label=c("75/215","31/288")), y=c(3,3) , x=c(2,1), color=c("white","black"),size=6) + theme_classic(base_size=28)

png("~/dosage_manuscript/figure_2/neuronal_score_gsea.png", width = 7, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p11)
dev.off()

###### venn diagram for all markers in the neuronal & differentiated clusters, S4C
p12 <- ggVennDiagram(gene_sets, label= "count", category.names = c("",""),label_size=12)+ 
    scale_fill_distiller(palette = "PiYG")+
    labs("Marker genes")+ 
    scale_x_continuous(expand = expansion(mult = .2))+
    geom_text(aes(label=c("Differentiated","Neuronal")), y=c(-3,7) , x=c(0,0), color=c("black","black"),size=7)+ 
    theme_void(16)

png("~/dosage_manuscript/figure_2/neuronal_score_venndiagram.png", width = 6, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p12)
dev.off()


######## markers selected for single cell validation, S6A

#P3F-specific 
p3f_2 <- c("GAS1")
p3f_1 <- c("CD9", "CXCR4")

progenitor_2 <- c("THBS1", "TAGLN")
progenitor_1 <- c("LRP1","ALCAM")

cycling_2 <- c("TOP2A","MKI67")
cycling_1 <- c("HMMR")

diff_2 <- c("TNNT2","MYH3")
diff_1 <- c("ERBB3","BCAM")

extracellular_scRNA_flow <- c(diff_1, p3f_1, cycling_1,progenitor_1)
intracellular_scRNA_flow <- c(diff_2, p3f_2, cycling_2,progenitor_2)

#marker dotplot
p13 <- DotPlot(dbt_2, features = extracellular_scRNA_flow, scale=T)+
        labs(x="Marker genes",y="Cluster")+
        theme(axis.text.x = element_text(angle=45, hjust=1)) +labs(title="Cell surface markers")

p14 <- DotPlot(dbt_2, features = intracellular_scRNA_flow, scale=T)+
        labs(x="Marker genes",y="Cluster")+
        theme(axis.text.x = element_text(angle=45, hjust=1)) +labs(title="Intracellular markers",y="")

png("~/dosage_manuscript/figure_2/scRNA_flow_markers_revisions.png", width = 12, height = 4, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p13 + p14)
dev.off()
