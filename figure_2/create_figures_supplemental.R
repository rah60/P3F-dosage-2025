rm(list = ls())
graphics.off()

#plots all supplemental figures

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

source("~/dosage_manuscript/figure_2/figure_2_functions.R")

tol_muted <- c('#332288','#88CCEE','#44AA99','#117733','#999933','#DDCC77','#CC6677','#882255','#AA4499')

#distribution of cells across 2 wk clusters  colors: colors from figure 3
col.v1 <- viridis(8)
col.v <- col.v1[1:7] # 0, 75, 150, 250, 500, 1000, common  so 0, 75, 500 is 1, 2, 5
col_3 <- c(col.v[1],col.v[2],col.v[5])

library(enrichR)
websiteLive <- getOption("enrichR.live")
options(enrichR.base.address="https://amp.pharm.mssm.edu/Enrichr/")
dbs_to_use <- c("MSigDB_Hallmark_2020", "ChEA_2022", "Cancer_Cell_Line_Encyclopedia")

##################################### short time points

dbt_1 <- readRDS("~/dosage_manuscript/rds/short_timepts_batchcorrection_rna_namedmergedclusters_2.rds")

#Extended Data Figure 2A

#cell cycle 
p1 <- DimPlot(dbt_1, cols=c(tol_muted[2],tol_muted[9]), pt.size=0.8 ) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4)))
p2 <- DimPlot(dbt_1, group.by="Phase", pt.size=0.8) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4))) + labs(title="")

p3 <- p1 + p2

png("~/dosage_manuscript/figure_2/short_timept_S2_dimplot.png", width = 16, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p3)
dev.off()


# Extended Data Figure 2E

markers <- FindAllMarkers(dbt_1, logfc.threshold = 0.25, min.pct=0.1,  min.diff.pct=0.1)

#saveRDS(markers, "~/dosage_manuscript/figure_2/markers_short_timepoint_S2.rds" )
markers <- readRDS("~/dosage_manuscript/figure_2/markers_short_timepoint_S2.rds")

markers <- markers[which(markers$avg_log2FC > 0),] 

markers_gs <- enrichr(markers[which(markers$cluster  == "ground state"),7], dbs_to_use)
plot_gs <- plot_gene_sets2(markers_gs[["MSigDB_Hallmark_2020"]],  plot_title = "Ground state", num_sets=5, fill_col=tol_muted[2])

labeled_gs <- plot_gs + geom_text(aes(label=markers_gs[["MSigDB_Hallmark_2020"]]$Overlap), y=0.75 , x=1, color="white",size=6) + theme_classic(base_size=20)

markers_cyc <- enrichr(markers[which(markers$cluster  == "cycling"),7], dbs_to_use)
plot_cyc <- plot_gene_sets2(markers_cyc[["MSigDB_Hallmark_2020"]],  plot_title = "Cycling", num_sets=5, fill_col=tol_muted[9])

labeled_cyc <- plot_cyc + geom_text(aes(label=markers_cyc[["MSigDB_Hallmark_2020"]]$Overlap[1:5]), y=c(20,20,45,45,45) , x=c(5,4,3,2,1), color=c("white","white","black","black","black"),size=6) + theme_classic(base_size=20)

p4 <- labeled_gs / labeled_cyc

png("~/dosage_manuscript/figure_2/short_timept_gsea_plot.png", width = 6, height = 10, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p4)
dev.off()

######### again with ChEA

# plot_cyc <- plot_gene_sets2(markers_cyc[["ChEA_2022"]],  plot_title = "Cycling", num_sets=5, fill_col=tol_muted[2])
# labeled_cyc <- plot_cyc + geom_text(aes(label=markers_cyc[["ChEA_2022"]]$Overlap[1:5]), y=c(15,15,15,15,15) , x=c(5,4,3,2,1), color=c("white"),size=4)

# plot_gs <- plot_gene_sets2(markers_gs[["ChEA_2022"]],  plot_title = "Ground state", num_sets=5, fill_col=tol_muted[9]) #NULL


#Extended Data Figure 2H, left

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#marker dotplot
p5 <- DotPlot(dbt_1, features = unique(top10$gene))+
        labs(x="Marker genes",y="Cluster")+
        theme(axis.text.x = element_text(angle=45, hjust=1))

png("~/dosage_manuscript/figure_2/short_timept_markers_plot_2.png", width = 6, height = 3, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p5)
dev.off()
        

#danielli scores by cluster

#Extended Data Figure 2C, left

dbt_1_scores <- FetchData(dbt_1, vars=c("clusters","sample_names","Progenitor_Score1","Proliferative_Score2","Differentiated_Score3")) 

p6 <- ggplot(dbt_1_scores,aes(y=Progenitor_Score1, x=clusters, fill=clusters))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 20)+
    scale_fill_manual(values=col.v)+
    labs(y="Progenitor Score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.title.y=element_text(size=16))

p7 <- ggplot(dbt_1_scores,aes(y=Proliferative_Score2, x=clusters, fill=clusters))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 20)+
    scale_fill_manual(values=col.v)+
    labs(y="Proliferative Score")+
    theme(axis.title.x = element_blank(),axis.text.x =element_blank(),axis.title.y=element_text(size=16))

p8 <- ggplot(dbt_1_scores,aes(y=Differentiated_Score3, x=clusters, fill=clusters))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 20)+
    scale_fill_manual(values=col.v)+
    labs(y="Differentiated Score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45, vjust=1, hjust=1),axis.title.y=element_text(size=16))

png("~/dosage_manuscript/figure_2/short_timept_score_boxplots.png", width = 4, height = 8, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p6 / p7 / p8)
dev.off()


#Extended Figure 2D, left 
# ##########
n_cells <- FetchData(dbt_1, #look at cells per sample in each cluster
                     vars=c("clusters","sample_names")) %>%
          count(clusters, sample_names) %>%
          tidyr::spread(clusters, n) %>%
          pivot_longer(cols=c("ground state", "cycling"), names_to="cluster")

p3 <- ggplot(n_cells, aes( x=sample_names, fill=sample_names, y=value ))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(x="Sample",y="Cell Count",fill="Sample")+
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),legend.position = "none")+
  facet_wrap(vars(cluster), nrow=3, )+
  scale_fill_manual(values = col_3)

png("~/dosage_manuscript/figure_2/short_timept_cell_cycle_distribution.png", width = 4, height = 7, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p3)
dev.off()

##################################### long time points

dbt_2 <- readRDS("~/dosage_manuscript/rds/dosages_2wks_batchcorrection_rna_subclustered_named1.rds")

col.v2 <- c(tol_muted[1],tol_muted[3:5],tol_muted[7:8])

metadata <- dbt_2$subclusters_2
metadata2 <- sapply(str_replace(metadata, "PAX3-FOXO1","PAX3::FOXO1"), "[",1)
names(metadata2) <- names(metadata)
metadata2 <- factor(metadata2, levels = c("Cycling","Progenitor-like", "Ground state","AP-1 expressing","Differentiated","PAX3::FOXO1-specific"))
dbt_2$subclusters_2 <- metadata2

Idents(dbt_2) <- "subclusters_2"

#Extended Data Figure 2B

#cell cycle 
p1 <- DimPlot(dbt_2, cols=col.v2, pt.size=0.8 ) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4)))
p2 <- DimPlot(dbt_2, group.by="Phase", pt.size=0.8) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4))) + labs(title="")

p3 <- p1 + p2

png("~/dosage_manuscript/figure_2/long_timept_S2_dimplot.png", width = 16, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p3)
dev.off()

#Extended Data figure 2D, left

metadata <- dbt_2$subclusters_2
metadata2 <- sapply(str_replace(metadata, "expressing","exp."), "[",1)
metadata2 <- sapply(str_replace(metadata2, "PAX3-FOXO1","P3F"), "[",1)
metadata2 <- names(metadata2)
dbt_2$abbrev_cluster_names <- metadata2

n_cells <- FetchData(dbt_2, #look at cells per sample in each cluster
                     vars=c("abbrev_cluster_names","sample_names")) %>%
          count(abbrev_cluster_names, sample_names) %>%
          tidyr::spread(abbrev_cluster_names, n) %>%
          pivot_longer(cols=c("Cycling","Progenitor-like", "Ground state","AP-1 exp.", "Differentiated", "P3F-specific"), names_to="cluster")

p3 <- ggplot(n_cells, aes( x=sample_names, fill=sample_names, y=value ))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(x="Sample",y="Cell Count",fill="Sample")+
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),legend.position = "none")+
  facet_wrap(vars(cluster), nrow=3, )+
  scale_fill_manual(values = col_3)

png("~/dosage_manuscript/figure_2/long_timept_cluster_distribution_newscale.png", width = 6, height = 8, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p3)
dev.off()

#Extended Data Figure 2F

markers <- FindAllMarkers(dbt_2, logfc.threshold = 0.25, min.pct=0.1,  min.diff.pct=0.1)

#saveRDS(markers, "~/dosage_manuscript/figure_2/markers_long_timepoint_S2.rds" )

markers <- readRDS("~/dosage_manuscript/figure_2/markers_long_timepoint_S2.rds" )

markers <- markers[which(markers$avg_log2FC > 0),] 

markers_cyc <- enrichr(markers[which(markers$cluster  == "Cycling"),7], dbs_to_use)
plot_cyc <- plot_gene_sets2(markers_cyc[["MSigDB_Hallmark_2020"]],  plot_title = "Cycling", num_sets=5, fill_col=col.v2[1])

labeled_cyc <- plot_cyc + geom_text(aes(label=markers_cyc[["MSigDB_Hallmark_2020"]]$Overlap[1:5]), y=c(25,25,50,50,50) , x=c(5,4,3,2,1), color=c("white","white","black","black","black"),size=6) + theme_classic(base_size=20)


markers_pro <- enrichr(markers[which(markers$cluster  == "Progenitor-like"),7], dbs_to_use)
plot_pro <- plot_gene_sets2(markers_pro[["MSigDB_Hallmark_2020"]],  plot_title = "Progenitor-like", num_sets=5, fill_col=col.v2[2])

labeled_pro <- plot_pro + geom_text(aes(label=markers_pro[["MSigDB_Hallmark_2020"]]$Overlap[1:5]), y=c(5,5,10,10,10) , x=c(5,4,3,2,1), color=c("white","white","black","black","black"),size=6) + theme_classic(base_size=20)


markers_gs <- enrichr(markers[which(markers$cluster  == "Ground state"),7], dbs_to_use)
plot_gs <- plot_gene_sets2(markers_gs[["MSigDB_Hallmark_2020"]],  plot_title = "Ground state", num_sets=5, fill_col=col.v2[3])

labeled_gs <- plot_gs + geom_text(aes(label=markers_gs[["MSigDB_Hallmark_2020"]]$Overlap[1:5]), y=c(2.5,2.5,4,4,4) , x=c(5,4,3,2,1), color=c("white","white","black","black","black"),size=6) + theme_classic(base_size=20)


markers_ap1 <- enrichr(markers[which(markers$cluster  == "AP-1 expressing"),7], dbs_to_use)
plot_ap1 <- plot_gene_sets2(markers_ap1[["MSigDB_Hallmark_2020"]],  plot_title = "AP-1 expressing", num_sets=5, fill_col=col.v2[4])

labeled_ap1 <- plot_ap1 + geom_text(aes(label=markers_ap1[["MSigDB_Hallmark_2020"]]$Overlap[1:5]), y=c(0.5,0.5,0.5,0.5,0.5) , x=c(5,4,3,2,1), color="white",size=6) + theme_classic(base_size=20)


markers_diff <- enrichr(markers[which(markers$cluster  == "Differentiated"),7], dbs_to_use)
plot_diff <- plot_gene_sets2(markers_diff[["MSigDB_Hallmark_2020"]],  plot_title = "Differentiated", num_sets=5, fill_col=col.v2[5])

labeled_diff <- plot_diff + geom_text(aes(label=markers_diff[["MSigDB_Hallmark_2020"]]$Overlap[1:5]), y=c(5,5,12,12,12) , x=c(5,4,3,2,1), color=c("white","white","black","black","black"),size=6) + theme_classic(base_size=20)


markers_p3f <- enrichr(markers[which(markers$cluster  == "PAX3-FOXO1-specific"),7], dbs_to_use)
plot_p3f <- plot_gene_sets2(markers_p3f[["MSigDB_Hallmark_2020"]],  plot_title = "PAX3::FOXO1-specific", num_sets=5, fill_col=col.v2[6])

labeled_p3f <- plot_p3f + geom_text(aes(label=markers_p3f[["MSigDB_Hallmark_2020"]]$Overlap[1:5]), y=c(1,1,1,1,1) , x=c(5,4,3,2,1), color="white",size=6) + theme_classic(base_size=20)

p4 <- (labeled_cyc + labeled_pro) / (labeled_gs + labeled_ap1) / (labeled_diff + labeled_p3f)

png("~/dosage_manuscript/figure_2/long_timept_gsea_plot.png", width = 14, height = 10, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p4)
dev.off()

#### again with chea

#Extended Data Figure 2G

plot_cyc <- plot_gene_sets2(markers_cyc[["ChEA_2022"]],  plot_title = "Cycling", num_sets=5, fill_col=col.v2[1])
labeled_cyc <- plot_cyc + geom_text(aes(label=markers_cyc[["ChEA_2022"]]$Overlap[1:5]), y=c(20,20,20,20,20) , x=c(5,4,3,2,1), color=c("white"),size=6) + theme_classic(base_size=20)

plot_pro <- plot_gene_sets2(markers_pro[["ChEA_2022"]],  plot_title = "Progenitor-like", num_sets=5, fill_col=col.v2[2])
labeled_pro <- plot_pro + geom_text(aes(label=markers_pro[["ChEA_2022"]]$Overlap[1:5]), y=c(2,2,2,2,2) , x=c(5,4,3,2,1), color="white",size=6) + theme_classic(base_size=20)

plot_gs <- plot_gene_sets2(markers_gs[["ChEA_2022"]],  plot_title = "Ground state", num_sets=5, fill_col=col.v2[3])
labeled_gs <- plot_gs + geom_text(aes(label=markers_gs[["ChEA_2022"]]$Overlap[1:5]), y=c(5,5,5,5,5) , x=c(5,4,3,2,1), color="white",size=6) + theme_classic(base_size=20)

plot_ap1 <- plot_gene_sets2(markers_ap1[["ChEA_2022"]],  plot_title = "AP-1 expressing", num_sets=5, fill_col=col.v2[4])
labeled_ap1 <- plot_ap1 + geom_text(aes(label=markers_ap1[["ChEA_2022"]]$Overlap[1:5]), y=c(10,10,10,10,10) , x=c(5,4,3,2,1), color="white",size=6) + theme_classic(base_size=20)

plot_diff <- plot_gene_sets2(markers_diff[["ChEA_2022"]],  plot_title = "Differentiated", num_sets=5, fill_col=col.v2[5])
labeled_diff <- plot_diff + geom_text(aes(label=markers_diff[["ChEA_2022"]]$Overlap[1:5]), y=c(12,12,12,12,12) , x=c(5,4,3,2,1), color="white",size=6) + theme_classic(base_size=20)

plot_p3f <- plot_gene_sets2(markers_p3f[["ChEA_2022"]],  plot_title = "PAX3::FOXO1-specific", num_sets=5, fill_col=col.v2[6])
labeled_p3f <- plot_p3f + geom_text(aes(label=markers_p3f[["ChEA_2022"]]$Overlap[1:5]), y=c(5,5,10,10,10) , x=c(5,4,3,2,1), color=c("white","white","black","black","black"),size=6) + theme_classic(base_size=20)

p6 <- (labeled_cyc + labeled_pro) / (labeled_gs + labeled_ap1) / (labeled_diff + labeled_p3f)

png("~/dosage_manuscript/figure_2/long_timept_gsea_plot_chea.png", width = 22, height = 10, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p6)
dev.off()


#Extended Data Figure 2H, right

crc <- c("MYOD1","PAX3","SOX8","RARA","PITX3","MYCN","MYOG","FOSL2","ZNF410","RELA","SNAI1","SREBF1","PKNOX2")

cycling <- c("TOP2A","CCNA2","KIF2C","MKI67","GTSE1")

progenitor <- c("TAGLN","COL3A1","COL1A2","THBS1","FN1")

ground_state <- c("PPDPF","UXT","LSM4","ROMO1","VAMP5")

ap_1 <- c("FOSB","JUN","MYCN","SOX4","RASD1")

differentiated <- c("MYBPH","MYH8","MYH3","TNNT3","LMOD2")

p3f_specific <- c("MYOD1","FGFR4","GADD45G","GAS1","PITX3")

other_markers <- unique( c( cycling, progenitor,ground_state, ap_1, differentiated, p3f_specific))

#marker dotplot
p5 <- DotPlot(dbt_2, features = other_markers, scale=T)+
        labs(x="Marker genes",y="Cluster")+
        theme(axis.text.x = element_text(angle=45, hjust=1))

png("~/dosage_manuscript/figure_2/long_timept_markers_plot_2.png", width = 12, height = 3, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p5)
dev.off()


#Extended Data Figure 2C, right

dbt_2_scores <- FetchData(dbt_2, vars=c("subclusters_2","sample_names","Progenitor_Score1","Proliferative_Score2","Differentiated_Score3")) 

p6 <- ggplot(dbt_2_scores,aes(y=Progenitor_Score1, x=subclusters_2, fill=subclusters_2))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 20)+
    scale_fill_manual(values=col.v2)+
    labs(y="Progenitor Score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.title.y=element_text(size=16))

p7 <- ggplot(dbt_2_scores,aes(y=Proliferative_Score2, x=subclusters_2, fill=subclusters_2))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 20)+
    scale_fill_manual(values=col.v2)+
    labs(y="Proliferative Score")+
    theme(axis.title.x = element_blank(),axis.text.x =element_blank(),axis.title.y=element_text(size=16))

p8 <- ggplot(dbt_2_scores,aes(y=Differentiated_Score3, x=subclusters_2, fill=subclusters_2))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 20)+
    scale_fill_manual(values=col.v2)+
    labs(y="Differentiated Score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45, vjust=1, hjust=1),axis.title.y=element_text(size=16))

png("~/dosage_manuscript/figure_2/long_timept_score_boxplots.png", width = 6, height = 10, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p6 / p7 / p8)
dev.off()
