rm(list = ls())
graphics.off()

#plots P3F withdrawal figures added during revision, in Figure 2 and Supplemental Figure 5

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
library(lemon)
library(aplot)
library(rstatix)
library(ggpubr)

source("~/dosage_manuscript/figure_2/figure_2_functions.R")

tol_muted <- c('#332288','#88CCEE','#44AA99','#117733','#999933','#DDCC77','#CC6677','#882255','#AA4499')

#distribution of cells across 2 wk clusters  colors: colors from figure 3
col.v1 <- viridis(8)
col.v <- col.v1[1:7] # 0, 75, 150, 250, 500, 1000, common  so 0, 75, 500 is 1, 2, 5
col_3 <- c(col.v[1],col.v[2],col.v[5])

col.v2 <- c(tol_muted[1],tol_muted[3],tol_muted[8],tol_muted[7])

#from Ambuj
dbt_3 <- qread("/home/gdstantonlab/lab/sharing/Ambuj/PAX3_FOXO1/2026/PAX3_FOXO1_Rachel_March13/Phase2/seurat_phase2_annotated.qs")

metadata <- dbt_3$cell_label
metadata2 <- sapply(str_replace(metadata, "PAX3-FOXO1","PAX3::FOXO1"), "[",1)
names(metadata2) <- names(metadata)
metadata2 <- factor( metadata2, levels= c("Cycling","Progenitor-like","PAX3::FOXO1-specific","Differentiated"))
dbt_3$subclusters_3 <- metadata2

Idents(dbt_3) <- "subclusters_3"


####### dimplots

####### main fig

p1 <- DimPlot(dbt_3, pt.size=0.8, cols = col.v2) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4))) + labs(title="")

png("~/dosage_manuscript/figure_2/p3f_withdrawal_dimplot.png", width = 9, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()

######## supplemental fig

p2 <- DimPlot(dbt_3, pt.size=0.8, group.by="Phase") + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4))) + labs(title="")

png("~/dosage_manuscript/figure_2/p3f_withdrawal_phase_dimplot.png", width = 7, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()

###### plot distribution across timepoints

metadata <- dbt_3$orig.ident
metadata <- factor(gsub("SEG0503", "0 hrs", metadata))
metadata <- factor(gsub("SEG0504", "48 hrs", metadata))
metadata <- factor(gsub("SEG0505", "96 hrs", metadata))
names(metadata) <- names(dbt_3$orig.ident)

dbt_3$sample_names <- metadata

metadata <- dbt_3$subclusters_3
metadata2 <- sapply(str_replace(metadata, "expressing","exp."), "[",1)
metadata2 <- sapply(str_replace(metadata2, "PAX3::FOXO1","P3F"), "[",1)
metadata2 <- names(metadata2)
dbt_3$abbrev_cluster_names <- metadata2

######## main fig

n_cells <- FetchData(dbt_3, #look at cells per sample in each cluster
                     vars=c("abbrev_cluster_names","sample_names")) %>%
          count(abbrev_cluster_names, sample_names, name="value") 


n_totals <- FetchData(dbt_3, #look at cells per sample
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
  labs(x="Hours after PAX3::FOXO1 removal",y="Fraction of sample",fill="Sample")+
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1,size=20),legend.position = "none")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = 0 ) )+ #don't expand y-axis on bottom, expand y-axis on top by 10%
  facet_rep_wrap(vars(abbrev_cluster_names), nrow=1,scales="free" )+ 
  scale_fill_manual(values = col_3)

png("~/dosage_manuscript/figure_2/withdrawal_cluster_distribution_percent.png", width = 10, height =3.33, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p3)
dev.off()

####### proliferative, progenitor, differentiated scores by sample & by cluster, supplemental figures

#set up for module scores
danielli_markers <- read_csv("~/Dbt_scRNA/danielli_2023_markers.csv")
progenitor_markers <- as.vector(danielli_markers$Progenitor_markers)
progenitor_markers  <- progenitor_markers[!is.na(progenitor_markers)]
proliferation_markers <- as.vector(danielli_markers$Proliferative_markers)
proliferation_markers  <- proliferation_markers[!is.na(proliferation_markers)]
differentiated_markers <- as.vector(danielli_markers$Differentiated_markers) 
differentiated_markers  <- differentiated_markers[!is.na(differentiated_markers)]

dbt_3 <- AddModuleScore(object=dbt_3, features=list(progenitor_markers, proliferation_markers,differentiated_markers), name=c("Progenitor_Score", "Proliferative_Score","Differentiated_Score"), assay="RNA", search=T)

dbt_3_scores <- FetchData(dbt_3, vars=c("ident","sample_names","Progenitor_Score1","Proliferative_Score2","Differentiated_Score3")) 

aov_res <- aov(Progenitor_Score1 ~ sample_names, data = dbt_3_scores) %>% tukey_hsd()

aov_res2 <- aov(Proliferative_Score2 ~ sample_names, data = dbt_3_scores) %>% tukey_hsd()

aov_res3 <- aov(Differentiated_Score3 ~ sample_names, data = dbt_3_scores) %>% tukey_hsd()

my_comparisons <- list(c("0 hrs", "48 hrs"), c("0 hrs","96 hrs"),c("48 hrs","96 hrs"))

p5 <- ggplot(dbt_3_scores,aes(y=Progenitor_Score1, x=sample_names, color=sample_names))+
    geom_boxplot(show.legend=F, notch=T)+
    theme_classic(base_size = 20)+
    scale_color_manual(values=col_3)+
    labs(y="Progenitor Score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.title.y=element_text(size=16))+
    scale_y_continuous(limits=c(min(dbt_3_scores$Progenitor_Score1),0.85))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=c("****","****","****"), y_position=c(0.55,0.65,0.8), show.legend=F, color = "black")

p6 <- ggplot(dbt_3_scores,aes(y=Proliferative_Score2, x=sample_names, color=sample_names))+
    geom_boxplot(show.legend=F, notch=T)+
    theme_classic(base_size = 20)+
    scale_color_manual(values=col_3)+
    labs(y="Proliferative Score")+
    theme(axis.title.x = element_blank(),axis.text.x =element_blank(),axis.title.y=element_text(size=16))+
    scale_y_continuous(limits=c(min(dbt_3_scores$Proliferative_Score2),1.5))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=c("****","****","****"), y_position=c(0.95,1.15,1.35), show.legend=F, color = "black")

p7 <- ggplot(dbt_3_scores,aes(y=Differentiated_Score3, x=sample_names, color=sample_names))+
    geom_boxplot(show.legend=F, notch=T)+
    theme_classic(base_size = 20)+
    scale_color_manual(values=col_3)+
    labs(y="Differentiated Score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45, vjust=1, hjust=1),axis.title.y=element_text(size=16))+
    scale_y_continuous(limits=c(min(dbt_3_scores$Differentiated_Score3),1.7))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=c("****","****","****"), y_position=c(1,1.2,1.4), show.legend=F, color = "black")

png("~/dosage_manuscript/figure_2/withdrawal_scores_bysample.png", width = 3.5, height = 8, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p5 / p6 / p7)
dev.off() 


##### by cluster
dbt_3_scores2 <- FetchData(dbt_3, vars=c("subclusters_3","sample_names","Progenitor_Score1","Proliferative_Score2","Differentiated_Score3")) 

aov_res <- aov(Progenitor_Score1 ~ subclusters_3, data = dbt_3_scores2) %>% tukey_hsd()

aov_res2 <- aov(Proliferative_Score2 ~ subclusters_3, data = dbt_3_scores2) %>% tukey_hsd()

aov_res3 <- aov(Differentiated_Score3 ~ subclusters_3, data = dbt_3_scores2) %>% tukey_hsd()

my_comparisons <- list(c("PAX3::FOXO1-specific","Differentiated"),c("Progenitor-like","Cycling"),
        c("Progenitor-like","PAX3::FOXO1-specific"),c("Cycling","PAX3::FOXO1-specific"),c("Progenitor-like","Differentiated"), 
        c("Cycling","Differentiated"))
sig_levels <- c("****","****","****","****","****","****")

p8 <- ggplot(dbt_3_scores2,aes(y=Progenitor_Score1, x=subclusters_3, fill=subclusters_3))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 20)+
    scale_fill_manual(values=col.v2)+
    labs(y="Progenitor Score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.title.y=element_text(size=16))+
    scale_y_continuous(limits=c(min(dbt_3_scores2$Progenitor_Score1),1))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=sig_levels, y_position=c(0.45,0.55,0.6,0.7,0.85,0.95), show.legend=F, color = "black",textsize = 3)


p9 <- ggplot(dbt_3_scores2,aes(y=Proliferative_Score2, x=subclusters_3, fill=subclusters_3))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 20)+
    scale_fill_manual(values=col.v2)+
    labs(y="Proliferative Score")+
    theme(axis.title.x = element_blank(),axis.text.x =element_blank(),axis.title.y=element_text(size=16))+
    scale_y_continuous(limits=c(min(dbt_3_scores2$Proliferative_Score2),1.6))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=sig_levels, y_position=c(0.8,0.9,1,1.15,1.35,1.5), show.legend=F, color = "black",textsize = 3)


p10 <- ggplot(dbt_3_scores2,aes(y=Differentiated_Score3, x=subclusters_3, fill=subclusters_3))+
    geom_boxplot(show.legend=F, notch=T, linewidth=1.2)+
    theme_classic(base_size = 20)+
    scale_fill_manual(values=col.v2)+
    labs(y="Differentiated Score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45, vjust=1, hjust=1),axis.title.y=element_text(size=16))+
    scale_y_continuous(limits=c(min(dbt_3_scores2$Differentiated_Score3),1.6))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=sig_levels, y_position=c(0.9,0.8,1,1.15,1.35,1.5), show.legend=F, color = "black",textsize = 3)


png("~/dosage_manuscript/figure_2/withdrawal_scores_byclusters.png", width = 6, height = 10, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p8 / p9 / p10)
dev.off()
