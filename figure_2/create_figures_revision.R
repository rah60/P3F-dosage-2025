
rm(list = ls())
graphics.off()

#creates Figure 2 panels

library(tidyverse)
library(Seurat)
library(Cairo)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(aplot)
library(rstatix)
library(lemon)
library(qs)

#load in short timepoints dataset
dbt_1 <- qread("/home/gdstantonlab/lab/sharing/Ambuj/PAX3_FOXO1/2026/PAX3_FOXO1_Rachel_March13/Phase_8_24_hrs/seurat_induction_8_24_hrs_annotated.qs")

#load in 2 week time points dataset
dbt_2 <- qread("/home/gdstantonlab/lab/sharing/Ambuj/PAX3_FOXO1/2026/PAX3_FOXO1_Rachel_March13/Phase1/seurat_phase1_annotated.qs")
dbt_2$sample_names <- factor(dbt_2$sample_names, levels= c("0 ng/mL, 2 wks","75 ng/mL, 2 wks","500 ng/mL, 2 wks"))

#Dimplots 8, 24 hours and 2 weeks   colors: muted tol color scheme
tol_muted <- c('#332288','#88CCEE','#44AA99','#117733','#999933','#DDCC77','#CC6677','#882255','#AA4499')

#distribution of cells across 2 wk clusters  colors: colors from figure 3
col.v1 <- viridis(8)
col.v <- col.v1[1:7] # 0, 75, 150, 250, 500, 1000, common  so 0, 75, 500 is 1, 2, 5
col_3 <- c(col.v[1],col.v[2],col.v[5])

#Figure 2A

p1 <- DimPlot(dbt_1, cols=c(tol_muted[2],tol_muted[9]), pt.size=0.8 ) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4))) 

#DimPlot(dbt_1, cols=c(tol_muted[2],tol_muted[9]), pt.size=0.8, split.by= "sample_names" ) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4))) 


png("~/dosage_manuscript/figure_2/short_timept_dimplot_revised.png", width = 9, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()


#############


#Figure 2B
Idents(dbt_2) <- "cell_label"

metadata <- dbt_2$cell_label
metadata2 <- sapply(str_replace(metadata, "PAX3-FOXO1","PAX3::FOXO1"), "[",1)
metadata2 <- sapply(str_replace(metadata2, "Ground-state","Ground state"), "[",1)
names(metadata2) <- names(metadata)
metadata2 <- factor( metadata2, levels= c("Progenitor-like","Cycling","Ground state","PAX3::FOXO1-specific","Differentiated"))
dbt_2$subclusters_3 <- metadata2

Idents(dbt_2) <- "subclusters_3"

#correct order of colors
cols_dbt2 <-  c(tol_muted[3], tol_muted[1],tol_muted[4],tol_muted[8],tol_muted[7])

p2 <- DimPlot(dbt_2, cols=cols_dbt2, pt.size=0.8 ) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4)))

#DimPlot(dbt_2, cols=cols_dbt2, pt.size=0.8, split.by="sample_names" ) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4)))


png("~/dosage_manuscript/figure_2/long_timept_dimplot_revised.png", width = 10, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()


#Figure 2D

metadata <- dbt_2$subclusters_3
metadata2 <- sapply(str_replace(metadata, "expressing","exp."), "[",1)
metadata2 <- sapply(str_replace(metadata2, "PAX3::FOXO1","P3F"), "[",1)
metadata2 <- names(metadata2)
dbt_2$abbrev_cluster_names <- metadata2

n_cells <- FetchData(dbt_2, #look at cells per sample in each cluster
                     vars=c("abbrev_cluster_names","sample_names")) %>%
          count(abbrev_cluster_names, sample_names, name="value") 


n_totals <- FetchData(dbt_2, #look at cells per sample
                     vars=c("sample_names")) %>%
          count(sample_names) 

n_cells$cluster_total <- NA
for(i in 1:nrow(n_cells)){

  row_num <- which(n_totals[,1] %in% n_cells[i,2])
  n_cells[i,4] <- n_totals[row_num,2]

}

n_cells$percent <- n_cells$value/n_cells$cluster_total

#for source data
write.csv(n_cells, "~/dosage_manuscript/csv/n_cells_figure_2d.csv")

p3 <- ggplot(n_cells, aes( x=sample_names, fill=sample_names, y=percent ))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(x="Doxycycline (ng/mL),\n2 week induction",y="Fraction of sample",fill="Sample")+
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1,size=18),legend.position = "none")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = 0 ) )+ #don't expand y-axis on bottom, expand y-axis on top by 10%
  facet_wrap(vars(abbrev_cluster_names), nrow=3,scales="free_y" )+ 
  scale_fill_manual(values = col_3)#+
  #theme(panel.spacing = unit(0, "lines"))

png("~/dosage_manuscript/figure_2/long_timept_cluster_distribution_percent_revised2.png", width = 6.4, height = 7.5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
#pdf("~/dosage_manuscript/figure_2/long_timept_cluster_distribution_percent_revised2.pdf", width = 8.5, height = 10, bg = "transparent")
print(p3)
dev.off()

#Figure 2C

#Cell cycle phase distribution at 2 weeks    colors: colors from figure 3
metadata <- dbt_2@meta.data
dbt_names <- as.vector(unique(metadata$sample_names))
cell_cycle.df <- tibble()
                                                             
for(i in 1:length(dbt_names)){
  S_count <- length(which(metadata$Phase == "S" & metadata$sample_names == dbt_names[i]))
  G2M_count <- length(which(metadata$Phase == "G2M" & metadata$sample_names == dbt_names[i]))
  G1_count <- length(which(metadata$Phase == "G1" & metadata$sample_names == dbt_names[i]))
  total <- sum(S_count, G2M_count, G1_count)
  cell_cycle.df <- rbind(cell_cycle.df, c(S_count/total, dbt_names[i], "S phase"), c(G2M_count/total,dbt_names[i], "G2/M phase"), c(G1_count/total,dbt_names[i], "G1 phase"))
  
}
colnames(cell_cycle.df) <- c("Fraction", "Sample","Phase")
cell_cycle.df$Fraction <- as.numeric(cell_cycle.df$Fraction)
cell_cycle.df$Sample <- factor(cell_cycle.df$Sample, levels= c("0 ng/mL, 2 wks","75 ng/mL, 2 wks","500 ng/mL, 2 wks"))
new_sample <- c("0 ng/mL, 2 wks" = "0","75 ng/mL, 2 wks" = "75","500 ng/mL, 2 wks" = "500")

#for source data
write.csv(cell_cycle.df, "~/dosage_manuscript/csv/cell_cycle_two_week.csv")

p4 <- cell_cycle.df %>%
  ggplot(aes( x=factor(Phase,level=c("G1 phase","S phase","G2/M phase")), fill=Sample, y=Fraction ))+
  geom_bar(stat = "identity", position = position_dodge())+
  labs(x="Cell cycle phase", y="Fraction of cells in phase")+
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(size=14))+
  scale_y_continuous(expand = c(0,0),limits=c(0,0.8))+
  scale_x_discrete(labels=c("G1","S","G2/M"))+
  scale_fill_manual(values=col_3)+
  facet_wrap(vars(Sample), ncol=3, labeller=labeller(Sample=new_sample))

png("~/dosage_manuscript/figure_2/long_timept_cell_cycle_distribution_revised_2.png", width = 7.5, height = 4, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p4)
dev.off()


#Figure 2E

#scores from Danielli et al    colors: colors from figure 3

#make a dataframe to pull out the data points for each score for each cell so I can test for significance
dbt_2_scores <- FetchData(dbt_2, vars=c("ident","sample_names","Progenitor_Score1","Proliferative_Score2","Differentiated_Score3")) 

#to check significance
aov_res <- aov(Progenitor_Score1 ~ sample_names, data = dbt_2_scores) %>% tukey_hsd()

aov_res2 <- aov(Proliferative_Score2 ~ sample_names, data = dbt_2_scores) %>% tukey_hsd()

aov_res3 <- aov(Differentiated_Score3 ~ sample_names, data = dbt_2_scores) %>% tukey_hsd()

my_comparisons <- list(c("0 ng/mL, 2 wks", "75 ng/mL, 2 wks"), c("0 ng/mL, 2 wks","500 ng/mL, 2 wks"),c("75 ng/mL, 2 wks","500 ng/mL, 2 wks"))

p5 <- ggplot(dbt_2_scores,aes(y=Progenitor_Score1, x=sample_names, color=sample_names))+
    geom_boxplot(show.legend=F, notch=T)+
    theme_classic(base_size = 20)+
    scale_color_manual(values=col_3)+
    labs(y="Progenitor Score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.title.y=element_text(size=16))+
    scale_y_continuous(limits=c(min(dbt_2_scores$Progenitor_Score1),1.25))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=c("****","****","****"), y_position=c(0.75,0.9,1.1), show.legend=F, color = "black")

p6 <- ggplot(dbt_2_scores,aes(y=Proliferative_Score2, x=sample_names, color=sample_names))+
    geom_boxplot(show.legend=F, notch=T)+
    theme_classic(base_size = 20)+
    scale_color_manual(values=col_3)+
    labs(y="Proliferative Score")+
    theme(axis.title.x = element_blank(),axis.text.x =element_blank(),axis.title.y=element_text(size=16))+
    scale_y_continuous(limits=c(min(dbt_2_scores$Proliferative_Score2),1.5))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=c("****","****","****"), y_position=c(0.95,1.15,1.40), show.legend=F, color = "black")

p7 <- ggplot(dbt_2_scores,aes(y=Differentiated_Score3, x=sample_names, color=sample_names))+
    geom_boxplot(show.legend=F, notch=T)+
    theme_classic(base_size = 20)+
    scale_color_manual(values=col_3)+
    labs(y="Differentiated Score")+
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45, vjust=1, hjust=1,size=18),axis.title.y=element_text(size=16))+
    scale_y_continuous(limits=c(min(dbt_2_scores$Differentiated_Score3),1.7))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=c("****","****","****"), y_position=c(1.1,1.3,1.55), show.legend=F, color = "black")

png("~/dosage_manuscript/figure_2/long_timept_scores_revised.png", width = 3.5, height = 8, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p5 / p6 / p7)
dev.off() 


