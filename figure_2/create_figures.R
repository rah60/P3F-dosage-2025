
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

#load in short timepoints dataset
dbt_1 <- readRDS("~/dosage_manuscript/rds/short_timepts_batchcorrection_rna_namedmergedclusters_2.rds")

#load in 2 week time points dataset
dbt_2 <- readRDS("~/dosage_manuscript/rds/dosages_2wks_batchcorrection_rna_subclustered_named1.rds")

#Dimplots 8, 24 hours and 2 weeks   colors: muted tol color scheme
tol_muted <- c('#332288','#88CCEE','#44AA99','#117733','#999933','#DDCC77','#CC6677','#882255','#AA4499')

#distribution of cells across 2 wk clusters  colors: colors from figure 3
col.v1 <- viridis(8)
col.v <- col.v1[1:7] # 0, 75, 150, 250, 500, 1000, common  so 0, 75, 500 is 1, 2, 5
col_3 <- c(col.v[1],col.v[2],col.v[5])

#Figure 2A

p1 <- DimPlot(dbt_1, cols=c(tol_muted[2],tol_muted[9]), pt.size=0.8 ) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4))) 

#DimPlot(dbt_1, cols=c(tol_muted[2],tol_muted[9]), pt.size=0.8, split.by= "sample_names" ) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4))) 


png("~/dosage_manuscript/figure_2/short_timept_dimplot.png", width = 9, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()


#############


#Figure 2B

metadata <- dbt_2$subclusters_2
metadata2 <- sapply(str_replace(metadata, "PAX3-FOXO1","PAX3::FOXO1"), "[",1)
names(metadata2) <- names(metadata)
metadata2 <- factor(metadata2, levels = c("Cycling","Progenitor-like", "Ground state","AP-1 expressing","Differentiated","PAX3::FOXO1-specific"))
dbt_2$subclusters_2 <- metadata2

Idents(dbt_2) <- "subclusters_2"

p2 <- DimPlot(dbt_2, cols=c(tol_muted[1],tol_muted[3:5],tol_muted[7:8]),pt.size=0.8 ) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4)))

DimPlot(dbt_2, cols=c(tol_muted[1],tol_muted[3:5],tol_muted[7:8]),pt.size=0.8, split.by="sample_names" ) + theme_classic(base_size = 25) + guides(colour = guide_legend(override.aes = list(size=4)))


png("~/dosage_manuscript/figure_2/long_timept_dimplot.png", width = 10, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()


#Figure 2D

metadata <- dbt_2$subclusters_2
metadata2 <- sapply(str_replace(metadata, "expressing","exp."), "[",1)
metadata2 <- sapply(str_replace(metadata2, "PAX3::FOXO1","P3F"), "[",1)
metadata2 <- names(metadata2)
dbt_2$abbrev_cluster_names <- metadata2

n_cells <- FetchData(dbt_2, #look at cells per sample in each cluster
                     vars=c("abbrev_cluster_names","sample_names")) %>%
          count(abbrev_cluster_names, sample_names, name="value") #%>% #not sure why I had this here, can remove?
          #tidyr::spread(abbrev_cluster_names, n) %>%
          #pivot_longer(cols=c("Cycling","Progenitor-like", "Ground state","AP-1 exp.", "Differentiated", "P3F-specific"), names_to="cluster")

# n_totals <- FetchData(dbt_2, #look at cells per sample
#                      vars=c("abbrev_cluster_names")) %>%
#           count(abbrev_cluster_names) 

#temp test version
n_totals <- FetchData(dbt_2, #look at cells per sample
                     vars=c("sample_names")) %>%
          count(sample_names) 

n_cells$cluster_total <- NA
for(i in 1:nrow(n_cells)){

  row_num <- which(n_totals[,1] %in% n_cells[i,2])
  n_cells[i,4] <- n_totals[row_num,2]

}

n_cells$percent <- n_cells$value/n_cells$cluster_total

#temp test figure
p3 <- ggplot(n_cells, aes( x=sample_names, fill=sample_names, y=percent ))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(x="Sample",y="Fraction of sample",fill="Sample")+
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),legend.position = "none")+
  scale_y_continuous(expand = c(0,0))+
  facet_rep_wrap(vars(abbrev_cluster_names), nrow=3,scales="free" )+ #take away scales="free" to get _2 version
  scale_fill_manual(values = col_3)

png("~/dosage_manuscript/figure_2/long_timept_cluster_distribution_percent_3.png", width = 6, height = 11, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p3)
dev.off()

p3 <- ggplot(n_cells, aes( x=sample_names, fill=sample_names, y=percent ))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(x="Sample",y="Fraction of cluster",fill="Sample")+
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),legend.position = "none")+
  scale_y_continuous(expand = c(0,0))+
  facet_rep_wrap(vars(cluster), nrow=3, )+
  scale_fill_manual(values = col_3)

png("~/dosage_manuscript/figure_2/long_timept_cluster_distribution_percent.png", width = 6, height = 8, units = "in", res = 200, bg = "transparent", type = "cairo-png")
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

p4 <- cell_cycle.df %>%
  ggplot(aes( x=factor(Phase,level=c("G1 phase","S phase","G2/M phase")), fill=Sample, y=Fraction ))+
  geom_bar(stat = "identity", position = position_dodge())+
  labs(x="Cell cycle phase", y="Fraction of cells in phase")+
  theme_classic(base_size = 16)+
  scale_y_continuous(expand = c(0,0),limits=c(0,0.8))+
  scale_fill_manual(values=col_3)

png("~/dosage_manuscript/figure_2/long_timept_cell_cycle_distribution.png", width = 7, height = 4, units = "in", res = 200, bg = "transparent", type = "cairo-png")
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
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=45, vjust=1, hjust=1),axis.title.y=element_text(size=16))+
    scale_y_continuous(limits=c(min(dbt_2_scores$Differentiated_Score3),1.7))+
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, annotation=c("****","****","****"), y_position=c(1.1,1.3,1.55), show.legend=F, color = "black")

png("~/dosage_manuscript/figure_2/long_timept_scores.png", width = 3.5, height = 8, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p5 / p6 / p7)
dev.off() 


#Figure 2F
#heatmap of aggregate CRTF expression by cluster

col <- rev(brewer.pal(11,"RdBu"))

core_circ_crispr <- c("MYOD1","PAX3","SOX8","RARA","PITX3","MYCN","MYOG","FOSL2","ZNF410","RELA","SNAI1","SREBF1","PKNOX2")

core_circuitry_reorder <- c("ZNF410","FOSL2","PITX3","SOX8","PAX3","RARA","MYCN","MYOG","SREBF1","RELA","MYOD1","SNAI1","PKNOX2")

plot_dep <- data.frame(cbind(core_circuitry_reorder, c(-0.25,-0.64,-1.05,-2.81,-5.86,-1.49,-1.04,-1.03, 0.08,-0.22,-7.64,-0.09,0.45)) ) #from Gryder et al
colnames(plot_dep) <- c("gene","RH4_CRISPR_score")
plot_dep$RH4_CRISPR_score <- as.numeric(plot_dep$RH4_CRISPR_score)

aggexp <- AggregateExpression(dbt_2,group.by="abbrev_cluster_names", return.seurat = T) %>%
          FetchData( vars=core_circ_crispr, layer="scale.data")

aggexp_long <- aggexp %>% 
  rownames_to_column(var="cluster_name") %>%
  pivot_longer(cols= all_of(core_circ_crispr), names_to="gene", values_to="agg_exp")

aggexp_long <- aggexp_long %>%
  left_join(plot_dep, by="gene") %>%
  replace

p8 <- ggplot(aggexp_long)+
  geom_raster(aes(x=cluster_name,y=reorder(gene, RH4_CRISPR_score, decreasing=T),fill=agg_exp))+
  scale_fill_gradientn(colors=col, limits = c(-2.1,2.1))+
  theme_classic(base_size=20)+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0), position="top")+
  theme(axis.text.x = element_text(angle=45, hjust=-0.05), plot.margin=unit(c(0,0,0,0),"inches"), axis.title.y=element_blank())+
  labs(x="Cluster", fill="", y="")


p9 <- ggplot(plot_dep)+
      geom_col(aes(x=RH4_CRISPR_score, y=reorder(gene, RH4_CRISPR_score, decreasing=T) ))+
      labs(y="",x="RH4 CRISPR score,\n Gryder et al 2019")+
      theme_classic(base_size=20) +
      scale_y_discrete(expand = c(0,0))+
      scale_x_continuous(limits = c(-8,1))+
      theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(),axis.line.y=element_blank(),plot.margin=unit(c(0,0,0,0),"inches"),axis.title.y=element_blank())+
      geom_vline(aes(xintercept=0))+
      ylim2(p8) #to align to y axis of p8

p10 <- p8 %>% insert_left(p9, width=0.6)

png("~/dosage_manuscript/figure_2/long_timept_crtf_heatmap.png", width = 10, height = 10, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p10)
dev.off()

