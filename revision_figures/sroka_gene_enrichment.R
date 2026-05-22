rm(list = ls())
graphics.off()

#analyzes DE genes in Sroka et al data for Supplemental Figures 13 C & E

library(DESeq2)
library(dplyr)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(fgsea)
library(msigdbr)

#from Sroka et al 2023 supplemental data
RH4_d7 <- read_csv("~/dosage_manuscript/csv/sroka_rh4_d7_sgP3F.csv") %>%
  mutate(across(c(log2FoldChange), as.numeric)) %>%
  dplyr::filter(!is.na(log2FoldChange))

RH41_d7 <- read_csv("~/dosage_manuscript/csv/sroka_RH41_d7_sgP3F.csv") %>%
  mutate(across(c(log2FoldChange), as.numeric)) %>%
  dplyr::filter(!is.na(log2FoldChange))

####### RH4

h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") #human hallmark gene sets: H. GO BP: C5
h_gene_sets2 <- split(x=h_gene_sets$gene_symbol, f=h_gene_sets$gs_name) 

markers.v <- RH4_d7$log2FoldChange #avglogFC

names(markers.v) <- RH4_d7$Gene
  
fgsea.l <- fgsea(pathways = h_gene_sets2, stats=markers.v , maxSize=500, minSize=5, nPerm=1000) %>% #or both & run correct FGSEA settings
  arrange(-NES)

# plot top up/down-regulated pathways, S13C
plot_pathways <- fgsea.l[which(fgsea.l$padj < 0.05),] %>%
                arrange(NES)

store_names <- strsplit(plot_pathways$pathway, "HALLMARK_") 
store_names <- sapply(store_names,"[[",2)
store_names <- str_replace_all(store_names, "_"," ")
plot_pathways$pathway <- store_names

plot_pathways$pathway <- factor(plot_pathways$pathway, levels=plot_pathways$pathway)

p1 <- ggplot(plot_pathways)+
  geom_count(aes(y=pathway,x=NES,color=padj))+
  labs(y="",x="Normalized Enrichment Score",title = "Enriched Hallmark Gene Sets\n in RH4 +sgP3F",color="Adjusted\np-value")+
  scale_color_viridis_c(begin=0.8, end=0, option="turbo")+
  theme_classic(base_size=12)+
  guides(size="none")

png( paste0("~/dosage_manuscript/revision_figures/sroka_rh4_gsea.png"), width = 6, height = 5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()

########## RH41
markers.v <- RH41_d7$log2FoldChange #avglogFC

names(markers.v) <- RH41_d7$Gene

fgsea.l <- fgsea(pathways = h_gene_sets2, stats=markers.v , maxSize=500, minSize=5, nPerm=1000) %>% #or both & run correct FGSEA settings
  arrange(-NES)

# plot top up/down-regulated pathways, S13E
plot_pathways <- fgsea.l[which(fgsea.l$padj < 0.05),] %>%
  arrange(NES)

store_names <- strsplit(plot_pathways$pathway, "HALLMARK_") 
store_names <- sapply(store_names,"[[",2)
store_names <- str_replace_all(store_names, "_"," ")
plot_pathways$pathway <- store_names

plot_pathways$pathway <- factor(plot_pathways$pathway, levels=plot_pathways$pathway)

p1 <- ggplot(plot_pathways)+
  geom_count(aes(y=pathway,x=NES,color=padj))+
  labs(y="",x="Normalized Enrichment Score",title = "Enriched Hallmark Gene Sets\n in RH41 +sgP3F",color="Adjusted\np-value")+
  scale_color_viridis_c(begin=0.8, end=0, option="turbo")+
  theme_classic(base_size=12)+
  guides(size="none")

png( paste0("~/dosage_manuscript/revision_figures/sroka_rh41_gsea.png"), width = 6, height = 5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()
