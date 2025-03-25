rm(list = ls())
graphics.off()

#documents the iterative manual process of annotating the single cell clusters

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

source("~/dosage_manuscript/figure_2/figure_2_functions.R")

#read in merged, processed Seurat object
#from analysis_scrna_RAH_modified.qmd, originally by Meng Wang
dbt <- readRDS("/gpfs0/home2/gdstantonlab/lab/scrna_analysis/all_samples_merged.rds")

#make sample names more readable
metadata <- dbt$condition
metadata2 <- sapply(strsplit(metadata, "Dbt-"), "[",2)
metadata2 <- sapply(str_replace(metadata2, "-"," ng/mL, "), "[",1)
metadata2 <- sapply(str_replace(metadata2, "hr"," hrs"), "[",1)
metadata2 <- sapply(str_replace(metadata2, "wk"," wks"), "[",1)
metadata2 <-  names(metadata2)
dbt$sample_names <- metadata2

#set up for module scores
danielli_markers <- read_csv("~/Dbt_scRNA/danielli_2023_markers.csv")
progenitor_markers <- as.vector(danielli_markers$Progenitor_markers)
progenitor_markers  <- progenitor_markers[!is.na(progenitor_markers)]
proliferation_markers <- as.vector(danielli_markers$Proliferative_markers)
proliferation_markers  <- proliferation_markers[!is.na(proliferation_markers)]
differentiated_markers <- as.vector(danielli_markers$Differentiated_markers) 
differentiated_markers  <- differentiated_markers[!is.na(differentiated_markers)]

#divide data into subsets
Idents(dbt) <- "sample_names"
dbt_1 <- subset(dbt, idents=c("0 ng/mL, 0 hrs","500 ng/mL, 8 hrs","500 ng/mL, 24 hrs") )

dbt_2 <- subset(dbt, idents=c("0 ng/mL, 2 wks","75 ng/mL, 2 wks","500 ng/mL, 2 wks") )

#uncomment below to run processing

#process_scrna_norm(dbt_1, "short_timepoints_new_aligned")
#dbt_1 <- readRDS("~/Dbt_scRNA/short_timepoints_new_aligned.rds")

#process_scrna_norm(dbt_2, "dosages_2wks_new_aligned")
#dbt_2 <- readRDS("~/Dbt_scRNA/dosages_2wks_new_aligned.rds")

#############
#short timepoints batch correction
dbt_1 <- readRDS("~/Dbt_scRNA/short_timepoints_new_aligned.rds")

#run module scores
dbt_1 <- AddModuleScore(object=dbt_1, features=list(progenitor_markers, proliferation_markers,differentiated_markers), name=c("Progenitor_Score", "Proliferative_Score","Differentiated_Score"), assay="RNA", search=T)
dbt_1$muscle_lineage_score <- dbt_1$Differentiated_Score3 - dbt_1$Progenitor_Score1

#factor names for plotting purposes
dbt_1$sample_names <- factor(dbt_1$sample_names, levels= c("0 ng/mL, 0 hrs","500 ng/mL, 8 hrs","500 ng/mL, 24 hrs"))

dbt_1 <- RunHarmony(dbt_1, #batch correct
    group.by.vars = "sample_names",
    dims.use = 1:30) %>%
    RunUMAP(reduction = "harmony", dims = 1:30) #and then re-run UMAP

dbt_1 <- FindNeighbors(dbt_1, reduction="harmony", dims=1:30) #re-run clustering on harmony reduction
dbt_1 <- FindClusters(dbt_1, resolution=c(0.1, 0.2, 0.3))

clustree(dbt_1) #examine resolutions to pick optimal one
Idents(dbt_1) <- "RNA_snn_res.0.2"

#stored intermediate
#dbt_1 <- readRDS("~/Dbt_scRNA/short_timepts_batchcorrection_rna.rds")

VlnPlot(dbt_1, features=c("Progenitor_Score1", "Proliferative_Score2","Differentiated_Score3"), pt.size=0, plasma(n=8, begin=0, end=0.8)) #look at scores across clusters

FeaturePlot(dbt_1, features=c("Progenitor_Score1"),order=T) + scale_color_viridis_c()
FeaturePlot(dbt_1, features=c("Proliferative_Score2"),order=T) + scale_color_viridis_c()
FeaturePlot(dbt_1, features=c("Differentiated_Score3"),order=T) + scale_color_viridis_c()

markers <- FindAllMarkers(dbt_1, logfc.threshold = 0.25, min.pct=0.1,  min.diff.pct=0.1) #call markers across clusters
#saveRDS(markers, "~/Dbt_scRNA/markers_short_timepts_batchcorrection_rna.rds" ) #saved intermediate
markers <- readRDS("~/Dbt_scRNA/markers_short_timepts_batchcorrection_rna.rds")

gsea_res <- run_fgsea(markers) #gsea on markers

DimPlot(dbt_1, split.by = "Phase") #look at clusters by cell cycle phase

new.cluster.ids <- c("ground state","cycling","cycling","?") #initial annotations of clusters
names(new.cluster.ids) <- levels(dbt_1)
dbt_1 <- RenameIdents(dbt_1, new.cluster.ids)
dbt_1$clusters <- Idents(dbt_1)
Idents(dbt_1) <- "clusters"

dbt_1 <- readRDS("~/Dbt_scRNA/short_timepts_batchcorrection_rna_tempclusters.rds")

#try subclustering the ? and ground state clusters
cluster_1 <- subset(dbt_1, ident="?")
cluster_2 <- subset(dbt_1, ident="ground state")

cluster_1 <- FindNeighbors(cluster_1, reduction="harmony", dims=1:30) #re-run clustering on harmony reduction
cluster_1 <- FindClusters(cluster_1, resolution=c(0.1, 0.2, 0.3))

clustree(cluster_1)
Idents(cluster_1) <- "RNA_snn_res.0.2"
DimPlot(cluster_1)
markers <- FindAllMarkers(cluster_1, logfc.threshold = 0.25, min.pct=0.1,  min.diff.pct=0.1)

#gsea on subclusters from ? cluster
gsea_res <- run_fgsea(markers)
#0 cycling
#1 NS
#2 TNFa/p53

#ground state subclustering
cluster_2 <- FindNeighbors(cluster_2, reduction="harmony", dims=1:30) #re-run clustering on harmony reduction
cluster_2 <- FindClusters(cluster_2, resolution=c(0.1, 0.2, 0.3))

clustree(cluster_2)
Idents(cluster_2) <- "RNA_snn_res.0.2"
DimPlot(cluster_2)
markers <- FindAllMarkers(cluster_2, logfc.threshold = 0.25, min.pct=0.1,  min.diff.pct=0.1)
gsea_res <- run_fgsea(markers)
#0: NS
#1: NS
#2: NS
#3: Glycolysis, TNFa, Kras

#split off cycling cluster from ? cluster, determined the rest are ground state/little significant gsea terms
new.cluster.ids <- c("cycling","ground state","ground state")
names(new.cluster.ids) <- levels(cluster_1)
cluster_1 <- RenameIdents(cluster_1, new.cluster.ids)

metadata <- dbt_1@meta.data
clusters <- droplevels(metadata$clusters)
names(clusters) <- rownames(metadata)

clusters[Cells(subset(cluster_1,ident="cycling"))] <- "cycling"
clusters[Cells(subset(cluster_1,ident="ground state"))] <- "ground state"

dbt_1$clusters <- clusters 
Idents(dbt_1) <- "clusters"

#saveRDS(dbt_1, "~/Dbt_scRNA/short_timepts_batchcorrection_rna_namedmergedclusters_2.rds") #saved intermediate

dbt_1 <- readRDS("~/Dbt_scRNA/short_timepts_batchcorrection_rna_namedmergedclusters_2.rds")

#find markers, run gsea on final clusters
markers <- FindAllMarkers(dbt_1, logfc.threshold = 0.25, min.pct=0.1,  min.diff.pct=0.1)
gsea_res <- run_fgsea(markers, perm_num = 10000)



#############
#dosage 2 wks batch correction

dbt_2 <- readRDS("~/Dbt_scRNA/dosages_2wks_batchcorrection_rna.rds")

#add module scores from Danielli et al
dbt_2 <- AddModuleScore(object=dbt_2, features=list(progenitor_markers, proliferation_markers,differentiated_markers), name=c("Progenitor_Score", "Proliferative_Score","Differentiated_Score"), assay="RNA", search=T)
dbt_2$muscle_lineage_score <- dbt_2$Differentiated_Score3 - dbt_2$Progenitor_Score1

#factor sample names
dbt_2$sample_names <- factor(dbt_2$sample_names, levels= c("0 ng/mL, 2 wks","75 ng/mL, 2 wks","500 ng/mL, 2 wks"))

dbt_2 <- RunHarmony(dbt_2, #batch correct
    group.by.vars = "sample_names",
    dims.use = 1:30) %>%
    RunUMAP(reduction = "harmony", dims = 1:30)

dbt_2 <- FindNeighbors(dbt_2, reduction="harmony", dims=1:30) #re-run clustering on harmony reduction
dbt_2 <- FindClusters(dbt_2, resolution=c(0.1, 0.2, 0.3))

#determine clustering resolution
clustree(dbt_2)
Idents(dbt_2) <- "RNA_snn_res.0.1"
DimPlot(dbt_2, split.by="sample_names")

#examine module scores across clusters
VlnPlot(dbt_2, features=c("Progenitor_Score1", "Proliferative_Score2","Differentiated_Score3"), pt.size=0, plasma(n=7, begin=0, end=0.8)) 
FeaturePlot(dbt_2, features=c("Differentiated_Score3"),order=T) + scale_color_viridis_c()

#cluster markers
markers <- FindAllMarkers(dbt_2, logfc.threshold = 0.25, min.pct=0.1,  min.diff.pct=0.1)
#saveRDS(markers, "~/Dbt_scRNA/markers_dosages_2wks_batchcorrection_rna_02.rds" ) #saved intermediate
markers <- readRDS("~/Dbt_scRNA/markers_dosages_2wks_batchcorrection_rna.rds")

top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

#heatmap
p1 <- DoHeatmap(subset(dbt_2, downsample=300), features = unique(top10$gene), group.by="ident", group.color=magma(7)[1:6]) +
      scale_fill_viridis()
p1

VlnPlot(dbt_2, features=c("JUN")) #examine AP1 TF expression across clusters

gsea_res <- run_fgsea(markers)

FeaturePlot(dbt_2, features="muscle_lineage_score",order=T) + scale_color_viridis_c()

FeaturePlot(dbt_2, features=c("Progenitor_Score1"),order=T) + scale_color_viridis_c()
FeaturePlot(dbt_2, features=c("Proliferative_Score2"),order=T) + scale_color_viridis_c()
FeaturePlot(dbt_2, features=c("Differentiated_Score3"),order=T) + scale_color_viridis_c()


#saveRDS(dbt_2, "~/Dbt_scRNA/dosages_2wks_batchcorrection_rna_subcluster.rds" )

#subcluster cluster 1 & 2

cluster_1 <- subset(dbt_2, ident=1)
cluster_2 <- subset(dbt_2, ident=2)

cluster_1 <- FindNeighbors(cluster_1, reduction="harmony", dims=1:30) #re-run clustering on harmony reduction
cluster_1 <- FindClusters(cluster_1, resolution=c(0.1, 0.2, 0.3))

clustree(cluster_1)
Idents(cluster_1) <- "RNA_snn_res.0.1"
DimPlot(cluster_1)

markers <- FindAllMarkers(cluster_1, logfc.threshold = 0.25, min.pct=0.1,  min.diff.pct=0.1)
#saveRDS(markers, "~/Dbt_scRNA/markers_dosages_2wks_batchcorrection_rna_cluster1.rds" ) #saved intermediate
markers <- readRDS("~/Dbt_scRNA/markers_dosages_2wks_batchcorrection_rna_cluster1.rds")

gsea_res <- run_fgsea(markers)

#0: NS? combine with 1 perhaps
#1: EMT (Progenitor-like)
#2: Myogenesis

#cluster 2 subcluster
cluster_2 <- FindNeighbors(cluster_2, reduction="harmony", dims=1:30) #re-run clustering on harmony reduction
cluster_2 <- FindClusters(cluster_2, resolution=c(0.1, 0.2, 0.3))

clustree(cluster_2)
Idents(cluster_2) <- "RNA_snn_res.0.1"
DimPlot(cluster_2)

markers <- FindAllMarkers(cluster_2, logfc.threshold = 0.25, min.pct=0.1,  min.diff.pct=0.1)
#saveRDS(markers, "~/Dbt_scRNA/markers_dosages_2wks_batchcorrection_rna_cluster2.rds" ) #saved intermediate
markers <- readRDS( "~/Dbt_scRNA/markers_dosages_2wks_batchcorrection_rna_cluster2.rds")

gsea_res <- run_fgsea(markers)
#0: EMT (Progenitor)
#1: NS

#splitting out subclusters from cluster 1
new.cluster.ids <- c("1","5","6")
names(new.cluster.ids) <- levels(cluster_1)
cluster_1 <- RenameIdents(cluster_1, new.cluster.ids)

#splitting out subclusters from cluster 2
new.cluster.ids <- c("2","7")
names(new.cluster.ids) <- levels(cluster_2)
cluster_2 <- RenameIdents(cluster_2, new.cluster.ids)

new.cluster.ids <- c("0","1","5","6","2","7","3","4")
metadata <- dbt_2@meta.data
clusters <- droplevels(metadata$RNA_snn_res.0.1)
names(clusters) <- rownames(metadata)
levels(clusters) <- c("0","1","2","3","4","5","6","7")

clusters[Cells(subset(cluster_1,ident=1))] <- 1
clusters[Cells(subset(cluster_1,ident=5))] <- 5
clusters[Cells(subset(cluster_1,ident=6))] <- 6
clusters[Cells(subset(cluster_2,ident=2))] <- 2
clusters[Cells(subset(cluster_2,ident=7))] <- 7

dbt_2$clusters <- clusters
Idents(dbt_2) <- "clusters"
DimPlot(dbt_2)
#saveRDS(dbt_2, "~/Dbt_scRNA/dosages_2wks_batchcorrection_rna_subcluster.rds") #saved intermediate

dbt_2 <- readRDS("~/Dbt_scRNA/dosages_2wks_batchcorrection_rna_subcluster.rds")
DimPlot(dbt_2)

#look at markers with subclusters split out
markers <- FindAllMarkers(dbt_2, logfc.threshold = 0.25, min.pct=0.1,  min.diff.pct=0.1)
#saveRDS( markers, "~/Dbt_scRNA/markers_dosages_2wks_batchcorrection_rna_subclustered.rds") #saved intermediate
markers <- readRDS("~/Dbt_scRNA/markers_dosages_2wks_batchcorrection_rna_subclustered.rds")

fgsea_res <- run_fgsea(markers, hallmarks=F)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

#heatmap
p1 <- DoHeatmap(subset(dbt_2, downsample=300), features = unique(top10$gene), group.by="ident", group.color=magma(7)[1:6]) +
      scale_fill_viridis()
p1

VlnPlot(dbt_2, features=c("Progenitor_Score1", "Proliferative_Score2","Differentiated_Score3"), pt.size=0, plasma(n=8, begin=0, end=0.8)) 

#labels for clusters, see supplemental figures for some more detail on these clusters
new.cluster.ids <- c("Cycling", "Progenitor-like" ,"Progenitor-like", "Ground state","AP-1 expressing","Progenitor-like", "Differentiated", "PAX3-FOXO1-specific")
names(new.cluster.ids) <- levels(dbt_2)
dbt_2 <- RenameIdents(dbt_2, new.cluster.ids)
dbt_2$subclusters_2 <- Idents(dbt_2)

#saveRDS( dbt_2, "~/Dbt_scRNA/dosages_2wks_batchcorrection_rna_subclustered_named1.rds")

dbt_2 <- readRDS("~/Dbt_scRNA/dosages_2wks_batchcorrection_rna_subclustered_named1.rds")

markers <- readRDS("~/Dbt_scRNA/markers_dosages_2wks_batchcorrection_rna_subclustered_tempnames.rds")

top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

#heatmap
p1 <- DoHeatmap(subset(dbt_2, downsample=300), features = unique(top10$gene), group.by="ident", group.color=magma(7)[1:6]) +
      scale_fill_viridis()
p1

#run markers and gsea on final clusters
markers <- FindAllMarkers(dbt_2, logfc.threshold = 0.25, min.pct=0.2,  min.diff.pct=0.2, only.pos=T)
markers <- readRDS("~/Dbt_scRNA/markers_dosages_2wks_batchcorrection_rna_subclustered_tempnames.rds")
gsea_res <- run_fgsea(markers)
