rm(list = ls()) 
graphics.off()

library(GenomicRanges) 
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(aplot)
library(lemon)
library(enrichR)
library(viridis)
library(cowplot)

#Figure 6B, Extended Data Figures 5A & B
#using clusters established in heatmap.sh

############process datasets & join to get info together in one dataframe

#clusters from heatmap
pulse_clusters <- read.table("~/dosage_manuscript/figure_6/P3F_pulse_heatmap_110624_rpgc_kmeans3_all_sites.bed", header=T)

colnames(pulse_clusters)[1] <- "seqnames"

#dosage peaks full list from figure_3/peak_set_overlaps.R, has peak_ids used in homer 
chip_seq <- read.table("~/dbt_dosage_two_reps/all_peaks_categories.bed", sep="\t", header=F)

colnames(chip_seq) <- c("seqnames", "start","end","peak_ids", "width","strand")

#dosage peaks from figure_3/peak_set_overlaps.R
chip_seq_temp <- read.table("~/dbt_dosage_two_reps/peak_categories_GR_090624.bed", header=F)

colnames(chip_seq_temp) <- c("seqnames", "start","end", "width","strand","peak_in_75","peak_in_150","peak_in_250","peak_in_500","peak_in_1000","peak_in_0","peak_category")

#join to add peak ID
chip_seq <- left_join(chip_seq, chip_seq_temp)

#annotated dosage peaks, annotated using homer
chip_seq_2 <- read.table("~/dbt_dosage_two_reps/all_peaks_categories_annotated.txt",  sep="\t", fill=T, header=T, quote="")

colnames(chip_seq_2)[1] <- "peak_ids"

#add in annotations via peak ID
chip_seq_merged <- left_join(chip_seq, chip_seq_2) 

#drop columns that aren't needed
chip_seq_merged <- subset(chip_seq_merged, select= - c( Chr,Start,End,Strand,Peak.Score, Focus.Ratio.Region.Size)) 

#add in clusters based on repliATAC signal
pulse <- left_join(pulse_clusters, chip_seq_merged, by=join_by(seqnames, start, end))

#keep only protein-coding gene annotations
pulse <- subset(pulse, Gene.Type == "protein-coding")

#get more readable name for annotations
pulse$annotation_simplified <- sapply(str_split(pulse$Annotation, fixed("(") ),"[[", 1)


######## plots

######### Unique dosage peak plot by cluster

#Extended Data Figure 5B

#look at peaks unique to one dosage
subset_category <- c("0","75","150","250","500","1000")

hello <- pulse[which(pulse$peak_category %in% subset_category),]

hello$peak_category <- factor(hello$peak_category, levels = subset_category)

hello_2 <- count(hello, peak_category, deepTools_group)

hello_totals <- count(hello, deepTools_group)

#normalize by number of peaks in cluster
hello_2$n_norm <- NA
for(i in 1:nrow(hello_2)){
    divide_by_index <- which(hello_totals$deepTools_group == hello_2[i, 2] )
    hello_2[i,4] <- hello_2[i, 3]/hello_totals[divide_by_index,2]
}

#plot
p1 <- ggplot(data = hello_2)+
    geom_col(aes(x=peak_category, y=n_norm, fill=deepTools_group), position="dodge")+
    theme_classic(base_size = 16)+
    labs(x= "Unique peak category\n(ng/ml doxycycline)", y="Percent of peaks in cluster", fill="Peak category")+
    scale_y_continuous(expand=c(0,0))+
    scale_fill_manual(values=c("#997700","#994455","#004488"),labels=c("Cluster 1", "Cluster 2", "Cluster3"))

png("~/dosage_manuscript/figure_5/unique_peak_category_clusters.png", width = 7, height = 5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()


######## annotation plot by cluster

#Figure 6B

#summarize annotation by cluster
pulse_sum <- count(pulse, annotation_simplified, deepTools_group)

#normalize by number of peaks in cluster
pulse_totals <- count(pulse, deepTools_group)

pulse_sum$n_norm <- NA
for(i in 1:nrow(pulse_sum)){
    divide_by_index <- which(pulse_totals$deepTools_group == pulse_sum[i, 2] )
    pulse_sum[i,4] <- pulse_sum[i, 3]/pulse_totals[divide_by_index,2]
}

#plot
p2 <- ggplot(data = pulse_sum)+
    geom_col(aes(x=annotation_simplified, y= n_norm, fill=deepTools_group), position="dodge")+
    scale_y_continuous(expand=c(0,0))+
    theme_classic(base_size = 16)+
    scale_fill_manual(values=c("#997700","#994455","#004488"),labels=c("Cluster 1", "Cluster 2", "Cluster 3"))+
    labs(x= "Peak annotation", y="Percent of peaks", fill="Peak category")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=16),legend.position = c(0.2, 0.7))

png("~/dosage_manuscript/figure_6/annotations_clusters.png", width = 5, height = 4, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()


########### set up gene sets by cluster plots
#shorten dataframe for focusing on genes
enrich_genes <- pulse[,c("deepTools_group","Gene.Name")]

cluster_num <- length(unique(enrich_genes$deepTools_group))

#set up for Enrichr
#is connection working? T/F
websiteLive <- getOption("enrichR.live")
options(enrichR.base.address="https://amp.pharm.mssm.edu/Enrichr/")
#View(listEnrichrDbs())

#gene set to use
dbs_to_use <- c("MSigDB_Hallmark_2020")

#function for plotting Enrichr output
plot_gene_sets <- function(gene_sets, num_sets = 10, fill_col="black", plot_title = "",max_val=15){

    #max_val = max y-axis value/max -log10 adj p-value for plot

    subset_set <- gene_sets[which(gene_sets$Adjusted.P.value < 0.05),]

    if(nrow(subset_set) == 0){
        return(NULL)
    }
    
    subset_set$log10_adj_pval <- -log10(subset_set$Adjusted.P.value)
    subset_set$log10_pval <- -log10(subset_set$Adjusted.P.value)
    
    #select top num_sets terms & plot
    if(nrow(subset_set) >= num_sets){
        select_20 <- subset_set[order(subset_set$log10_adj_pval, decreasing=T),][1:num_sets,]
    }else{
        select_20 <- subset_set
    }
    
    order.v <- select_20[order(select_20$log10_adj_pval),1]
    
    p1 <- ggplot(select_20) +
        geom_col(aes(x=factor(Term, level=order.v), y=log10_adj_pval, fill="blank"))+
        coord_flip()+
        theme_classic(base_size=16)+
        scale_fill_manual(values=fill_col)+
        scale_y_continuous(expand = c(0,0), limits=c(0,max_val))+
        labs(x= "Gene sets" , y="-Log10 adjusted p-value", title=plot_title)+
        guides(fill="none")

    return(p1)
}

#Extended Data Figure 5A

#run Enrichr analysis
enriched_clust_1 <- enrichr(enrich_genes[which(enrich_genes$deepTools_group == "cluster_1"),2], dbs_to_use)
enriched_clust_2 <- enrichr(enrich_genes[which(enrich_genes$deepTools_group == "cluster_2"),2], dbs_to_use)
enriched_clust_3 <- enrichr(enrich_genes[which(enrich_genes$deepTools_group == "cluster_3"),2], dbs_to_use)

#set up max value for plots 
max_val <- max(max(-log10(enriched_clust_1[["MSigDB_Hallmark_2020"]]$Adjusted.P.value)) , max(-log10(enriched_clust_2[["MSigDB_Hallmark_2020"]]$Adjusted.P.value)) , max(-log10(enriched_clust_3[["MSigDB_Hallmark_2020"]]$Adjusted.P.value)) )

#plot for each cluster
clust_1 <- plot_gene_sets(enriched_clust_1[["MSigDB_Hallmark_2020"]], fill_col="#997700", plot_title = "Cluster 1", num_sets=10, max_val=max_val)

clust_2 <- plot_gene_sets(enriched_clust_2[["MSigDB_Hallmark_2020"]], fill_col="#994455", plot_title = "Cluster 2", num_sets=10, max_val=max_val)

clust_3 <- plot_gene_sets(enriched_clust_3[["MSigDB_Hallmark_2020"]], fill_col="#004488", plot_title = "Cluster 3", num_sets=10, max_val=max_val)

#combine plots and save

png( paste0("~/dosage_manuscript/figure_6/hallmarks_enrichr_clusters_plot_v2.png"), width = 12, height = 8, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(
plot_grid(clust_1, clust_2, clust_3, ncol=2)
)
dev.off()

