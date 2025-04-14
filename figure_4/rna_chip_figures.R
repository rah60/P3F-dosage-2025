rm(list = ls())
graphics.off()

#plot Figure 4B, also establishes clusters used in extended data figure 3

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(fgsea)
library(msigdbr)
library(cowplot)
library(enrichR)
 
#get overlap in genes near P3F ChIP-seq peaks (start w/ all genes, w/ annotation of their categories) and RNA-seq genes that are DE in any way
#make sure to have normalized avg signal for these genes across dosages
rna_seq <- read.table("~/dosage_manuscript/csv/DE_result.txt", sep="\t", as.is=T, header=T, quote="\"")
rna_seq <- rna_seq[which(rna_seq$padj < 0.05),]

#need to annotate this table with nearest gene for each peak--take peaks to bed file, annotate w/ homer, add back to this dataframe
chip_seq <- read.table("~/dosage_manuscript/csv/peak_categories_GR_090624.tsv", header=T, sep="\t")

peak_ids <- c()
for(i in 1:nrow(chip_seq)){
    new_name <- paste0("peak_",i)
    peak_ids <- c(peak_ids , new_name)
}

chip_seq$peak_ids <- peak_ids

#for homer annotation
new_chip_seq <- cbind(chip_seq[,c(1:3)], peak_ids,chip_seq[,c(4,5)] )
write.table(new_chip_seq, "~/dbt_dosage_two_reps/all_peaks_categories.bed", sep="\t", quote=F, row.names=F, col.names=F, eol='\r\n') #run this thru homer

#annotated by homer, in submit_homer_annotate.sh
chip_seq_2 <- read.table("~/dbt_dosage_two_reps/all_peaks_categories_annotated.txt",  sep="\t", fill=T, header=T, quote="")

colnames(chip_seq_2)[1] <- "peak_ids"

chip_seq_merged <- merge(chip_seq, chip_seq_2) 

chip_seq_merged <- subset(chip_seq_merged, select= - c( Chr,Start,End,Strand,Peak.Score, Focus.Ratio.Region.Size))

chip_seq_merged <- subset(chip_seq_merged, Gene.Type == "protein-coding")

#narrow to 50kb from TSS
tss_subset <- which(abs(chip_seq_merged$Distance.to.TSS) <= 50000)
chip_seq_merged <- chip_seq_merged[tss_subset,]

#select overlapping genes
get_these <- which(chip_seq_merged$Nearest.Ensembl %in% rna_seq$ensembl) #which chip-seq peaks nearest genes are in the rna-seq

get_these_subset <- chip_seq_merged[get_these,] #subset those

exp_genes <- unique(get_these_subset$Nearest.Ensembl) #unique the subset of genes

get_these2 <- which(rna_seq$ensembl %in% exp_genes) #figure out indices of genes subset in rna-seq

plot_these <- rna_seq[get_these2,] #get rna-seq subset

#make a heatmap w/ the rna-seq signal, kmeans cluster

#matrix with rownames & colnames
plot_these2 <- plot_these[,c(13,12,9,10,11,8)] %>%
    as.matrix() 
plot_these2 <- log2(plot_these2 + 1)

plot_these2 <- t(scale(t(plot_these2))) #scale gene expression

rownames(plot_these2) <- plot_these$gene_name

breakList <- seq(-2, 2, by=0.4)
col <- rev(colorRampPalette((brewer.pal(n = 9, name ="RdBu")))(length(breakList) - 1))

p1 <- pheatmap::pheatmap(plot_these2, color=col,cluster_cols=F, cluster_rows=T, show_rownames=F, #k-means clusters heatmap
    breaks = breakList,kmeans_k=6 )[["kmeans"]]

#saveRDS(p1, "~/dbt_dosage_two_reps/rna_chip_kmeans_6_tss50kb_cutoff.rds") #save clusters info
p1 <- readRDS("~/dbt_dosage_two_reps/clusters/rna_chip_kmeans_6_tss50kb_cutoff.rds")

#pheatmap::pheatmap(p1[[2]])

#get cluster info into dataframe with genes & expression data (plot_these, pre-scaling)
genes_clusters <- p1[["cluster"]]
genes_clusters <- genes_clusters %>%
                    as.data.frame() %>%
                    rownames_to_column()

colnames(genes_clusters) <- c("gene_name", "cluster")

clustered_df <- left_join(plot_these, genes_clusters)

#normalize and scale gene expression data
clustered_df_2 <- clustered_df[,c(13,12,9,10,11,8)] %>%
    as.matrix() 
clustered_df_2 <- log2(clustered_df_2 + 1) %>%
                t() %>%
                scale() %>%
                t()

#prepare for plotting heatmap
clustered_df_3 <- clustered_df_2[order(clustered_df$cluster),]

clustered_df_4 <- clustered_df[order(clustered_df$cluster),]

rownames(clustered_df_3) <- clustered_df_4$gene_name

genes_clusters1 <- dplyr::select(clustered_df_4, c(gene_name, cluster)) %>%
                    remove_rownames() %>%
                    column_to_rownames( var="gene_name")


colnames(genes_clusters1) <- c("Cluster")
breakList <- seq(-2, 2, by=0.4)
col <- rev(colorRampPalette((brewer.pal(n = 9, name ="RdBu")))(length(breakList) - 1))

col_names <- c("0 ng/ml","75 ng/ml","150 ng/ml","250 ng/ml","500 ng/ml","1000 ng/ml")
genes_clusters1$Cluster <- as.character(genes_clusters1$Cluster)

col_clust <- c(brewer.pal(n = length(unique(genes_clusters1$Cluster)), name ="Paired"))
names(col_clust) <- c("1","2","3","4","5","6")

ann_col <- list( Cluster = col_clust )

library(ComplexHeatmap)

#plot Figure 4B heatmap
heatmap_var <- pheatmap(clustered_df_3, color=col,cluster_cols=F, cluster_rows=F, show_rownames=F,
       breaks = breakList, annotation_row = genes_clusters1, annotation_colors = ann_col, annotation_names_row=F,
       labels_col = as.character(col_names), angle_col = "315", heatmap_legend_param = list(title=as.character("Scaled\nexpression"), legend_height = unit(4,"cm")) , 
       fontsize = 14, row_split = genes_clusters1$Cluster,
       cellwidth=40 )

png( paste0("~/dosage_manuscript/figure_4/heatmap_kmeans6_subset50kb_cutoff.png"), width = 5, height = 10, units = "in", res = 200, bg = "transparent", type = "cairo-png")
draw(heatmap_var)
dev.off()

#Figure 4B gene set enrichment
#examine enriched gene sets in each cluster using Enrichr

enrich_genes <- clustered_df[,c("cluster","gene_name")]

cluster_num <- length(unique(enrich_genes$cluster))

#is connection working? T/F
websiteLive <- getOption("enrichR.live")
#View(listEnrichrDbs())

dbs_to_use <- c("MSigDB_Hallmark_2020")

#produces the tables needed to make plots, just select which gene set is of interest

plot_gene_sets <- function(gene_sets, num_sets = 10, fill_col="black", plot_title = ""){

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
        scale_y_continuous(expand = c(0,0))+
        labs(x= "Gene sets" , y="-Log10 adjusted p-value", title=plot_title)+
        guides(fill="none")

    return(p1)
}

#establish plot for each cluster after running Enrichr for genes in each cluster
col <- brewer.pal(n = length(unique(genes_clusters1$Cluster)), name ="Paired")

enriched_clust_1 <- enrichr(enrich_genes[which(enrich_genes$cluster == 1),2], dbs_to_use)

clust_1 <- plot_gene_sets(enriched_clust_1[["MSigDB_Hallmark_2020"]], fill_col=col[1], plot_title = "Cluster 1")

enriched_clust_2 <- enrichr(enrich_genes[which(enrich_genes$cluster == 2),2], dbs_to_use)

clust_2 <- plot_gene_sets(enriched_clust_2[["MSigDB_Hallmark_2020"]], fill_col=col[2], plot_title = "Cluster 2")

enriched_clust_3 <- enrichr(enrich_genes[which(enrich_genes$cluster == 3),2], dbs_to_use)

clust_3 <- plot_gene_sets(enriched_clust_3[["MSigDB_Hallmark_2020"]], fill_col=col[3], plot_title = "Cluster 3")

enriched_clust_4 <- enrichr(enrich_genes[which(enrich_genes$cluster == 4),2], dbs_to_use)

clust_4 <- plot_gene_sets(enriched_clust_4[["MSigDB_Hallmark_2020"]], fill_col=col[4], plot_title = "Cluster 4")

enriched_clust_5 <- enrichr(enrich_genes[which(enrich_genes$cluster == 5),2], dbs_to_use)

clust_5 <- plot_gene_sets(enriched_clust_5[["MSigDB_Hallmark_2020"]], fill_col=col[5], plot_title = "Cluster 5")

enriched_clust_6 <- enrichr(enrich_genes[which(enrich_genes$cluster == 6),2], dbs_to_use)

clust_6 <- plot_gene_sets(enriched_clust_6[["MSigDB_Hallmark_2020"]],fill_col=col[6], plot_title = "Cluster 6")

#plot in two pngs so they can be easily arranged around heatmap
png( paste0("~/dosage_manuscript/figure_4/hallmarks_plot_1_kmeans6_subset50kb_cutoff.png"), width = 6, height = 12, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(
plot_grid(clust_2, clust_3, clust_4, ncol=1)
)
dev.off()

png( paste0("~/dosage_manuscript/figure_4/hallmarks_plot_2_kmeans6_subset50kb_cutoff.png"), width = 6, height = 12, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(
plot_grid(clust_1, clust_5, clust_6, ncol=1)
)
dev.off()


