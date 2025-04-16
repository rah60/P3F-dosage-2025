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

#quantifying the plot2DO fragment length distributions in each cluster/dosage
#Extended Data Figure 5D

# file names 

file_path <- "~/plot2DO-master/output/2D_occ_cluster_"

cluster_prefix <- c( "1_repliATAC_pulse",  "2_repliATAC_pulse",  "3_repliATAC_pulse") 

bam_files_prefix <- c("0_pulse_merged", "75_pulse_merged", "500_pulse_merged") 

# read in RData for each sample

cluster_list <- vector("list")
total_dist <- vector("list")

#get RData from plot2DO for each cluster & sample
for(f in 1:length(cluster_prefix)){
   
    dosage_list <- vector("list")

    for(i in 1:length(bam_files_prefix)){
     
       name_1 <- paste0(file_path, cluster_prefix[f],"/OCC_matrix.cluster_",cluster_prefix[f],".30_300.",bam_files_prefix[i],".RData")

        #matrix used to plot plot2DO aggregate heatmaps
       load(name_1)

        #rowSums to summarize across the entire window of plot
        temp1 <- rowSums(occMatrix)

        #fragment length distribution of whole sample
        temp3 <- lengthHist

        #normalize based on window size
        dosage_list[[i]] <- temp1/1001

        total_dist[[i]] <- lengthHist

        rm(occMatrix)

    } #end for

    cluster_list[[f]] <- dosage_list

} #end for

#get data in a format plottable by ggplot
cluster_1 <- c( cluster_list[[1]][[1]], cluster_list[[1]][[2]], cluster_list[[1]][[3]] )

cluster_2 <- c( cluster_list[[2]][[1]], cluster_list[[2]][[2]], cluster_list[[2]][[3]] )

cluster_3 <- c( cluster_list[[3]][[1]], cluster_list[[3]][[2]], cluster_list[[3]][[3]] )

frag_len <- rep(30:300, times = 3) 

dosage <- c(rep(0,times=271), rep(75,times=271), rep(500,times=271))

rows1 <- cbind(cluster_1, frag_len, dosage, rep("cluster_1",times=271) )

rows2 <- cbind(cluster_2, frag_len, dosage, rep("cluster_2",times=271) )

rows3 <- cbind(cluster_3, frag_len, dosage, rep("cluster_3",times=271) )

fragment.df <- as.data.frame(rbind(rows1, rows2, rows3))

colnames(fragment.df) <- c("coverage","frag_len","dosage","cluster")
fragment.df$coverage <- as.numeric(fragment.df$coverage)
fragment.df$frag_len <- as.numeric(fragment.df$frag_len)
fragment.df$dosage <- factor(fragment.df$dosage, levels=c("0","75","500"))
fragment.df$cluster <- factor(fragment.df$cluster, levels=c("cluster_1","cluster_2","cluster_3"), labels=c("Cluster 1","Cluster 2","Cluster 3"))

col.v1 <- viridis(8)
col.v <- col.v1[1:7] # 0, 75, 150, 250, 500, 1000, common  so 0, 75, 500 is 1, 2, 5
col_3 <- c(col.v[1],col.v[2],col.v[5])

#plot fragment length distribution at p3f sites in each cluster
p1 <- ggplot(fragment.df, aes(x=frag_len, y=coverage, col=dosage, shape=dosage))+
    geom_point(alpha=0.2)+
    geom_smooth(method="loess", span=0.1)+ #geom_smooth, line also works
    theme_classic(base_size=18)+
    coord_flip()+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0), limits=c(0,NA))+
    scale_color_manual(values=col_3)+
    labs(x="", y="Average relative coverage (%)",color="Dosage\n (ng/ml)",shape="Dosage\n (ng/ml)", title="By cluster")+
    facet_wrap(vars(cluster),scale="free")+
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# total dist for each sample

total.v <- c( total_dist[[1]], total_dist[[2]], total_dist[[3]] )
frag_len <- rep(30:300, times = 3) 
dosage <- c(rep(0,times=271), rep(75,times=271), rep(500,times=271))

total.df <- as.data.frame(cbind(total.v, frag_len, dosage))
colnames(total.df) <- c("coverage","frag_len","dosage")

total.df$coverage <- as.numeric(total.df$coverage)
total.df$frag_len <- as.numeric(total.df$frag_len)
total.df$dosage <- factor(total.df$dosage, levels=c("0","75","500"))

#plot overall fragment length distribution for sample
p2 <- ggplot(total.df, aes(x=frag_len, y=coverage, col=dosage, shape=dosage))+
    geom_point(alpha=0.2)+
    geom_smooth(method="loess", span=0.1)+ #geom_smooth, line also works
    theme_classic(base_size=18)+
    coord_flip()+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0), limits=c(0,1.3))+
    scale_color_manual(values=col_3)+
    labs(x="Fragment length (bp)", y="Coverage (%)",color="Dosage\n (ng/ml)",shape="Dosage\n (ng/ml)",title="Whole sample")+
    theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    ylim2(p1)

#combine plots
p3 <- p1 %>% insert_left(p2,width=0.28)

png("~/dosage_manuscript/figure_6/combined_frag_len_plots_pulse.png", width = 13, height = 7, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p3)
dev.off()