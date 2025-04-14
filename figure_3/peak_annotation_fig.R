
rm(list = ls())
graphics.off()

#plots Figure 3D, Extended Data Figure 3C using output from submit_homer.sh

library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(lemon)

#set up for dataframe to store results
dosage <- c("75","150","250","500","1000")
unique_df <- vector("list")
data_names <- c("annotation","dosage", "distance_to_TSS", "gene_name", "peak_score")
plot_annotations <- data.frame()

col.v1 <- viridis(8)
col.v <- col.v1[2:6]
names(col.v) <- dosage

for(i in 1:length(dosage)){
#read in table with annotations
    unique_df[[i]] <- read.table(paste0("~/dbt_dosage_two_reps/homer_results/",dosage[i],"_2w_50bp_peaks_annotated.txt"), header = T, sep="\t",fill=T, quote = "")

    annotation.v <- sapply(strsplit(unique_df[[i]]$Annotation, "[(]"), "[",1) #process annotations in homer to make more readable/plottable

    temp <- cbind(annotation.v, rep(dosage[i],times=length(annotation.v)), #create df
             unique_df[[i]]$Distance.to.TSS,unique_df[[i]]$Gene.Name, unique_df[[i]]$Peak.Score)

    colnames(temp) <- data_names

    temp <- as.data.frame(temp)
    plot_annotations <- rbind(plot_annotations, temp) #combine temp df into overall results df

}

new.labs <- c("intergenic" , "intron", "TTS","3' UTR","TSS","non-coding","exon","5' UTR") #correct labels
names(new.labs) <- unique(plot_annotations$annotation)

#Extended Data Figure 3C

p1 <- ggplot(plot_annotations, aes(x=dosage,fill=dosage) )+
  geom_bar(show.legend = FALSE)+
  labs(y="Peak count",x="Doxycycline (ng/mL)", title="")+
  scale_x_discrete(limits = dosage)+
  scale_y_continuous(expand = c(0,0))+
   theme_classic(base_size=20)+
   scale_fill_manual(values = col.v)+
   theme(plot.title=element_text(hjust=0.5))+
   facet_wrap(vars(annotation), labeller = labeller(annotation = new.labs),nrow=2)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

png("~/dosage_manuscript/figure_2/total_peak_annotation_color_by_dosage_newscale.png", width = 12, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()


#Figure 3D

#to plot percent of each category, summarize data
new_ann <- count(plot_annotations, annotation, dosage)

#total_ann <- count(plot_annotations, annotation)

#temp test
total_ann <- count(plot_annotations, dosage)

#establish total for each category and add as a new column
new_ann$cluster_total <- NA
for(i in 1:nrow(new_ann)){

  row_num <- which(total_ann[,1] %in% new_ann[i,2]) #change to 1 for prev way
  new_ann[i,4] <- total_ann[row_num,2]

}

new_ann$percent <- new_ann$n/new_ann$cluster_total #calculate percent

p2 <- ggplot(new_ann, aes(x=dosage,y=percent, fill=dosage) )+
  geom_bar(show.legend = FALSE,stat = "identity", position = "dodge")+
  labs(y="Fraction of sites in sample",x="Doxycycline (ng/mL)", title="")+
  scale_x_discrete(limits = dosage)+
  scale_y_continuous(expand = c(0,0))+
   theme_classic(base_size=20)+
   scale_fill_manual(values = col.v)+
   theme(plot.title=element_text(hjust=0.5))+
   facet_rep_wrap(vars(annotation), labeller = labeller(annotation = new.labs),nrow=2, scales="free")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

png("~/dosage_manuscript/figure_2/total_peak_annotation_color_by_dosage_percent_3.png", width = 12, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()