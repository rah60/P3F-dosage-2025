rm(list = ls())
graphics.off()

library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggrepel)

#Extended Data Figure 4B
#annotations from submit_homer_annotate.sh

#set up for dataframe
dosages <- c("Dox0_pulse","Dox75_pulse","Dox500_pulse","Dox0_chase", "Dox75_chase", "Dox500_chase")
unique_df <- vector("list")
data_names <- c("annotation","dosage", "distance_to_TSS", "gene_name", "peak_score")
plot_annotations <- data.frame()

#set up for plot
dosage <- c("0","75","500")
col.v1 <- viridis(8)
col.v <- col.v1[c(1,2,5)]
names(col.v) <- dosage

#read in & process plot annotations from homer into dataframe
for(i in 1:length(dosages)){
  #read in table with annotations
  unique_df[[i]] <- read.table(paste0("~/repliATAC/homer_results/",dosages[i],"_annotated.txt"), header = T, sep="\t",fill=T, quote = "")
  
  annotation.v <- sapply(strsplit(unique_df[[i]]$Annotation, "[(]"), "[",1)
  
  temp <- cbind(annotation.v, rep(dosages[i],times=length(annotation.v)),
                unique_df[[i]]$Distance.to.TSS,unique_df[[i]]$Gene.Name, unique_df[[i]]$Peak.Score)
  
  colnames(temp) <- data_names
  
  temp <- as.data.frame(temp)
  plot_annotations <- rbind(plot_annotations, temp)
  
} 

new.labs <- c("promoter-TSS", "intron", "intergenic", "exon", "5' UTR", "3' UTR", "non-coding", "TTS"  )
names(new.labs) <- unique(plot_annotations$annotation)

#pulse
pulse_cat <- c("Dox0_pulse","Dox75_pulse","Dox500_pulse")

pulse_df <- plot_annotations[which(plot_annotations$dosage %in% pulse_cat),]

temp1 <- str_replace_all(pulse_df$dosage, "Dox", "") 
pulse_df$dosage <- factor(str_replace_all(temp1, "_pulse", "") , levels = c("0","75","500"))

#Extended Data Figure 4B, top

#actual #
p1 <- ggplot(pulse_df, aes(x=dosage,fill=dosage) )+
  geom_bar(show.legend = FALSE)+
  labs(y="Peak count",x="Doxycycline (ng/mL)", title="")+
  scale_y_continuous(expand = c(0,0))+
  theme_classic(base_size=20)+
  scale_fill_manual(values = col.v)+
  theme(plot.title=element_text(hjust=0.5))+
  facet_wrap(vars(annotation),labeller = labeller(annotation = new.labs),nrow=2)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

png("~/dosage_manuscript/figure_5/total_peak_annotation_pulse.png", width = 12, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()

#Extended Data Figure 4B, bottom

#chase
chase_cat <- c("Dox0_chase","Dox75_chase","Dox500_chase")

chase_df <- plot_annotations[which(plot_annotations$dosage %in% chase_cat),]

temp1 <- str_replace_all(chase_df$dosage, "Dox", "") 
chase_df$dosage <- factor(str_replace_all(temp1, "_chase", "") , levels = c("0","75","500"))

#actual #
p2 <- ggplot(chase_df, aes(x=dosage,fill=dosage) )+
  geom_bar(show.legend = FALSE )+
  labs(y="Peak count",x="Doxycycline (ng/mL)", title="")+
  scale_y_continuous(expand = c(0,0))+
  theme_classic(base_size=20)+
  scale_fill_manual(values = col.v)+
  theme(plot.title=element_text(hjust=0.5))+
  facet_wrap(vars(annotation), labeller = labeller(annotation = new.labs),nrow=2)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

png("~/dosage_manuscript/figure_5/total_peak_annotation_chase.png", width = 12, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()