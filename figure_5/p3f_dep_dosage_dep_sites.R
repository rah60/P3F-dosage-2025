
rm(list = ls())  
graphics.off()

library(GenomicRanges) 
library(tidyverse) 
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(aplot)
library(lemon)
library(cowplot)
library(enrichR)
library(ggsignif)
library(tidyverse)

### plots Supplemental Figure 14A

#chase
P3F_independent <- c("0c","0c_75c","0c_75c_500c","0c_500c")
P3F_dependent <- c("75c","75c_500c","500c")

chase_inde_dep <- read.table("~/dosage_manuscript/csv/peak_categories_GR_103024_chase.tsv",sep="\t")

chase_inde_dep$p3f_status <- NA

chase_inde_dep[which(chase_inde_dep$category_string %in% P3F_dependent), 10 ] <- "P3F_dependent"
chase_inde_dep[which(chase_inde_dep$category_string %in% P3F_independent), 10 ] <- "P3F_independent"

plot_chase <- chase_inde_dep %>% count(category_string, p3f_status)
plot_chase$n <- as.numeric(plot_chase$n)
plot_chase$category_string <- factor(plot_chase$category_string, levels = c("75c", "500c" ,"75c_500c" ,"0c" ,"0c_75c","0c_500c", "0c_75c_500c"  ))

#pulse

P3F_independent2 <- c("0p","0p_75p","0p_75p_500p","0p_500p")
P3F_dependent2 <- c("75p","75p_500p","500p")

pulse_inde_dep <-  read.table("~/dosage_manuscript/csv/peak_categories_GR_103024_pulse.tsv",sep="\t")

pulse_inde_dep$p3f_status <- NA

pulse_inde_dep[which(pulse_inde_dep$category_string %in% P3F_dependent2), 10 ] <- "P3F_dependent"
pulse_inde_dep[which(pulse_inde_dep$category_string %in% P3F_independent2), 10 ] <- "P3F_independent"

plot_pulse <- pulse_inde_dep %>% count(category_string, p3f_status)
plot_pulse$category_string <- factor(plot_pulse$category_string, levels = c("75p", "500p" ,"75p_500p" ,"0p" ,"0p_75p","0p_500p", "0p_75p_500p"  ))

p1 <- ggplot(data=plot_pulse)+
    geom_col(aes(x=p3f_status,fill=category_string, y=n), position="dodge")+
    scale_y_continuous(expand =  c(0,0),limits = c(0,6000))+
    scale_x_discrete(labels=c("P3F dependent","P3F independent"))+
    labs(x="Category",y="Number of repliATAC sites",fill="Doxycycline (ng/mL)", title="Nascent chromatin")+
    theme_classic(base_size=20)+
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_viridis_d(option = "H",end=0.8, labels=c("75","500","75, 500","0","0, 75","0, 500","0, 75, 500"))

p2 <- ggplot(data=plot_chase)+
    geom_col(aes(x=p3f_status,fill=category_string, y=n), position="dodge")+
    scale_y_continuous(expand =  c(0,0),limits = c(0,10000))+
    scale_x_discrete(labels=c("P3F dependent","P3F independent"))+
    labs(x="Category",y="Number of repliATAC sites",fill="Doxycycline (ng/mL)", title="Mature chromatin")+
    theme_classic(base_size=20)+
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_viridis_d(option = "H",end=0.8, labels=c("75","500","75, 500","0","0, 75","0, 500","0, 75, 500"))

 png(paste0("~/dosage_manuscript/figure_5/dosage_dep_sites.png" ), width = 14, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
  print(p1 + p2)
  dev.off()
