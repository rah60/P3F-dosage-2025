rm(list = ls())
graphics.off()

library(ggplot2)
library(plotly)
library(dplyr)
library(reshape2)
library(plotly)
library(stringr)
library(ggrepel)
library(DESeq2)
library(readr)
library(tidyr)

########### examine individual genes in SMS-CTR/SMS-CTR-iP3F, Supplemental Figure 11
#uses output from sms_ip3f_rna_Seq.R

de_result <- read_csv( "~/dosage_manuscript/csv/SMS_DE_result.csv")

temp <- de_result %>%
  dplyr::select(gene_name, SEG0585, SEG0586, SEG0577, SEG0581, SEG0578, SEG0582,SEG0579,SEG0583,SEG0580, SEG0584) %>% 
  pivot_longer(cols=starts_with("SEG0"), names_to = "sample" , values_to = "expression" ) 

key_table <- data.frame(sample = c("SEG0585", "SEG0586", "SEG0577", "SEG0581", "SEG0578", "SEG0582","SEG0579","SEG0583","SEG0580", "SEG0584"),
                        dosage = c("0","0","0","0","75","75","250","250","500","500"),
                        rep = c(1,2,1,2,1,2,1,2,1,2),
                        cell_line = c("SMS-CTR","SMS-CTR","SMS-CTR-iP3F","SMS-CTR-iP3F","SMS-CTR-iP3F","SMS-CTR-iP3F","SMS-CTR-iP3F","SMS-CTR-iP3F","SMS-CTR-iP3F","SMS-CTR-iP3F"),
                        name = c("SMS-CTR, 0","SMS-CTR, 0","SMS-CTR-iP3F, 0","SMS-CTR-iP3F, 0","SMS-CTR-iP3F, 75",
                                 "SMS-CTR-iP3F, 75","SMS-CTR-iP3F, 250","SMS-CTR-iP3F, 250","SMS-CTR-iP3F, 500","SMS-CTR-iP3F, 500"))

gene_exp_table <- left_join(temp, key_table)

gene_exp_table$name <- factor(gene_exp_table$name, levels = c("SMS-CTR, 0","SMS-CTR-iP3F, 0","SMS-CTR-iP3F, 75","SMS-CTR-iP3F, 250","SMS-CTR-iP3F, 500"))

########## cell cycle

genes_to_plot <- c("BUB1","PCNA","CCNE1")

p1 <- ggplot(subset(gene_exp_table, gene_exp_table$gene_name %in% genes_to_plot))+
  geom_boxplot(aes(x=name, y= expression, color=name))+
  geom_point(aes(x=name, y=expression, color=name))+ 
  theme_classic(base_size = 25)+
  labs(x="Cell line, doxycycline (ng/mL)",y="Expression", title="SMS-CTR")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  facet_wrap(vars(gene_name),scales="free")+
  theme(legend.position = "none",plot.title=element_text(hjust = 0.5))

genes_to_plot_string <- paste0(genes_to_plot,collapse="_")

png(paste0("~/dosage_manuscript/revision_figures/SMS_",genes_to_plot_string,".png"), width=9, height=5.3, units = "in",res=200)
print(p1)
dev.off()

######## EMT

genes_to_plot <- c("CDH2","CDH11","COL11A1")

p1 <- ggplot(subset(gene_exp_table, gene_exp_table$gene_name %in% genes_to_plot))+
  geom_boxplot(aes(x=name, y= expression, color=name))+
  geom_point(aes(x=name, y=expression, color=name))+ 
  theme_classic(base_size = 25)+
  labs(x="Cell line, doxycycline (ng/mL)",y="Expression", title="SMS-CTR")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  facet_wrap(vars(gene_name),scales="free")+
  theme(legend.position = "none",plot.title=element_text(hjust = 0.5))

genes_to_plot_string <- paste0(genes_to_plot,collapse="_")

png(paste0("~/dosage_manuscript/revision_figures/SMS_",genes_to_plot_string,".png"), width=9, height=5.3, units = "in",res=200)
print(p1)
dev.off()


######## myogenic

genes_to_plot <- c("MYOG","MYBPH","MYLPF")

p1 <- ggplot(subset(gene_exp_table, gene_exp_table$gene_name %in% genes_to_plot))+
  geom_boxplot(aes(x=name, y= expression, color=name))+
  geom_point(aes(x=name, y=expression, color=name))+ 
  theme_classic(base_size = 25)+
  labs(x="Cell line, doxycycline (ng/mL)",y="Expression", title="SMS-CTR")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  facet_wrap(vars(gene_name),scales="free")+
  theme(legend.position = "none",plot.title=element_text(hjust = 0.5))

genes_to_plot_string <- paste0(genes_to_plot,collapse="_")

png(paste0("~/dosage_manuscript/revision_figures/SMS_",genes_to_plot_string,".png"), width=9, height=5.3, units = "in",res=200)
print(p1)
dev.off()

######## FP-RMS

genes_to_plot <- c("FGFR4","FGF8","MYOD1")

p1 <- ggplot(subset(gene_exp_table, gene_exp_table$gene_name %in% genes_to_plot))+
  geom_boxplot(aes(x=name, y= expression, color=name))+
  geom_point(aes(x=name, y=expression, color=name))+ 
  theme_classic(base_size = 25)+
  labs(x="Cell line, doxycycline (ng/mL)",y="Expression", title="SMS-CTR")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  facet_wrap(vars(gene_name),scales="free")+
  theme(legend.position = "none",plot.title=element_text(hjust = 0.5))

genes_to_plot_string <- paste0(genes_to_plot,collapse="_")

png(paste0("~/dosage_manuscript/revision_figures/SMS_",genes_to_plot_string,".png"), width=9, height=5.3, units = "in",res=200)
print(p1)
dev.off()

