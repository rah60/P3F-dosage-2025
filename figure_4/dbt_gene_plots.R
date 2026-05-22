rm(list = ls())
graphics.off()

#plots expression of key genes for Supplemental Figure 8C

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
library(viridis)

de_result <- read.table('~/dosage_manuscript/csv/DE_result.txt', sep="\t", as.is=T, header=T, quote="\"")

temp <- de_result %>%
  select(gene_name, paste('H', 25:48, sep="")) %>% 
  pivot_longer(cols=starts_with("H"), names_to = "sample" , values_to = "expression" ) 

key_table <- data.frame(sample = paste('H', 25:48, sep=""),
                        dosage = c("0","0","0","0","75","75","75","75","150","150","150","150","250","250","250","250","500","500","500","500","1000","1000","1000","1000"),
                        rep = c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4),
                        name = c(rep("0 ng/mL",times=4),rep("75 ng/mL",times=4),rep("150 ng/mL",times=4),rep("250 ng/mL",times=4),rep("500 ng/mL",times=4),rep("1000 ng/mL",times=4)))

gene_exp_table <- left_join(temp, key_table)

gene_exp_table$name <- factor(gene_exp_table$name, levels = c("0 ng/mL","75 ng/mL","150 ng/mL","250 ng/mL","500 ng/mL","1000 ng/mL"))
gene_exp_table$dosage <- factor(gene_exp_table$dosage, levels = c("0","75","150","250","500","1000"))

genes_to_plot <-c("CCNE1","BUB1","CDH2","MYBPH","MYOG","CDH11","MYOD1","FGFR4")

col.v1 <- viridis(8)
col.v <- col.v1[1:7]
plot_me <- subset(gene_exp_table, gene_exp_table$gene_name %in% genes_to_plot)
plot_me$gene_name <- factor(plot_me$gene_name, levels=c("CCNE1","BUB1","CDH2","MYBPH","MYOG","CDH11","MYOD1","FGFR4"))

p1 <- ggplot(plot_me)+
  geom_boxplot(aes(x=dosage, y= expression, color=dosage))+
  geom_point(aes(x=dosage, y=expression, color=dosage))+ 
  scale_color_manual(values = col.v[1:7])+ 
  theme_classic(base_size = 20)+
  labs(x="Doxycycline dosage (ng/mL)",y="Expression")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  facet_wrap(vars(gene_name),scales="free_y",axes="all",axis.labels = "margins",nrow=2)+
  theme(legend.position = "none")

genes_to_plot_string <- paste0(genes_to_plot,collapse="_")

png(paste0("~/dosage_manuscript/figure_4/Dbt_",genes_to_plot_string,".png"), width=10, height=4.5, units = "in",res=200)
print(p1)
dev.off()