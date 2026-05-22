rm(list = ls())
graphics.off()

#analyzes DE genes in Sroka et al data
# plots supplemental figures 13D, F

library(DESeq2)
library(dplyr)
library(tidyverse)
library(viridis)

#Dbt/MYCN/iP3F
de_result <- read.table('~/dosage_manuscript/csv/DE_result.txt', sep="\t", as.is=T, header=T, quote="\"")

#from Sroka et al 2023 supplemental data
RH4_d7 <- read_csv("~/dosage_manuscript/csv/sroka_rh4_d7_sgP3F.csv"  ) %>%
          mutate(across(c(log2FoldChange), as.numeric)) %>%
          dplyr::filter(!is.na(log2FoldChange)) %>%
          dplyr::filter(abs(log2FoldChange) >= 0.5) %>%
          dplyr::filter(padj < 0.05) 
write_csv(RH4_d7, "~/dosage_manuscript/csv/sroka_rh4_d7_selected.csv")

RH41_d7 <- read_csv("~/dosage_manuscript/csv/sroka_RH41_d7_sgP3F.csv") %>%
            mutate(across(c(log2FoldChange), as.numeric)) %>%
            dplyr::filter(!is.na(log2FoldChange)) %>%
            dplyr::filter(abs(log2FoldChange) >= 0.5) %>%
            dplyr::filter(padj < 0.05)
write_csv(RH41_d7, "~/dosage_manuscript/csv/sroka_rh41_d7_selected.csv")                                                                                                                 

#compare to up & downregulated genes in dosage data

##### RH4

##### annotate genes as up or down regulated in rh4
RH4_d7$regulation <- NA
RH4_d7[which(RH4_d7$log2FoldChange > 0),'regulation'] <- "up"
RH4_d7[which(RH4_d7$log2FoldChange < 0),'regulation'] <- "down"

#select genes in RH4 filtered set in de_result

plot_these <- which(de_result$gene_name %in% RH4_d7$Gene) 

temp <- de_result[plot_these,]

#eliminate duplicate gene name for RH4
index_temp <- which(temp$gene_name =="IDS")[1]
index_temp <- c(index_temp, which(temp$gene_name =="RBL1")[1])
temp <- temp[-index_temp,]

temp <- left_join(temp, RH4_d7, by=join_by("gene_name" == "Gene")) %>%
    mutate()

plot_these2 <- temp[,c(13,12,9,10,11,8)] %>%
  as.matrix() 

plot_these2 <- log2(plot_these2 + 1)

plot_these2 <- t(scale(t(plot_these2))) #scale gene expression

rownames(plot_these2) <- temp[,62]

up_indices <- which(rownames(plot_these2) %in% RH4_d7[which(RH4_d7$regulation == "up"),]$Gene )
down_indices <- which(rownames(plot_these2) %in% RH4_d7[which(RH4_d7$regulation == "down"),]$Gene )

temp_row <- rownames(plot_these2)
plot_these3 <- as_tibble(plot_these2)
plot_these3$gene_name <- temp_row

plot_these4 <- left_join(plot_these3, dplyr::select(RH4_d7, Gene, regulation), by=join_by(gene_name==Gene))
plot_these5 <- pivot_longer(plot_these4, cols=c(Control, X75ng.ml,X150ng.ml,X250ng.ml,X500ng.ml,X1000ng.ml),names_to = "dosage", values_to ="scaled_expression" )
plot_these5$dosage <- factor(plot_these5$dosage, levels=unique(plot_these5$dosage))

test1 <- subset(plot_these5, plot_these5$dosage == "Control")
t.test( scaled_expression ~ regulation, data= test1) #p <2.2e-16

test2 <- subset(plot_these5, plot_these5$dosage == "X75ng.ml")
t.test( scaled_expression ~ regulation, data= test2) #p < 2.2e-16

test3 <- subset(plot_these5, plot_these5$dosage == "X150ng.ml")
t.test( scaled_expression ~ regulation, data= test3) #p < 2.2e-16

test4 <- subset(plot_these5, plot_these5$dosage == "X250ng.ml")
t.test( scaled_expression ~ regulation, data= test4) #p < 2.2e-16

test5 <- subset(plot_these5, plot_these5$dosage == "X500ng.ml")
t.test( scaled_expression ~ regulation, data= test5) #p < 2.2e-16

test6 <- subset(plot_these5, plot_these5$dosage == "X1000ng.ml")
t.test( scaled_expression ~ regulation, data= test6) #p < 2.2e-16

#S13D

p1 <- ggplot(plot_these5, aes(x=dosage,y=scaled_expression,color=regulation))+
  geom_boxplot()+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  scale_x_discrete(labels=c("0","75","150","250","500","1000"))+
  theme_classic(base_size = 20)+
  labs(y="Scaled expression",x="Doxycycline (ng/mL)", 
       title="RH4 +sgP3F differentially\nexpressed genes in Dbt/MYCN/iP3F",color="Expression\n+sgP3F")+
  geom_signif(y_position = c(2.2,2.2,2.2,2.2,2.2,2.2), xmin=c(0.6,1.6,2.6,3.6,4.6,5.6), xmax=c(1.4,2.4,3.4,4.4,5.4,6.4), annotations = c("****","****","****","****","****","****"),col="black")

png( paste0("~/dosage_manuscript/revision_figures/sroka_rh4_up_down_genes.png"), width = 7, height = 5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()


############ RH41

##### annotate genes as up or down regulated in rh41
RH41_d7$regulation <- NA
RH41_d7[which(RH41_d7$log2FoldChange > 0),'regulation'] <- "up"
RH41_d7[which(RH41_d7$log2FoldChange < 0),'regulation'] <- "down"

#select genes in RH41 filtered set in de_result

plot_these <- which(de_result$gene_name %in% RH41_d7$Gene) 

temp <- de_result[plot_these,]

#eliminate duplicate gene name for RH41
index_temp <- which(temp$gene_name =="CRYBG3")[1]
index_temp <- c(index_temp, which(temp$gene_name =="RBL1")[1])
temp <- temp[-index_temp,]

temp <- left_join(temp, RH41_d7, by=join_by("gene_name" == "Gene")) %>%
  mutate()

plot_these2 <- temp[,c(13,12,9,10,11,8)] %>%
  as.matrix() 

plot_these2 <- log2(plot_these2 + 1)

plot_these2 <- t(scale(t(plot_these2))) #scale gene expression

rownames(plot_these2) <- temp[,62]

up_indices <- which(rownames(plot_these2) %in% RH41_d7[which(RH41_d7$regulation == "up"),]$Gene )
down_indices <- which(rownames(plot_these2) %in% RH41_d7[which(RH41_d7$regulation == "down"),]$Gene )

temp_row <- rownames(plot_these2)
plot_these3 <- as_tibble(plot_these2)
plot_these3$gene_name <- temp_row

plot_these4 <- left_join(plot_these3, dplyr::select(RH41_d7, Gene, regulation), by=join_by(gene_name==Gene))
plot_these5 <- pivot_longer(plot_these4, cols=c(Control, X75ng.ml,X150ng.ml,X250ng.ml,X500ng.ml,X1000ng.ml),names_to = "dosage", values_to ="scaled_expression" )
plot_these5$dosage <- factor(plot_these5$dosage, levels=unique(plot_these5$dosage))



test1 <- subset(plot_these5, plot_these5$dosage == "Control")
t.test( scaled_expression ~ regulation, data= test1) #p <2.2e-16

test2 <- subset(plot_these5, plot_these5$dosage == "X75ng.ml")
t.test( scaled_expression ~ regulation, data= test2) #p < 2.2e-16

test3 <- subset(plot_these5, plot_these5$dosage == "X150ng.ml")
t.test( scaled_expression ~ regulation, data= test3) #p < 2.2e-16

test4 <- subset(plot_these5, plot_these5$dosage == "X250ng.ml")
t.test( scaled_expression ~ regulation, data= test4) #p 3.011e-16

test5 <- subset(plot_these5, plot_these5$dosage == "X500ng.ml")
t.test( scaled_expression ~ regulation, data= test5) #p < 2.2e-16

test6 <- subset(plot_these5, plot_these5$dosage == "X1000ng.ml")
t.test( scaled_expression ~ regulation, data= test6) #p < 2.2e-16

#S13F

p2 <- ggplot(plot_these5, aes(x=dosage,y=scaled_expression,color=regulation))+
  geom_boxplot()+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  scale_x_discrete(labels=c("0","75","150","250","500","1000"))+
  theme_classic(base_size = 20)+
  labs(y="Scaled expression",x="Doxycycline (ng/mL)", 
       title="RH41 +sgP3F differentially\nexpressed genes in Dbt/MYCN/iP3F",color="Expression\n+sgP3F")+
  geom_signif(y_position = c(2.2,2.2,2.2,2.2,2.2,2.2), xmin=c(0.6,1.6,2.6,3.6,4.6,5.6), xmax=c(1.4,2.4,3.4,4.4,5.4,6.4), annotations = c("****","****","****","****","****","****"),col="black")

png( paste0("~/dosage_manuscript/revision_figures/sroka_rh41_up_down_genes.png"), width = 7, height = 5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()