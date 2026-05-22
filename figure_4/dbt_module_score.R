
library(ggplot2)
library(plotly)
library(dplyr)
library(reshape2)
library(plotly)
library(stringr)
library(ggrepel)
library(cowplot)

### calculates module scores for Dbt/MYCN/iP3F bulk RNA-seq, Supplemental Figure 8D and part of Supplemental Table 4

######## check module scores for dosages
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(ensembldb)

data <- read.csv("~/dosage_manuscript/csv/gene_count.csv", header = T, row.names = 1) 
data$ensembl <- rownames(data)

mislabled_index <- grep("^[0-9]*-[A-Z]", data$gene_name)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = data$ensembl[mislabled_index], keytype = "GENEID", columns = c("SYMBOL", "GENEID"))

for (i in 1:length(mislabled_index)) {
  data$gene_name[mislabled_index[i]] <- geneIDs1$SYMBOL[i]
}

data <- round(data[,c(25:48)])
bulk_rna <- CreateSeuratObject(counts = data)
bulk_rna$sample_name <- rownames(bulk_rna@meta.data)
bulk_rna$dosage <- rep(c("0", "75", "150", "250", "500", "1000"), each=4)

bulk_rna <- NormalizeData(bulk_rna)
#bulk_rna <- ScaleData(bulk_rna)

#set up for module scores
danielli_markers <- read_csv("~/Dbt_scRNA/danielli_2023_markers.csv")
neural_markers <- read_csv("~/neural_markers.csv")

progenitor_markers <- as.vector(danielli_markers$Progenitor_markers)
progenitor_markers  <- progenitor_markers[!is.na(progenitor_markers)]
progenitor_markers <- ensembldb::select(EnsDb.Hsapiens.v86, keys = progenitor_markers, keytype = "SYMBOL", columns = c( "GENEID"))$GENEID
  
proliferation_markers <- as.vector(danielli_markers$Proliferative_markers)
proliferation_markers  <- proliferation_markers[!is.na(proliferation_markers)]
proliferation_markers <- ensembldb::select(EnsDb.Hsapiens.v86, keys = proliferation_markers, keytype = "SYMBOL", columns = c( "GENEID"))$GENEID

differentiated_markers <- as.vector(danielli_markers$Differentiated_markers) 
differentiated_markers  <- differentiated_markers[!is.na(differentiated_markers)]
differentiated_markers <- ensembldb::select(EnsDb.Hsapiens.v86, keys = differentiated_markers, keytype = "SYMBOL", columns = c( "GENEID"))$GENEID

neural_markers2 <- as.vector(neural_markers$gene)
neural_markers2 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = neural_markers2, keytype = "SYMBOL", columns = c( "GENEID"))$GENEID

bulk_rna <- AddModuleScore(object=bulk_rna, features=list(progenitor_markers, proliferation_markers,differentiated_markers), name=c("Progenitor_Score", "Proliferative_Score","Differentiated_Score"), assay="RNA", search=T)

bulk_rna$muscle_lineage_score <- bulk_rna$Differentiated_Score3 - bulk_rna$Progenitor_Score1

neural_markers2 <- neural_markers2[which(neural_markers2 %in% rownames(bulk_rna))]
bulk_rna <- AddModuleScore(bulk_rna, features=list(neural_markers2), name=c("Neuronal_Score"),assay="RNA")

plot_scores <- bulk_rna@meta.data
plot_scores$dosage <- factor(plot_scores$dosage, levels=c("0","75","150","250","500","1000"))

prog <- aov(Progenitor_Score1 ~ dosage, data=plot_scores) %>% TukeyHSD() %>% tidy() #****
prolif <- aov(Proliferative_Score2 ~ dosage, data=plot_scores) %>% TukeyHSD() %>% tidy() #****
diff <- aov(Differentiated_Score3 ~ dosage, data=plot_scores) %>% TukeyHSD() %>% tidy() #****
neur <- aov(Neuronal_Score1 ~ dosage, data=plot_scores) %>% TukeyHSD() %>% tidy() #****

###make pairwise comp table
prog$score <- "Progenitor_Score"
prolif$score <- "Proliferative_Score"
diff$score <- "Differentiated_Score"
neur$score <- "Neuronal_Score"

save_table <- rbind(prog, prolif, diff, neur)
write_csv(save_table, "~/dosage_manuscript/figure_4/dbt_module_score_comp.csv") #for Supplemental Table 4

col.v1 <- viridis(8)
col.v <- col.v1[1:6]

p1 <- ggplot(plot_scores)+
  geom_boxplot(aes(x=dosage,y=Neuronal_Score1,color=dosage),show.legend=F)+
  geom_point(aes(x=dosage,y=Neuronal_Score1))+
  scale_y_continuous(limits = c(0,max(plot_scores$Neuronal_Score1)))+
  labs(y="Neuronal Score",x="Doxycycline (ng/mL)")+
  scale_color_manual(values = col.v)+
  theme_classic(base_size=16)

p2 <- ggplot(plot_scores)+
  geom_boxplot(aes(x=dosage,y=Progenitor_Score1, color=dosage),show.legend=F)+
  geom_point(aes(x=dosage,y=Progenitor_Score1))+
  scale_y_continuous(limits = c(0,max(plot_scores$Progenitor_Score1)))+
  labs(y="Progenitor Score",x="Doxycycline (ng/mL)")+
  scale_color_manual(values = col.v)+
  theme_classic(base_size=16)

p3 <- ggplot(plot_scores)+
  geom_boxplot(aes(x=dosage,y=Proliferative_Score2,color=dosage),show.legend=F)+
  geom_point(aes(x=dosage,y=Proliferative_Score2))+
  scale_y_continuous(limits = c(0,max(plot_scores$Proliferative_Score2)))+
  labs(y="Proliferative Score",x="Doxycycline (ng/mL)")+
  scale_color_manual(values = col.v)+
  theme_classic(base_size=16)

p4 <- ggplot(plot_scores)+
  geom_boxplot(aes(x=dosage,y=Differentiated_Score3,color=dosage),show.legend=F)+
  geom_point(aes(x=dosage,y=Differentiated_Score3))+
  scale_y_continuous(limits = c(0,max(plot_scores$Differentiated_Score3)))+
  labs(y="Differentiated Score",x="Doxycycline (ng/mL)")+
  scale_color_manual(values = col.v)+
  theme_classic(base_size=16)

png(paste0("~/dosage_manuscript/figure_4/Dbt_scores_bulk_rna.png"), width=12, height=3, units = "in",res=200)
plot_grid(p1,p2,p3,p4, nrow=1)
dev.off()
