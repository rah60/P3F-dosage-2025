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
library(cowplot)

## analysis of RD/RD-iP3F RNA-seq data, determines module scores for RD & SMS-CTR bulk RNA-seq
# plots Supplemental Figures 10B (left), 10C, parts of Supplemental Table 4

#add gene info
gene_info <- read.table("~/dosage_manuscript/tsv/G_list.txt", header=T,sep="\t",quote="\"")

#information about samples
sample_design <- read.csv("~/dosage_manuscript/csv/sample_design_rd_sms.csv")
rownames(sample_design) <- sample_design[,4]

#gene expression matrix
data <- read.table("~/dosage_manuscript/tsv/rd_sms_salmon.merged.gene_counts.tsv", header=T)

data1 <- left_join(data, gene_info, by= join_by(gene_id == hgnc_symbol))
data2 <- data1[-grep("GL|MT", data$gene_id), ]
data2 <- data2[-which(data2$gene_name == "POLR2J4")[2],]
rownames(data2) <- data2$gene_name

y <- data2[, 3:22]

reorder_indices <- match(rownames(sample_design),colnames(y))
y <- y[,reorder_indices]
y <- round(y)


#RD cell lines only
subset_indices <- which(sample_design$Cell.line %in% c("RD","RD-iP3F") ) 
y <- y[, subset_indices]
sample_design <- sample_design[subset_indices, ] 
sample_design <- sample_design[,2:4]

sample_design$Dosage <- factor(sample_design$Dosage, levels = c("RD, 0 ng/mL", "RD-iP3F, 0 ng/mL", "RD-iP3F, 75 ng/mL", "RD-iP3F, 250 ng/mL", "RD-iP3F, 500 ng/mL"))
sample_design$Rep <- as.numeric(sample_design$Rep)

dds <- DESeqDataSetFromMatrix( #set up DESeqDataSet with RNA counts & design df, specify independent variable
  countData = y, colData = sample_design,
  design = ~Dosage, ignoreRank = TRUE
)

#pca
vsd <- vst(dds, blind = F)
sampleDists <- dist(t(assay(vsd)))

pca_df <- plotPCA(vsd, intgroup = c("Dosage"), returnData = T)
plotPCA(vsd, intgroup = c("Dosage"), returnData = F)

pca_df$Dosage <- factor(pca_df$Dosage, levels = unique(pca_df$Dosage))

# S10B

p1 <- ggplot(pca_df)+
  geom_point(aes(x=PC1, y= PC2, color= Dosage), size=3, alpha=0.6)+
  theme_classic(base_size= 16)+
  labs(color="Condition")+
  xlab("PC1: 81% variance")+
  ylab("PC2: 9% variance")

png(paste0("~/dosage_manuscript/revision_figures/pca_RD_dosage.png" ), width = 6, height = 3, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()

rm(dds)


dds <- DESeqDataSetFromMatrix(
  countData = y, colData = sample_design,
  design = ~Dosage, ignoreRank = TRUE
)
dds <- estimateSizeFactors(dds)
se <- SummarizedExperiment(log2(counts(dds, normalized = TRUE) + 1),
                           colData = colData(dds)
)

# LRT
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds <- DESeq(dds, test = "LRT", reduced = ~1)
res_ds <- results(dds)
res_ds <- data.frame(res_ds)
res_ds$gene_name <- rownames(res_ds)

# calculate average expression per gene
counts <- counts(dds, normalized = T)
counts <- data.frame(counts)
counts$gene <- rownames(counts)

# the tidyverse way of doing this is very complicated. loop will be more straightforward
avg_exp <- counts %>%
  data.frame() %>%
  melt() 

colnames(avg_exp) <- c( "gene","sample","expression"  )

avg_exp$dose <- "Control"
for (x in unique(rownames(sample_design))) {
  avg_exp$dose[which(avg_exp$sample == x)] <- as.character(sample_design$Dosage[match(x, rownames(sample_design))])
}

avg_exp <- avg_exp %>%
  group_by(dose, gene) %>%
  summarize(ae = mean(expression)) %>%
  ungroup() %>%
  dcast(formula = gene ~ dose)


x <- merge(res_ds, avg_exp, by.x = "gene_name", by.y = "gene")

# remove NA
res_ds <- x[!is.na(x$padj), ]

data$ensembl <- rownames(data)
temp <- merge(x, data, by.x = "gene_name", by.y = "gene_name")
colnames(temp)[8:12] <- c("RD-iP3F 0 ng/mL" ,  "RD-iP3F 250 ng/mL" ,"RD-iP3F 500 ng/mL", "RD-iP3F 75 ng/mL" , "RD 0 ng/mL") 

#save RD DE analysis for use in other scripts
write.csv(temp, "~/dosage_manuscript/csv/RD_SMS_DE_result.csv", quote = F, row.names = F)


############# calculate module scores for RD & SMS
library(Seurat)

sample_design <- read.csv("~/dosage_manuscript/csv/sample_design_rd_sms.csv")
rownames(sample_design) <- sample_design[,4]

#gene expression matrix
data <- read.table("~/dosage_manuscript/tsv/rd_sms_salmon.merged.gene_counts.tsv", header=T)

rownames(data) <- data$gene_id 
data <- round(data[,c(3:22)])
bulk_rna <- CreateSeuratObject(counts = data)
bulk_rna$sample_name <- rownames(bulk_rna@meta.data)

reorder_indices <- match(rownames(bulk_rna@meta.data),sample_design$sample)
sample_design <- sample_design[reorder_indices,]

bulk_rna$cell_line <- sample_design$Cell.line
bulk_rna$dosage <- sample_design$Dosage
bulk_rna$replicate <- sample_design$Rep

bulk_rna <- NormalizeData(bulk_rna)

#set up for module scores
danielli_markers <- read_csv("~/Dbt_scRNA/danielli_2023_markers.csv")
progenitor_markers <- as.vector(danielli_markers$Progenitor_markers)
progenitor_markers  <- progenitor_markers[!is.na(progenitor_markers)]
proliferation_markers <- as.vector(danielli_markers$Proliferative_markers)
proliferation_markers  <- proliferation_markers[!is.na(proliferation_markers)]
differentiated_markers <- as.vector(danielli_markers$Differentiated_markers) 
differentiated_markers  <- differentiated_markers[!is.na(differentiated_markers)]

bulk_rna <- AddModuleScore(object=bulk_rna, features=list(progenitor_markers, proliferation_markers,differentiated_markers), name=c("Progenitor_Score", "Proliferative_Score","Differentiated_Score"), assay="RNA", search=T)

bulk_rna$muscle_lineage_score <- bulk_rna$Differentiated_Score3 - bulk_rna$Progenitor_Score1

plot_scores <- bulk_rna@meta.data
plot_scores$dosage <- factor(plot_scores$dosage, levels=unique(plot_scores$dosage))

dosage_num <- str_split(plot_scores$dosage, pattern = " ")
dosage_num <- sapply(dosage_num,"[[",2)
plot_scores$dosage_num <- dosage_num
plot_scores$cell_background <- rep(c( "RD","SMS-CTR"), each=10)
plot_scores$name2 <- paste(plot_scores$cell_line, plot_scores$dosage_num, sep=", ")
plot_scores$name2 <- factor(plot_scores$name2, levels= c("RD, 0" ,"SMS-CTR, 0" ,"RD-iP3F, 0" ,"SMS-CTR-iP3F, 0","RD-iP3F, 75" ,"SMS-CTR-iP3F, 75","RD-iP3F, 250" ,"SMS-CTR-iP3F, 250","RD-iP3F, 500" ,"SMS-CTR-iP3F, 500" ) )

#RD stats
plot_scores_RD <- subset(plot_scores, cell_background == "RD")
prog <- aov(Progenitor_Score1 ~ dosage, data=plot_scores_RD ) %>% TukeyHSD() %>% tidy() #**
prolif <- aov(Proliferative_Score2 ~ dosage, data=plot_scores_RD ) %>% TukeyHSD() %>% tidy() #***
diff <- aov(Differentiated_Score3 ~ dosage, data=plot_scores_RD ) %>% TukeyHSD() %>% tidy() #****

#SMS stats
plot_scores_SMS <- subset(plot_scores, cell_background == "SMS-CTR")
prog2 <- aov(Progenitor_Score1 ~ dosage, data=plot_scores_SMS) #0.0748
prolif2 <- aov(Proliferative_Score2 ~ dosage, data=plot_scores_SMS) #0.816
diff2 <- aov(Differentiated_Score3 ~ dosage, data=plot_scores_SMS) %>% TukeyHSD() %>% tidy() #**

#make pairwise comparison table
prog$score <- "Progenitor_Score"
prog$cell_background <- "RD"
prolif$score <- "Proliferative_Score"
prolif$cell_background <- "RD"
diff$score <- "Differentiated_Score"
diff$cell_background <- "RD"
diff2$score <- "Differentiated_Score"
diff2$cell_background <- "SMS-CTR"

save_table <- rbind(prog, prolif, diff, diff2)
write_csv(save_table,"~/dosage_manuscript/revision_figures/module_score_comp.csv") #for Supplemental Table 4

#S10C

p1 <- ggplot(plot_scores)+
  geom_boxplot(aes(x=name2,y=Progenitor_Score1, color=cell_line))+
  geom_point(aes(x=name2,y=Progenitor_Score1, color=cell_line))+
  theme_classic(base_size=16)+
  labs(y="Progenitor Score", x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  facet_wrap(vars(cell_background),scales="free")+
  theme(legend.position = "none")

p2 <- ggplot(plot_scores)+
  geom_boxplot(aes(x=name2,y=Differentiated_Score3, color=cell_line))+
  geom_point(aes(x=name2,y=Differentiated_Score3, color=cell_line))+
  theme_classic(base_size=16)+
  labs(y="Differentiated Score", x="Cell line, doxycycline (ng/mL)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  facet_wrap(vars(cell_background),scales="free")+
  theme(legend.position = "none")

p3 <- ggplot(plot_scores)+
  geom_boxplot(aes(x=name2,y=Proliferative_Score2, color=cell_line))+
  geom_point(aes(x=name2,y=Proliferative_Score2, color=cell_line))+
  theme_classic(base_size=16)+
  labs(y="Proliferative Score", x="Cell line, doxycycline (ng/mL)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  facet_wrap(vars(cell_background),scales="free")+
  theme(legend.position = "none")

p5 <- plot_grid(p1, p2, p3, ncol=2)

png(paste0("~/dosage_manuscript/revision_figures/RD_module_scores.png"), width=12, height=7, units = "in",res=200,bg = "transparent", type = "cairo-png")
print(p5)
dev.off()




