rm(list = ls())
graphics.off()

#analyzes DE genes from Zhang et al paper between 72 hr dTAG and DMSO
# plots Supplemental Figures 13A & B

library(DESeq2)
library(dplyr)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
library(fgsea)
library(msigdbr)
library(reshape2)
library(ggsignif)

###### first identify DE genes/analyze w/ DE-seq
timepoints <- c(0,72)

gene_info <- read.table("~/dosage_manuscript/tsv/G_list.txt", header=T,sep="\t",quote="\"")

sample_sheet <- read.csv("~/dosage_manuscript/csv/sample_design_revisions_Zhang.csv")
sample_sheet <- sample_sheet[which(sample_sheet$technique=="RNA-seq"),]
sample_sheet <- sample_sheet[which(sample_sheet$Cell.line=="RH30_dTAG"),]
sample_sheet <- sample_sheet[which(sample_sheet$Time %in% timepoints),]
rownames(sample_sheet) <- sample_sheet[,4]

#from GSE183297 & then nf-core RNA-seq pipeline analysis
data <- read.table("~/dosage_manuscript/tsv/zhang_2022_salmon.merged.gene_counts.tsv", header=T)

data1 <- left_join(data, gene_info, by= join_by(gene_id == hgnc_symbol))
data2 <- data1[-grep("GL|MT", data1$gene_id), ]
data2 <- data2[-which(data2$gene_name == "POLR2J4")[2],]
rownames(data2) <- data2$gene_name

y <- data2[, 3:ncol(data2)]

reorder_indices <- match(rownames(sample_sheet),colnames(y))
y <- y[,reorder_indices]
y <- round(y)

sample_design <- sample_sheet[,2:4]

sample_design$Rep <- as.numeric(sample_design$Rep)
sample_design$Time <- factor(sample_design$Time, levels=unique(sample_design$Time) )

dds <- DESeqDataSetFromMatrix(
  countData = y, colData = sample_design,
  design = ~Time, ignoreRank = TRUE
)
dds <- estimateSizeFactors(dds)
se <- SummarizedExperiment(log2(counts(dds, normalized = TRUE) + 1),
                           colData = colData(dds)
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds <- DESeq(dds, test = "Wald")
res_ds <- results(dds)
res_ds <- data.frame(res_ds)
res_ds$gene_name <- rownames(res_ds)

counts <- counts(dds, normalized = T)
counts <- data.frame(counts)
counts$gene <- rownames(counts)

avg_exp <- counts %>%
  data.frame() %>%
  melt() 

colnames(avg_exp) <- c( "gene","sample","expression"  )

avg_exp$Time <- "0"
for (x in unique(rownames(sample_design))) {
  avg_exp$Time[which(avg_exp$sample == x)] <- as.character(sample_design$Time[match(x, rownames(sample_design))])
}

avg_exp <- avg_exp %>%
  group_by(Time, gene) %>%
  summarize(ae = mean(expression)) %>%
  ungroup() %>%
  dcast(formula = gene ~ Time)


x <- merge(res_ds, avg_exp, by.x = "gene_name", by.y = "gene")

# remove NA
x <- x[!is.na(x$padj), ]

data2$gene_name <- rownames(data2)
temp <- merge(x, data2, by.x = "gene_name", by.y = "gene_name")

write_csv(temp, "~/dosage_manuscript/csv/zhang_RNA_RH30_dTAG_72_DE_result.csv")

de_result <- read_csv("~/dosage_manuscript/csv/zhang_RNA_RH30_dTAG_72_DE_result.csv")

###### run GSEA for S13A
h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") #human hallmark gene sets: H. GO BP: C5
h_gene_sets2 <- split(x=h_gene_sets$gene_symbol, f=h_gene_sets$gs_name) 

markers.v <- de_result$log2FoldChange #avglogFC

names(markers.v) <- de_result$gene_name

fgsea.l <- fgsea(pathways = h_gene_sets2, stats=markers.v , maxSize=500, minSize=5, nPerm=1000) %>% #or both & run correct FGSEA settings
  arrange(-NES)

plot_pathways <- fgsea.l[which(fgsea.l$padj < 0.05),] %>%
  arrange(NES)

store_names <- strsplit(plot_pathways$pathway, "HALLMARK_") 
store_names <- sapply(store_names,"[[",2)
store_names <- str_replace_all(store_names, "_"," ")
plot_pathways$pathway <- store_names

plot_pathways$pathway <- factor(plot_pathways$pathway, levels=plot_pathways$pathway)

p1 <- ggplot(plot_pathways)+
  geom_count(aes(y=pathway,x=NES,color=padj))+
  labs(y="",x="Normalized Enrichment Score",title = "Enriched Hallmark Gene Sets\n 72 hrs P3F dTAG",color="Adjusted\np-value")+
  scale_color_viridis_c(begin=0.8, end=0, option="turbo")+
  theme_classic(base_size=12)+
  guides(size="none")

png( paste0("~/dosage_manuscript/revision_figures/zhang_72_gsea.png"), width = 6, height = 3, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()

########## aggregate gene expression for S13B
de_result2 <- read.table('~/dosage_manuscript/csv/DE_result.txt', sep="\t", as.is=T, header=T, quote="\"")

RH30_DE <- de_result %>% dplyr::filter(abs(log2FoldChange) >= 0.5) %>%
  dplyr::filter(padj < 0.05) 
RH30_DE$regulation <- NA
RH30_DE$regulation[which(RH30_DE$log2FoldChange >0)] <- "up"
RH30_DE$regulation[which(RH30_DE$log2FoldChange <0)] <- "down"

plot_these <- which(de_result2$gene_name %in% RH30_DE$gene_name)

temp <- de_result2[plot_these,]

temp <- left_join(temp, RH30_DE, by=join_by("gene_name" == "gene_name")) 

plot_these2 <- temp[,c(13,12,9,10,11,8)] %>%
  as.matrix() 

plot_these2 <- log2(plot_these2 + 1)

plot_these2 <- t(scale(t(plot_these2))) #scale gene expression

rownames(plot_these2) <- temp[,62]

up_indices <- which(rownames(plot_these2) %in% RH30_DE[which(RH30_DE$regulation == "up"),]$gene_name )
down_indices <- which(rownames(plot_these2) %in% RH30_DE[which(RH30_DE$regulation == "down"),]$gene_name )

temp_row <- rownames(plot_these2)
plot_these3 <- as_tibble(plot_these2)
plot_these3$gene_name <- temp_row

plot_these4 <- left_join(plot_these3, dplyr::select(RH30_DE, gene_name, regulation), by=join_by(gene_name==gene_name))
plot_these5 <- pivot_longer(plot_these4, cols=c(Control, X75ng.ml,X150ng.ml,X250ng.ml,X500ng.ml,X1000ng.ml),names_to = "dosage", values_to ="scaled_expression" )
plot_these5$dosage <- factor(plot_these5$dosage, levels=unique(plot_these5$dosage))

test1 <- subset(plot_these5, plot_these5$dosage == "Control")
t.test( scaled_expression ~ regulation, data= test1) #p 0.7118

test2 <- subset(plot_these5, plot_these5$dosage == "X75ng.ml")
t.test( scaled_expression ~ regulation, data= test2) #p 4.94e-06

test3 <- subset(plot_these5, plot_these5$dosage == "X150ng.ml")
t.test( scaled_expression ~ regulation, data= test3) #p 2.011e-10

test4 <- subset(plot_these5, plot_these5$dosage == "X250ng.ml")
t.test( scaled_expression ~ regulation, data= test4) #p 0.7138

test5 <- subset(plot_these5, plot_these5$dosage == "X500ng.ml")
t.test( scaled_expression ~ regulation, data= test5) #p 0.004145

test6 <- subset(plot_these5, plot_these5$dosage == "X1000ng.ml")
t.test( scaled_expression ~ regulation, data= test6) #p 0.05613

p2 <- ggplot(plot_these5, aes(x=dosage,y=scaled_expression,color=regulation))+ 
  geom_boxplot()+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  scale_x_discrete(labels=c("0","75","150","250","500","1000"))+
  theme_classic(base_size = 20)+
  labs(y="Scaled expression",x="Doxycycline (ng/mL)", 
       title="RH30 +P3F-dTAG differentially\nexpressed genes in Dbt/MYCN/iP3F",color="Expression\n+dTAG-47")+
  geom_signif(y_position = c(1.4,1.8,1), xmin=c(1.6,2.6,4.6), xmax=c(2.4,3.4,5.4), annotations = c("****","****","**"),col="black")

png( paste0("~/dosage_manuscript/revision_figures/zhang_72_up_down_genes.png"), width = 7, height = 5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()

