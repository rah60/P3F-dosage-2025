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

# analyze SMS-CTR/SMS-CTR-iP3F bulk RNA-seq data
# plot Supplemental Figure 10B (right)

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

#SMS cell lines only
subset_indices <- which(sample_design$Cell.line %in% c("SMS-CTR","SMS-CTR-iP3F") ) 
y <- y[, subset_indices]
sample_design <- sample_design[subset_indices, ] 
sample_design <- sample_design[,2:4]

sample_design$Dosage <- factor(sample_design$Dosage, levels = unique(sample_design$Dosage))
sample_design$Rep <- as.numeric(sample_design$Rep)

dds <- DESeqDataSetFromMatrix( #set up DESeqDataSet with RNA counts & design df, specify independent variable
  countData = y, colData = sample_design,
  design = ~Dosage, ignoreRank = TRUE
)

#S10B
#pca
vsd <- vst(dds, blind = F)
sampleDists <- dist(t(assay(vsd)))

pca_df <- plotPCA(vsd, intgroup = c("Dosage"), returnData = T)
plotPCA(vsd, intgroup = c("Dosage"), returnData = F)

pca_df$Dosage <- factor(pca_df$Dosage, levels = unique(pca_df$Dosage))


p1 <- ggplot(pca_df)+
  geom_point(aes(x=PC1, y= PC2, color= Dosage), size=3, alpha=0.6)+
  theme_classic(base_size= 16)+
  labs(color="Condition")+
  xlab("PC1: 72% variance")+
  ylab("PC2: 12% variance")

png(paste0("~/dosage_manuscript/revision_figures/pca_SMS_dosage.png" ), width = 6, height = 3, units = "in", res = 200, bg = "transparent", type = "cairo-png")
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
res_ds$gene_name<- rownames(res_ds)

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
colnames(temp)[8:12] <- c("SMS-CTR-iP3F 0 ng/mL",   "SMS-CTR-iP3F 250 ng/mL", "SMS-CTR-iP3F 500 ng/mL", "SMS-CTR-iP3F 75 ng/mL" , "SMS-CTR 0 ng/mL" ) 

#save SMS DE analysis for use in downstream scripts
write.csv(temp, "~/dosage_manuscript/csv/SMS_DE_result.csv", quote = F, row.names = F)

