rm(list = ls())
graphics.off()

#analyzes dosage RNA-seq data, makes Extended Data Figure 4A, Figure 4A

library(DESeq2)
library(dplyr)
library(GSVA)  #need to use R 4.4.0 for this package to work
library(msigdbr)
library(tidyverse)
library(viridis)

#analysis following https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html and RNA-seq script from Meng Wang

#information about samples
sample_sheet <- read.csv("~/dosage_manuscript/csv/Sample.Information.csv")

#gene expression matrix
data <- read.csv("~/dosage_manuscript/csv/gene_count.csv", header = T, row.names = 1)  #list of genes, ensembl format name, and expression by sample, some gene info 

# remove mitocondria and GL genes, select protein-coding genes
data <- data[-grep("GL|MT", data$gene_chr), ]
data <- data[-grep("pseudogene", data$gene_description), ]
data <- data[data$gene_biotype == "protein_coding", ]

y <- data[, 1:48] #remove metadata columns

# re-factor sample_design
sample_design <- strsplit(sample_sheet$Category, "-") %>% #split category on dash & put in a dataframe
  unlist() %>%
  matrix(, nrow = 48, byrow = T) %>%
  data.frame()

rownames(sample_design) <- colnames(y)

# convert the sample design spreadsheet
sample_design <- sample_design[, -1] #get rid of first column
colnames(sample_design) <- c("Treatment", "Dosage", "Rep")

# remove excessive strings
sample_design$Treatment <- sub("DBT MYCN/", "", sample_design$Treatment) #simplify names
sample_design$Treatment <- sub("#D6", "", sample_design$Treatment)
sample_design$Treatment <- sub("#2", "", sample_design$Treatment)

# make dosage more readable
sample_design$Dosage <- sub("C", "0ng/ml", sample_design$Dosage) #simplify more names
sample_design$Dosage <- sub("Dox75 ng/ml", "75ng/ml", sample_design$Dosage)
sample_design$Dosage <- sub("Dox150 ng/ml", "150ng/ml", sample_design$Dosage)
sample_design$Dosage <- sub("Dox250 ng/ml", "250ng/ml", sample_design$Dosage)
sample_design$Dosage <- sub("Dox500 ng/ml", "500ng/ml", sample_design$Dosage)
sample_design$Dosage <- sub("Dox1000 ng/ml", "1000ng/ml", sample_design$Dosage)

# change into factor #factor as needed & convert some variable types
sample_design$Dosage <- factor(sample_design$Dosage, levels = c("0ng/ml", "75ng/ml", "150ng/ml", "250ng/ml", "500ng/ml", "1000ng/ml"))
sample_design$Rep <- as.numeric(sample_design$Rep)

dds <- DESeqDataSetFromMatrix( #set up DESeqDataSet with RNA counts & design df, specify independent variable
  countData = y, colData = sample_design,
  design = ~Dosage, ignoreRank = TRUE
)

#pca
vsd <- vst(dds, blind = F)
sampleDists <- dist(t(assay(vsd)))

pca_df <- plotPCA(vsd, intgroup = c("Dosage", "Treatment"), returnData = T)

pca_df$Dosage <- str_replace(pca_df$Dosage, "ng/ml"," ng/ml")
pca_df$Dosage <- factor(pca_df$Dosage, levels = unique(pca_df$Dosage))
pca_df$Treatment <- str_replace(pca_df$Treatment, "iP3F", "Dbt/MYCN/iP3F")
pca_df$Treatment <- str_replace(pca_df$Treatment, "iEV", "Dbt/MYCN/iEV")

# Extended Data Figure 4A

col.v1 <- viridis(8)
col.v <- col.v1[1:6]

p1 <- ggplot(pca_df)+
  geom_point(aes(x=PC1, y= PC2, color= Dosage, shape = Treatment), size=3, alpha=0.6)+
  theme_classic(base_size= 16)+
  scale_color_manual(values = col.v)+
  labs(shape = "Cell line")+
  xlab("PC1: 89.1% variance")+
  ylab("PC2: 6.3% variance")

png(paste0("~/dosage_manuscript/figure_4/pca_rnaseq_dosage.png" ), width = 8, height = 6, units = "in", res = 200, bg = "white", type = "cairo-png")
p1
dev.off()

rm(dds)

# restrict to Dbt/MYCN/iP3F only
y <- y[, which(sample_design$Treatment %in% "iP3F")] #this is count data for samples of interest only ("iP3F" only)
sample_design <- sample_design[which(sample_design$Treatment %in% "iP3F"), ] #narrow to MYCN/iP3F only, not iEV
sample_design <- sample_design[, -1] #remove treatment column since it's all the same

dds <- DESeqDataSetFromMatrix( #set up DESeqDataSet with RNA counts & design df, specify independent variable
  countData = y, colData = sample_design,
  design = ~Dosage, ignoreRank = TRUE
)
dds <- estimateSizeFactors(dds) #estimate size factors (purpose?)
se <- SummarizedExperiment(log2(counts(dds, normalized = TRUE) + 1), #summarize experiments 
                           colData = colData(dds)
)

norm_data <- assay(se)

#gene set
h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") #human hallmark gene sets: H

h_gene_sets2 <- split(x=h_gene_sets$ensembl_gene, f=h_gene_sets$gs_name)

gsvaPar <- gsvaParam(norm_data, h_gene_sets2, kcdf = "Gaussian", minSize = 15, maxSize = 500) #for some reason, had to update to new version of R to get this to work

res <- gsva(gsvaPar)

#DE
library(limma)

sample_design <- sample_design[,1:2]

test <- model.matrix(~0+sample_design$Dosage )
colnames(test) <- c("Control", "75ng/ml", "150ng/ml", "250ng/ml", "500ng/ml", "1000ng/ml")

fit <- lmFit(res, test)

fit <- eBayes(fit, trend = geneSetSizes(res))
res2 <- decideTests(fit, p.value=0.01) #picked more stringent cut-off
summary(res2)
save_me <- which(rowSums(abs(res2)) > 0)

#saveRDS(save_me, "~/dosage_manuscript/rds/dmp_dosage_sig_rows_dosages.rds")

tt <- topTable(fit,  number=Inf)

log10pval <- as.data.frame(cbind(-log10(tt$P.Value), rownames(tt) ))
log10pval$V1 <- as.numeric(log10pval$V1)

#heatmap

#Figure 4A

#rename rows & columns to something more readable
store_names <- strsplit(rownames(res), "HALLMARK_") 
store_names <- sapply(store_names,"[[",2)
store_names <- str_replace_all(store_names, "_"," ")
rownames(res) <- store_names

#set up for plotting heatmap
name_var <- paste0(sample_design$Dosage,"_",sample_design$Rep) %>%
            str_replace_all(pattern="ng/ml", replacement = " ng/ml") %>%
            str_replace_all(pattern="_", replacement = " ")

sample_design$sample_names <- name_var

colnames(res) <- sample_design$sample_names

ann_col1 <- sample_design$Dosage %>% str_replace_all(pattern="ng/ml", replacement = " ng/ml") %>% as.data.frame()

reps <- rep(c(1,2,3,4), times= length(unique(ann_col1$.)) )

ann_col1 <- as.data.frame(cbind(reps, ann_col1))

colnames(ann_col1) <- c( "Replicate", "Dosage")

ann_col1$Dosage <-  factor(ann_col1$Dosage, levels = c("0 ng/ml","75 ng/ml","150 ng/ml","250 ng/ml","500 ng/ml","1000 ng/ml"))

col_clust <- col.v
names(col_clust) <- factor(unique(ann_col1$Dosage), levels = c("0 ng/ml","75 ng/ml","150 ng/ml","250 ng/ml","500 ng/ml","1000 ng/ml"))

col_clust2 <- c(brewer.pal(n = length(unique(ann_col1$Replicate)), name ="Dark2"))
names(col_clust2) <- factor(unique(ann_col1$Replicate), levels = c(1,2,3,4))

ann_col <- list( Replicate = col_clust2, Dosage = col_clust )

library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)

col <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(50)
p2 <- pheatmap(res[save_me,], color=col,cluster_cols = F, angle_col = "315",  
               annotation_col = ann_col1, annotation_colors = ann_col, annotation_names_col=T,show_colnames = F,
          legend_breaks = c(-0.6,-0.4, -0.2, 0, 0.2,0.4,0.6), legend_labels = c("-0.6","-0.4","-0.2","0","0.2", "0.4","0.6"),main="",
          heatmap_legend_param = list(title=as.character("Scaled\nexpression"), direction="horizontal",title_position = "topcenter"),
          cellwidth = 15
)

#save for use in S4B
# saveRDS(row_order, "~/dosage_manuscript/rds/dmp_dosage_row_order_dosages.rds")

png(paste0("~/dosage_manuscript/figure_4/kmeans_6_heatmap_102224.png" ), width = 10, height = 8, units = "in", res = 200, bg = "white", type = "cairo-png")
draw(p2, heatmap_legend_side = "bottom", annotation_legend_side = "left")
dev.off()
