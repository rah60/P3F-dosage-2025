rm(list = ls())
graphics.off()

#creates Extended Data Figure 4B
#also saves expressed genes at each dosage for use in other scripts

library(DESeq2)
library(dplyr)
library(GSVA)
library(msigdbr)
library(tidyverse)
library(viridis)

#analysis following https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html and RNA-seq script from Meng Wang

col.v1 <- viridis(8)
col.v <- col.v1[1:6]

#sample info
sample_sheet <- read.csv("~/Documents/Sample.Information.csv")

#gene expression matrix
data <- read.csv("~/Documents/gene_count.csv", header = T, row.names = 1)  #list of genes, ensembl format name, and expression by sample, some gene info 

# remove mitocondria and GL genes
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

# restrict to iP3F
y <- y[, which(sample_design$Treatment %in% "iEV")] #this is count data for samples of interest only ("iP3F" only)
sample_design <- sample_design[which(sample_design$Treatment %in% "iEV"), ] #narrow to MYCN/iP3F only, not iEV
sample_design <- sample_design[, -1]

dds <- DESeqDataSetFromMatrix( #set up DESeqDataSet with RNA counts & design df, specify independent variable
  countData = y, colData = sample_design,
  design = ~Dosage, ignoreRank = TRUE
)
dds <- estimateSizeFactors(dds) #estimate size factors (purpose?)
se <- SummarizedExperiment(log2(counts(dds, normalized = TRUE) + 1), #summarize experiments (purpose?)
                           colData = colData(dds)
)

norm_data <- assay(se)

#save dataframes for other uses
# counts_long <- as.data.frame(counts) %>%
#               rownames_to_column(var="gene_name") %>%
#               pivot_longer(cols=starts_with("H"), names_to = c("condition"), values_to = c("expression"))
# 
# sample_design_long <- sample_design %>%
#                      rownames_to_column(var="condition")
# 
# data_long <- data %>%
#               select( gene_name_other = contains("gene_name")) %>%
#               rownames_to_column(var="gene_name")
# 
# counts_labeled <- left_join(counts_long, sample_design_long)
# 
# counts_labeled_2 <- left_join(counts_labeled, data_long)
# 
# #write.csv(counts_labeled_2, "~/Documents/dbt_dosage_rna_seq_mod.csv", quote=F, row.names = F)


#gene set
h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") #human hallmark gene sets: H. GO BP: C5

h_gene_sets2 <- split(x=h_gene_sets$ensembl_gene, f=h_gene_sets$gs_name)

gsvaPar <- gsvaParam(norm_data, h_gene_sets2, kcdf = "Gaussian", minSize = 15, maxSize = 500) #for some reason, had to update to new version of R to get this to work

#gsvaPar <- gsvaParam(counts, h_gene_sets2, kcdf = "Poisson") #for some reason, seems to cluster better with this...though counts are not integers, just not log2. generally trends look very similar

res <- gsva(gsvaPar)

pheatmap(res, cluster_cols = F)

#heatmap

#rename rows & columns something more readable
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

#select rows & put in same order as in Figure 4A for comparison
save_me <- readRDS("~/dosage_manuscript/rds/dmp_dosage_sig_rows_dosages.rds")
row_ordering <- readRDS("~/dosage_manuscript/rds/dmp_dosage_row_order_dosages.rds")

res <- res[save_me,]
res <- res[row_ordering,]

library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)

col <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(50)
p2 <- pheatmap(res, color=col,cluster_cols = F,cluster_rows = F, angle_col = "315",  
               annotation_col = ann_col1, annotation_colors = ann_col, annotation_names_col=T,show_colnames = F,
               legend_breaks = c(-0.6,-0.4, -0.2, 0, 0.2,0.4,0.6), legend_labels = c("-0.6","-0.4","-0.2","0","0.2", "0.4","0.6"),main="",
          heatmap_legend_param = list(title=as.character("Scaled\nexpression"), direction="horizontal",title_position = "topcenter"),
          cellwidth = 15, cellheight = 9 
)


png(paste0("~/Documents/rnaseq_iEV_heatmap_ip3f_sig_order.png" ), width =10, height =6, units = "in", res = 200, bg = "white", type = "cairo-png")
draw(p2, heatmap_legend_side = "bottom", annotation_legend_side = "left")
dev.off()


