#clear workspace
rm(list = ls())
graphics.off()

#identify genes expressed in each dosage for RD/RD-iP3F & SMS-CTR/SMS-CTR-iP3F RNA-seq and save for use in other scripts (motif_plot.R)

library(tidyverse)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

#import bulk RNA-seq gene expression
rna_seq <- read.table("~/nextflow_rna_seq/RD_SMS_iP3F/star_salmon/salmon.merged.gene_counts.tsv",header=T) 

rna_seq2 <- rna_seq %>%
            pivot_longer(cols=starts_with("SEG0"), names_to = "sample", values_to = "counts")

sample_info <- data.frame( sample= unique(rna_seq2$sample), 
  dosage = c(0,75,250,500,0,75,250,500,0,0,0,75,250,500,0,75,250,500,0,0), 
  rep= c("rep1","rep1","rep1","rep1","rep2","rep2","rep2","rep2","rep1","rep2","rep1","rep1","rep1","rep1","rep2","rep2","rep2","rep2","rep1","rep2"), 
  cell_line = c("RD-iP3F","RD-iP3F","RD-iP3F","RD-iP3F","RD-iP3F","RD-iP3F","RD-iP3F","RD-iP3F","RD","RD",
        "SMS-iP3F","SMS-iP3F","SMS-iP3F","SMS-iP3F","SMS-iP3F","SMS-iP3F","SMS-iP3F","SMS-iP3F","SMS","SMS") )

rna_seq3 <- left_join(rna_seq2, sample_info)

data_sum <- rna_seq3 %>% #summarize to get mean expression at each dosage across replicates
  group_by(dosage, cell_line, gene_name) %>%
  summarise( mean=mean(counts))

data_sum2 <- data_sum[which(data_sum$mean > 100),] #keep only genes with count of at least 100

#make abbreviated sample for samples
data_sum2$abbrev_name <- NA

data_sum2[which(data_sum2$cell_line %in% c("RD")),5] <- "RD"
data_sum2[which(data_sum2$cell_line %in% c("SMS")),5] <- "SMS"
data_sum2[which(data_sum2$cell_line %in% c("RD-iP3F")),5] <- paste("RD", data_sum2$dosage[which(data_sum2$cell_line %in% c("RD-iP3F"))], sep="-")
data_sum2[which(data_sum2$cell_line %in% c("SMS-iP3F")),5] <- paste("SMS", data_sum2$dosage[which(data_sum2$cell_line %in% c("SMS-iP3F"))], sep="-")


RD_sum <- data_sum2[which(data_sum2$cell_line %in% c("RD","RD-iP3F")),]

SMS_sum <- data_sum2[which(data_sum2$cell_line %in% c("SMS","SMS-iP3F")),]

#save object
saveRDS(RD_sum, "~/RD_SMS_iP3F/ATAC_seq/RD_gene_expressed_at_dosages.rds")

saveRDS(SMS_sum, "~/RD_SMS_iP3F/ATAC_seq/SMS_gene_expressed_at_dosages.rds")

saveRDS(data_sum2, "~/RD_SMS_iP3F/ATAC_seq/RD_SMS_gene_expressed_at_dosages.rds")