#clear workspace
rm(list = ls())
graphics.off()

#identify genes expressed in each dosage and save for use in other scripts

library(tidyverse)

#import bulk RNA-seq gene expression
rna_seq <- read_csv("~/dosage_manuscript/csv/dbt_dosage_rna_seq_mod.csv")

data_sum <- rna_seq %>% #summarize to get mean expression at each dosage across replicates
  group_by(Dosage, gene_name) %>%
  summarise( mean=mean(expression))

data_sum_2 <- left_join(data_sum, unique(rna_seq[,c(1,6)]), relationship="many-to-one") #combine back into dataframe

data_sum_2 <- data_sum_2[which(data_sum_2$mean > 100),] #keep only genes with count of at least 100

new_dosages <-  str_split(data_sum_2$Dosage, pattern="n")
new_dosages_2 <- sapply(new_dosages,"[[",1) #get only the number from dosage for ease of use downstream

data_sum_2$Dosage <- new_dosages_2

#save object
saveRDS(data_sum_2, "~/dosage_manuscript/rds/gene_expressed_at_dosages.rds")

