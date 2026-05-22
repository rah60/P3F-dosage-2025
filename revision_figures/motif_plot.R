#clear workspace
rm(list = ls())
graphics.off()

#compare motifs in the homer output file knownResults.txt via p-values/motif enrichment
# uses motif analysis from submit_homer_RD_chip.sh and submit_homer_SMS_chip.sh
#Supplemental Figure 12C, parts of Supplemental Table 2

library(gplots)
library(ggplot2)
library(ggrepel)
library(viridis)
library(tidyverse)
library(lemon)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(RColorBrewer)

#list prefixes for categories of peaks for which to plot motifs--within same vector if you want them plotted together, each item in list a separate plot
plots_list <- list(
          c("RD","RD-0","RD-75","RD-250","RD-500"),
          c("SMS","SMS-0","SMS-75","SMS-250","SMS-500")
          )

#number of categories in each plot
new_names <- c("RD", "SMS")

#plot name prefix
file_var <- c("RD_chip_overlap","SMS_chip_overlap")

path_prefix <- "~/RD_SMS_iP3F/ATAC_seq/peak_categories_chip"

#functions

#make function to take dosage & return vector of gene names
#### then I can split strings on "_" to get dosages, and get as many gene vectors as needed for given subset

#read in motif data table
read_motif_data <- function(x, path_prefix = "~/RD_SMS_iP3F/ATAC_seq/peak_categories_chip/"){ #x is data name prefix
    results_path <- paste(path_prefix,x,"/knownResults.txt",sep="") 
    results <- read.table(results_path, header = T, sep="\t", comment.char = "") 
    results
  }

#get motif types from dataframe, unlisted
get_homer_motif_types_df <- function(x){ #x is dataframe
    temp <- x$Motif.Name #modified to not select top 50 based on prior filtering steps
    temp1 <- strsplit(temp, "/")
    temp2 <- strsplit(sapply(temp1,"[[",1),"[(]")
    temp3 <- strsplit(sapply(temp2, "[",2),"[)]")
    unlist(sapply(temp3, "[[",1))
    }

#get vector of motif types with names that are the actual factors
get_list_motif_factor_names_df <- function(x){ #x is dataframe
    temp <- x$Motif.Name
    temp1 <- strsplit(temp, "/")
    temp2 <- strsplit(sapply(temp1,"[[",1),"[(]")
    temp3 <- strsplit(sapply(temp2, "[",1),"[)]")
    temp4 <- unlist(sapply(temp3, "[[",1))
    temp5 <- strsplit(temp4, "\\+")
    part_1 <- unlist(sapply(temp5, "[[",1))

    temp6 <- strsplit(sapply(temp2, "[",2),"[)]")
    part_2 <- unlist(sapply(temp6, "[[",1))

    names(part_2) <- part_1
    part_2   
    }

#get motif types from vector   
get_homer_motif_types <- function(x){ #x is vector of homer motif names
    temp1 <- strsplit(x, "/")
    temp2 <- strsplit(sapply(temp1,"[[",1),"[(]")
    temp3 <- strsplit(sapply(temp2, "[",2),"[)]")
    sapply(temp3, "[[",1)
    }

#get motif names from vector
get_homer_motif_names <- function(x){ #x is vector of homer motif names
    temp1 <- strsplit(x, "/")
    temp2 <- strsplit(sapply(temp1,"[[",1),"[(]")
    sapply(temp2, "[",1)
    }

#import list of expressed genes from bulk RNA-seq generated in make_expressed_gene_list.R
#and get genes expressed at a dosage  
get_dosage_gene_vector <- function(x, genes_rds = "~/RD_SMS_iP3F/ATAC_seq/RD_SMS_gene_expressed_at_dosages.rds"){ #x is one dosage
  genes_df <- readRDS(genes_rds)
  subset_x <- genes_df[which(genes_df$abbrev_name == x),]
  subset_x$gene_name
}

for(t in 1:length(plots_list)){
  
  #vector of categories for plot
  plots.v <- plots_list[[t]]

  #read in motif data for each category/dosage
  motifs_list <- lapply(plots.v, read_motif_data)

#get all potential motif names across samples
  all_factor_names_motif_type <- lapply(motifs_list, get_list_motif_factor_names_df) %>%
                      unlist() 
 
if("ZNF652" %in% names(all_factor_names_motif_type)){
 all_factor_names_motif_type[which(names(all_factor_names_motif_type) == "ZNF652")] <- "Zf"
}

 #unique factor names across all samples                     
  all_factor_names <- unique(names(all_factor_names_motif_type))

if("NF1" %in% all_factor_names){
  all_factor_names[which(all_factor_names == "NF1")] <- "NFIA" } #fixing homer naming issue

names(all_factor_names) <- toupper(all_factor_names) #to avoid case issues
P3F_index <- which(all_factor_names == "PAX3:FKHR-fusion") 

#get dosages
temp1 <- strsplit(plots.v, "_")
dosages.v <- unique(unlist(temp1))

#compare to the list of expressed genes in RNA seq
all_exp_genes <- lapply(dosages.v, get_dosage_gene_vector) %>%
                      unlist() %>%
                      unique()

#get gene symbol of expressed genes and factor names overlap
get_these <- all_factor_names[which(all_factor_names %in% all_exp_genes)]

#narrow based on expressed factors
keep_names <- c(get_these, all_factor_names[P3F_index]) #make sure P3F remains in list, b/c no ensembl id for it

#narrow motif families based on expressed factors, ultimately I am just filtering out motifs for which there are no expressed factors in that family (as opposed to filtering expressed factors themselves)
keep_motifs <- unique(all_factor_names_motif_type[keep_names])

#for any motif names that don't have a match in the expressed gene set, drop that from the dataset 
 for(w in 1:length(motifs_list)){

  keep_index <- which(get_list_motif_factor_names_df(motifs_list[[w]]) %in% keep_motifs)
  motifs_list[[w]] <- motifs_list[[w]][keep_index,]
  motifs_list[[w]] <- motifs_list[[w]][1:50,]

 }

#get all potential motif types across samples

  all_motif_types <- lapply(motifs_list, get_homer_motif_types_df) %>%
                      unlist() %>%
                      unique()

if("TEA" %in% all_motif_types){ #using TEAD for motifs that Homer labels as TEA and TEAD
  all_motif_types <- all_motif_types[-which(all_motif_types == "TEA")]}

if("Zf" %in% all_motif_types){
  all_motif_types <- all_motif_types[-which(all_motif_types == "Zf")] #using more readable name
  all_motif_types <- c(all_motif_types, "Zinc finger")} #to make this clearer
  
  col.v1 <- magma(length(all_motif_types)+1)
  col.v <- col.v1[1:(length(all_motif_types))]
  names(col.v) <- all_motif_types

  results_df <- data.frame()

    for(i in 1:length(motifs_list)){ 

    #get motif types and names in dataframe, along with pval, logpval
    
    subset_motifs <- motifs_list[[i]][,c(1,3,4)] 
    subset_motifs$name_only <- get_homer_motif_names(subset_motifs$Motif.Name)
    subset_motifs$type_only <- get_homer_motif_types(subset_motifs$Motif.Name)
    subset_motifs$minus_Log.P.value <- as.numeric(subset_motifs$Log.P.value)*-1

    #consolidate motif categories
    subset_motifs[which(subset_motifs$type_only == "TEA"),5] <- "TEAD"
    subset_motifs[which(subset_motifs$type_only == "Zf"),5] <- "Zinc finger"

    subset_motifs$type_only <- factor(subset_motifs$type_only, levels=all_motif_types)

  #add extra rows to make sure to plot all categories even if there are no motifs in category for that sample.
    add_motifs <- all_motif_types[which(!(all_motif_types %in% subset_motifs$type_only))]
    
    if(length(add_motifs) > 0){
    start_add <- nrow(subset_motifs)+1
    end_add <- nrow(subset_motifs)+length(add_motifs)

    subset_motifs[c(start_add:end_add), ] <- NA 
    subset_motifs[c(start_add:end_add), 5 ] <- add_motifs
    }

    subset_motifs$sample_name <- rep(plots.v[i], rep = nrow(subset_motifs))

    results_df <- rbind(results_df, subset_motifs)

    }

    
#plot pval vs motif type,
  compare_dosage_0 <- plots.v 

  #labels based on variable established at top of script
  
  if(new_names[t] == "RD"){
    new.labs <- c("RD, 0 ng/mL dox","RD-iP3F, 0 ng/mL dox","RD-iP3F, 75 ng/mL dox","RD-iP3F, 250 ng/mL dox","RD-iP3F, 500 ng/mL dox")
  }
  if(new_names[t] == "SMS"){
    new.labs <- c("SMS-CTR, 0 ng/mL dox","SMS-CTR-iP3F, 0 ng/mL dox","SMS-CTR-iP3F, 75 ng/mL dox","SMS-CTR-iP3F, 250 ng/mL dox"," SMS-CTR-iP3F, 500 ng/mL dox")
  }
  
  
  names(new.labs) <- compare_dosage_0
  
  results_df$sample_name <- factor(results_df$sample_name, levels=compare_dosage_0) #this fixes ordering of plots
  results_df <- results_df[-which(is.na(results_df$Motif.Name)),]

  #boxplots without factor names
  p2 <- ggplot(results_df, aes(x=type_only, y=minus_Log.P.value))+
    geom_boxplot(aes(color=type_only))+
    geom_jitter(aes(color=type_only), alpha=0.5)+
    scale_color_manual(values= col.v , drop=F)+
    theme_classic(base_size=25)+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),plot.title = element_text(hjust=0.5), strip.text = element_text(size=20), panel.spacing.y = unit(rep(0.1,length(plots.v)-1), "lines"))+
    labs(x="Motif type",y="-log10(p-value)",title="")+
     guides(color = "none") +
    facet_wrap(vars(sample_name), labeller=labeller(sample_name = new.labs), ncol=1, scales="free_y") 

  #plot dimensions based on number of plots
  if(length(plots.v) == 6){
    height_var2 <- 12
  }
  if(length(plots.v) == 5){
    height_var2 <- 13
  }
  if(length(plots.v) == 4){
    height_var2 <- 8
  }
  if(length(plots.v) == 3){
    height_var2 <- 6
  }

  png(paste0("~/dosage_manuscript/revision_figures/",file_var[t], "_pval_vs_motif_type.png" ), width = 5.5, height = height_var2, units = "in", res = 200, bg = "transparent", type = "cairo-png")
  print(p2)
  dev.off()

  #for Supplemental Table 2
  write.table(results_df, paste0("~/dosage_manuscript/revision_figures/motif_tables/",file_var[t],".tsv"), quote=F, row.names=F, sep="\t" )

} #end for loop over plots_list