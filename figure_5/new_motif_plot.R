#clear workspace
rm(list = ls())
graphics.off()

library(gplots)
library(ggplot2)
library(ggrepel)
library(viridis)
library(tidyverse)
library(lemon)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(RColorBrewer)

#Figure 5E, Extended Data Figure 5C, Figure 6C, Extended Data Figure 6C
#see figure_3/motif_plot.R for commented version of highly similar script

plots_list <- list(c("0p_peaks_103024_pulse", "75p_peaks_103024_pulse", "500p_peaks_103024_pulse", "0p_75p_500p_peaks_103024_pulse"),
         c("0c_peaks_103024_chase", "75c_peaks_103024_chase", "500c_peaks_103024_chase", "0c_75c_500c_peaks_103024_chase"),
         c("cluster_1", "cluster_2", "cluster_3")
          )

new_names <- c("four", "four", "three")

file_var <- c("pulse", "chase","pulse_clusters")

dosages <- c("0_75_500","0_75_500","0_75_500")

path_prefix <- "~/repliATAC/homer_results/"

#functions

#make function to take dosage & return vector of gene names
#### then I can split strings on "_" to get dosages, and get as many gene vectors as needed for given subset

read_motif_data <- function(x, path_prefix = "~/repliATAC/homer_results/"){ #x is data name prefix
    results_path <- paste(path_prefix,x,"/knownResults.txt",sep="") 
    results <- read.table(results_path, header = T, sep="\t", comment.char = "") 
    results
  }

get_homer_motif_types_df <- function(x){ #x is dataframe
    temp <- x$Motif.Name #modified to not select top 50 based on prior filtering steps
    temp1 <- strsplit(temp, "/")
    temp2 <- strsplit(sapply(temp1,"[[",1),"[(]")
    temp3 <- strsplit(sapply(temp2, "[",2),"[)]")
    unlist(sapply(temp3, "[[",1))
    }

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
    
get_homer_motif_types <- function(x){ #x is vector of homer motif names
    temp1 <- strsplit(x, "/")
    temp2 <- strsplit(sapply(temp1,"[[",1),"[(]")
    temp3 <- strsplit(sapply(temp2, "[",2),"[)]")
    sapply(temp3, "[[",1)
    }

get_homer_motif_names <- function(x){ #x is vector of homer motif names
    temp1 <- strsplit(x, "/")
    temp2 <- strsplit(sapply(temp1,"[[",1),"[(]")
    sapply(temp2, "[",1)
    }
  
#import list of expressed genes from bulk RNA-seq generated in make_expressed_gene_list.R
#and get genes expressed at a dosage  
get_dosage_gene_vector <- function(x, genes_rds = "~/dosage_manuscript/rds/gene_expressed_at_dosages.rds"){ #x is one dosage
  genes_df <- readRDS(genes_rds)
  subset_x <- genes_df[which(genes_df$Dosage == x),]
  subset_x$gene_name
}

for(t in 1:length(plots_list)){
  
  plots.v <- plots_list[[t]]

  motifs_list <- lapply(plots.v, read_motif_data)

#get all potential motif names across samples
  all_factor_names_motif_type <- lapply(motifs_list, get_list_motif_factor_names_df) %>%
                      unlist() 
                      
  all_factor_names <- unique(names(all_factor_names_motif_type))

if("NF1" %in% all_factor_names){
  all_factor_names[which(all_factor_names == "NF1")] <- "NFIA" } #fixing homer naming issue

names(all_factor_names) <- toupper(all_factor_names) #to avoid case issues
P3F_index <- which(all_factor_names == "PAX3:FKHR-fusion")

#have to use uppercase here #need to check what doesn't get converted to make sure it's ok/makes sense
factors_geneid <- ensembldb::select(EnsDb.Hsapiens.v86, keys=toupper(all_factor_names), keytype="SYMBOL", columns=c("SYMBOL", "GENEID"), ignore.case=T)

#get dosages
temp1 <- strsplit(dosages[t], "_")
dosages.v <- unique(unlist(temp1))

#compare to the list of expressed genes in scRNA/RNA seq (get with function for that)
all_exp_genes <- lapply(dosages.v, get_dosage_gene_vector) %>%
                      unlist() %>%
                      unique()

get_these <- factors_geneid[which(factors_geneid$GENEID %in% all_exp_genes),1]

keep_names <- all_factor_names[get_these]
keep_names <- c(keep_names, all_factor_names[P3F_index]) #make sure P3F remains in list, b/c no ensembl id for it

keep_motifs <- unique(all_factor_names_motif_type[keep_names])

#for any motif names that don't have a match in the expressed gene set, drop that from the dataset (or save a list here & drop later in the script)
 for(w in 1:length(motifs_list)){

  keep_index <- which(get_list_motif_factor_names_df(motifs_list[[w]]) %in% keep_motifs)
  motifs_list[[w]] <- motifs_list[[w]][keep_index,]
  motifs_list[[w]] <- motifs_list[[w]][1:50,]

 }

#get all potential motif types across samples

  all_motif_types <- lapply(motifs_list, get_homer_motif_types_df) %>%
                      unlist() %>%
                      unique()

if("TEA" %in% all_motif_types){
  all_motif_types <- all_motif_types[-which(all_motif_types == "TEA")]}

if("Zf" %in% all_motif_types){
  all_motif_types <- all_motif_types[-which(all_motif_types == "Zf")]
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

    
#plot pval vs motif type, and then add text with label as specific motif name
  compare_dosage_0 <- plots.v 

  if(new_names[t] == "four"){
    new.labs <- c("0 ng/ml","75 ng/ml","500 ng/ml","common")
  }
   if(new_names[t] == "three"){
    new.labs <- c("Cluster 1", "Cluster 2", "Cluster 3")
  }
  
  names(new.labs) <- compare_dosage_0
  
  results_df$sample_name <- factor(results_df$sample_name, levels=compare_dosage_0) #this fixes ordering of plots

  #plot motif names
  p1 <- ggplot(results_df, aes(x=type_only, y=minus_Log.P.value))+
    geom_point(aes(color=type_only), alpha=0.5)+
    scale_color_manual(values= col.v , drop=F)+
    theme_classic(base_size=20)+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),plot.title = element_text(hjust=0.5), strip.text = element_text(size=20))+
    labs(x="Motif type",y="-log10(p-value)",title="")+
    geom_text_repel(aes(color=type_only, label=name_only), box.padding = 0.25, show.legend = F, max.overlaps=nrow(subset_motifs)) +
    guides(color = "none") +
    facet_rep_wrap(vars(sample_name), labeller=labeller(sample_name = new.labs), ncol=1, scales="free_y") #need to convert this to facet_wrap? to switch labels to top at least for some plots

  #boxplot, no motif names
  p2 <- ggplot(results_df, aes(x=type_only, y=minus_Log.P.value))+
    geom_boxplot(aes(color=type_only))+
    geom_jitter(aes(color=type_only), alpha=0.5)+
    scale_color_manual(values= col.v , drop=F)+
    theme_classic(base_size=14)+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),plot.title = element_text(hjust=0.5), strip.text = element_text(size=12))+
    labs(x="Motif type",y="-log10(p-value)",title="")+
     guides(color = "none") +
    facet_rep_wrap(vars(sample_name), labeller=labeller(sample_name = new.labs), ncol=1, scales="free_y") 

  if(length(plots.v) == 6){
    height_var1 <- 24
    height_var2 <- 12
  }
  if(length(plots.v) == 5){
    height_var1 <- 20
    height_var2 <- 10
  }
  if(length(plots.v) == 4){
    height_var1 <- 16
    height_var2 <- 8
  }
  if(length(plots.v) == 3){
    height_var1 <- 12
    height_var2 <- 6.5
  }
  
  png(paste0("~/dosage_manuscript/figure_5/",file_var[t], "_pval_vs_motif_type_filtered_type_toplabel.png" ), width = 9, height = height_var1, units = "in", res = 200, bg = "transparent", type = "cairo-png")
  print(p1)
  dev.off()

  png(paste0("~/dosage_manuscript/figure_5/",file_var[t], "_pval_vs_motif_type_nolabel_filtered_type_toplabel.png" ), width = 4, height = height_var2, units = "in", res = 200, bg = "transparent", type = "cairo-png")
  print(p2)
  dev.off()

} #end for loop over plots_list 