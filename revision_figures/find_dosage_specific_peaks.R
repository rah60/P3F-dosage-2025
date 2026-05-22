# finds overlaps between dosage peaks

rm(list = ls())
graphics.off()

library(GenomicRanges)
library(tidyverse)

## identifies dosage-specific peaks for RD, SMS, RD-iP3F, SMS-iP3F ATAC-seq for Supplemental Figure 12

#file prefixes of peak files
sample_prefix <- c("RD", "RD_iP3F_0", "RD_iP3F_75", "RD_iP3F_250","RD_iP3F_500", "SMS", "SMS_iP3F_0", "SMS_iP3F_75", "SMS_iP3F_250","SMS_iP3F_500")
file_path <- "~/RD_SMS_iP3F/ATAC_seq/"
sample_suffix <- c("/peak/rep1_vs_rep2/rep1_vs_rep2.idr0.05.bfilt.narrowPeak.gz")

results_list <- vector("list")

#read in bed files to GR objects, store in list
for(i in 1:length(sample_prefix)){

    file_name <- paste0(file_path, sample_prefix[i], sample_suffix)
    temp_table <- read.table(file_name, sep="\t",header=F)
    temp_table$relative_peak_center <- temp_table$V10
    results_list[[i]] <- makeGRangesFromDataFrame(temp_table, seqnames.field="V1",start.field="V2", end.field="V3", keep.extra.columns=T, starts.in.df.are.0based=T)

}

#make a GR object of all peaks from all input files
all_peaks <- GenomicRanges::reduce(unlist(GRangesList(results_list)), min.gapwidth=0L) #min gap is to combine overlaps that are immediately adjacent, as bedtools appears to do

#### examined peaks across RD & SMS, but decided this didn't make sense ultimately
#establish metadata slots in GR object for each dosage
all_peaks$peak_in_RD <- F
all_peaks$peak_in_RD_0 <- F
all_peaks$peak_in_RD_75 <- F
all_peaks$peak_in_RD_250 <- F
all_peaks$peak_in_RD_500 <- F
all_peaks$peak_in_SMS <- F
all_peaks$peak_in_SMS_0 <- F
all_peaks$peak_in_SMS_75 <- F
all_peaks$peak_in_SMS_250 <- F
all_peaks$peak_in_SMS_500 <- F

#go through each dosage
for(j in 1:length(sample_prefix)){

    #find overlaps with dosage peaks, get the indices of the overlaps in all_peaks

   overlap_indices <- queryHits(findOverlaps(all_peaks, results_list[[j]])) 

    #increment metadata in all_peaks based on which dosage

    if((j == 1) ){
    all_peaks$peak_in_RD[overlap_indices] <- T
   }
   if((j == 2 )){ 
    all_peaks$peak_in_RD_0[overlap_indices] <- T
   }
   if((j == 3)){
    all_peaks$peak_in_RD_75[overlap_indices] <- T
   }
   if((j == 4 )){
    all_peaks$peak_in_RD_250[overlap_indices] <- T
   }
   if((j == 5)){
    all_peaks$peak_in_RD_500[overlap_indices] <- T
   }
   if((j == 6)){
    all_peaks$peak_in_SMS[overlap_indices] <- T
   }
   if((j == 7)){
    all_peaks$peak_in_SMS_0[overlap_indices] <- T
   }
   if((j == 8)){
    all_peaks$peak_in_SMS_75[overlap_indices] <- T
   }
   if((j == 9)){
    all_peaks$peak_in_SMS_250[overlap_indices] <- T
   }
   if((j == 10)){
    all_peaks$peak_in_SMS_500[overlap_indices] <- T
   }


}

#make strings to represent all the dosages each peak was called in

all_peaks$category_string <- "no category"
for(h in 1:length(all_peaks)){
    
    cat(h, "\n")
    
    build_category_string <- c()
    the_peak <- all_peaks[h]

    if(all_peaks[h]$peak_in_RD){
    build_category_string <- c(build_category_string, "RD")
   }
   if(all_peaks[h]$peak_in_RD_0){
    build_category_string <- c(build_category_string, "RD-0")
   }
    if(all_peaks[h]$peak_in_RD_75){
    build_category_string <- c(build_category_string, "RD-75")
   }
    if(all_peaks[h]$peak_in_RD_250){
    build_category_string <- c(build_category_string, "RD-250")
   }
   if(all_peaks[h]$peak_in_RD_500){
    build_category_string <- c(build_category_string, "RD-500")
   }
    if(all_peaks[h]$peak_in_SMS){
    build_category_string <- c(build_category_string, "SMS")
   }
   if(all_peaks[h]$peak_in_SMS_0){
    build_category_string <- c(build_category_string, "SMS-0")
   }
   if(all_peaks[h]$peak_in_SMS_75){
    build_category_string <- c(build_category_string, "SMS-75")
   }
   if(all_peaks[h]$peak_in_SMS_250){
    build_category_string <- c(build_category_string, "SMS-250")
   }
   if(all_peaks[h]$peak_in_SMS_500){
    build_category_string <- c(build_category_string, "SMS-500")
   }

   category_string <- paste(build_category_string, collapse="_")

  all_peaks[h]$category_string <- category_string

  rm(category_string, build_category_string)

}

#save GR object
saveRDS(all_peaks,"~/RD_SMS_iP3F/ATAC_seq/dosage_specific_peak_categories.rds")

#convert to dataframe
save_GR <- as.data.frame(all_peaks)

#save as a tsv
write.table(save_GR, file="~/RD_SMS_iP3F/ATAC_seq/dosage_specific_peak_categories.tsv",sep="\t")

categories <- unique(all_peaks$category_string)

#save bed files for each category
for(g in 1:length(categories)){

    subset_peaks <- all_peaks[which(all_peaks$category_string == categories[g]),]

    write.table(subset_peaks, file=paste0("~/RD_SMS_iP3F/ATAC_seq/peak_categories/",categories[g],"_peaks_040226.bed"),sep="\t", col.names=F, row.names=F,quote=F)


}

############# RD & SMS only peak categories, ultimately the peak categories I used for S12

all_peaks <- readRDS("~/RD_SMS_iP3F/ATAC_seq/dosage_specific_peak_categories.rds")

#### RD
all_peaks$category_string_RD <- "no category"
all_peaks$category_string_SMS <- "no category"

for(h in 1:length(all_peaks)){
    
    cat(h, "\n")
    
    build_category_string_RD <- c()
    build_category_string_SMS <- c()
    the_peak <- all_peaks[h]

    if(all_peaks[h]$peak_in_RD){
    build_category_string_RD <- c(build_category_string_RD, "RD")
   }
   if(all_peaks[h]$peak_in_RD_0){
    build_category_string_RD <- c(build_category_string_RD, "RD-0")
   }
    if(all_peaks[h]$peak_in_RD_75){
    build_category_string_RD <- c(build_category_string_RD, "RD-75")
   }
    if(all_peaks[h]$peak_in_RD_250){
    build_category_string_RD <- c(build_category_string_RD, "RD-250")
   }
   if(all_peaks[h]$peak_in_RD_500){
    build_category_string_RD <- c(build_category_string_RD, "RD-500")
   }
    if(all_peaks[h]$peak_in_SMS){
    build_category_string_SMS <- c(build_category_string_SMS, "SMS")
   }
   if(all_peaks[h]$peak_in_SMS_0){
    build_category_string_SMS <- c(build_category_string_SMS, "SMS-0")
   }
   if(all_peaks[h]$peak_in_SMS_75){
    build_category_string_SMS <- c(build_category_string_SMS, "SMS-75")
   }
   if(all_peaks[h]$peak_in_SMS_250){
    build_category_string_SMS <- c(build_category_string_SMS, "SMS-250")
   }
   if(all_peaks[h]$peak_in_SMS_500){
    build_category_string_SMS <- c(build_category_string_SMS, "SMS-500")
   }

   category_string_RD <- paste(build_category_string_RD, collapse="_")
   category_string_SMS <- paste(build_category_string_SMS, collapse="_")

  all_peaks[h]$category_string_RD <- category_string_RD
  all_peaks[h]$category_string_SMS <- category_string_SMS

  rm(category_string_RD, build_category_string_RD, category_string_SMS, build_category_string_SMS)

}

# saveRDS(all_peaks,"~/RD_SMS_iP3F/ATAC_seq/dosage_specific_peak_categories_plus.rds")
all_peaks <- readRDS("~/RD_SMS_iP3F/ATAC_seq/dosage_specific_peak_categories_plus.rds")

#convert to dataframe
#save_GR <- as.data.frame(all_peaks)

#save as a tsv
#write.table(save_GR, file="~/RD_SMS_iP3F/ATAC_seq/dosage_specific_peak_categories_plus.tsv",sep="\t")
save_GR <- read.table("~/RD_SMS_iP3F/ATAC_seq/dosage_specific_peak_categories_plus.tsv",sep="\t",header=T)

categories_RD <- unique(all_peaks$category_string_RD)
categories_SMS <- unique(all_peaks$category_string_SMS)

#save bed files for each category
for(g in 1:length(categories_RD)){

    subset_peaks <- all_peaks[which(all_peaks$category_string_RD == categories_RD[g]),]

    write.table(subset_peaks, file=paste0("~/RD_SMS_iP3F/ATAC_seq/peak_categories_RD/",categories_RD[g],"_peaks_040626.bed"),sep="\t", col.names=F, row.names=F,quote=F)


}

for(g in 1:length(categories_SMS)){

    subset_peaks <- all_peaks[which(all_peaks$category_string_SMS == categories_SMS[g]),]

    write.table(subset_peaks, file=paste0("~/RD_SMS_iP3F/ATAC_seq/peak_categories_SMS/",categories_SMS[g],"_peaks_040626.bed"),sep="\t", col.names=F, row.names=F,quote=F)


}


