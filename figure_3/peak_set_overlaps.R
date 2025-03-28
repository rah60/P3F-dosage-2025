# finds overlaps between dosage peaks
#peak files were extended from ENCODE output using add_50bp_to_peaks.sh prior to this script

rm(list = ls())
graphics.off()

library(GenomicRanges)
library(tidyverse)

#file prefixes of peak files
sample_prefix <- c("75_2w_50bp_peaks", "150_2w_50bp_peaks", "250_2w_50bp_peaks", "500_2w_50bp_peaks", "1000_2w_50bp_peaks","0_2w_50bp_peaks")
file_path <- "~/dbt_dosage_two_reps/"

results_list <- vector("list")

#read in bed files to GR objects, store in list
for(i in 1:length(sample_prefix)){

    file_name <- paste0(file_path, sample_prefix[i], ".bed")
    temp_table <- read.table(file_name, sep="\t",header=F)
    temp_table$relative_peak_center <- temp_table$V10
    results_list[[i]] <- makeGRangesFromDataFrame(temp_table, seqnames.field="V1",start.field="V2", end.field="V3", keep.extra.columns=T, starts.in.df.are.0based=T)

}

#make a GR object of all peaks from all input files
all_peaks <- GenomicRanges::reduce(unlist(GRangesList(results_list)), min.gapwidth=0L) #min gap is to combine overlaps that are immediately adjacent, as bedtools appears to do

#establish metadata slots in GR object for each dosage
all_peaks$peak_in_75 <- F
all_peaks$peak_in_150 <- F
all_peaks$peak_in_250 <- F
all_peaks$peak_in_500 <- F
all_peaks$peak_in_1000 <- F
all_peaks$peak_in_0 <- F

#go through each dosage
for(j in 1:length(sample_prefix)){

    #find overlaps with dosage peaks, get the indices of the overlaps in all_peaks

   overlap_indices <- queryHits(findOverlaps(all_peaks, results_list[[j]])) 

    #increment metadata in all_peaks based on which dosage

    if((j == 1) ){
    all_peaks$peak_in_75[overlap_indices] <- T
   }
   if((j == 2 )){ 
    all_peaks$peak_in_150[overlap_indices] <- T
   }
   if((j == 3)){
    all_peaks$peak_in_250[overlap_indices] <- T
   }
   if((j == 4 )){
    all_peaks$peak_in_500[overlap_indices] <- T
   }
   if((j == 5)){
    all_peaks$peak_in_1000[overlap_indices] <- T
   }
   if((j == 6)){
    all_peaks$peak_in_0[overlap_indices] <- T
   }


}

#make strings to represent all the dosages each peak was called in

all_peaks$category_string <- "no category"
for(h in 1:length(all_peaks)){
    
    cat(h, "\n")
    
    build_category_string <- c()
    the_peak <- all_peaks[h]

    if(all_peaks[h]$peak_in_75){
    build_category_string <- c(build_category_string, "75")
   }
   if(all_peaks[h]$peak_in_150){
    build_category_string <- c(build_category_string, "150")
   }
    if(all_peaks[h]$peak_in_250){
    build_category_string <- c(build_category_string, "250")
   }
    if(all_peaks[h]$peak_in_500){
    build_category_string <- c(build_category_string, "500")
   }
   if(all_peaks[h]$peak_in_1000){
    build_category_string <- c(build_category_string, "1000")
   }
    if(all_peaks[h]$peak_in_0){
    build_category_string <- c(build_category_string, "0")
   }

   category_string <- paste(build_category_string, collapse="_")

  all_peaks[h]$category_string <- category_string

  rm(category_string, build_category_string)

}

#save GR object
saveRDS(all_peaks,"~/dosage_manuscript/rds/peak_categories_GR_090624.rds")

#convert to dataframe
save_GR <- as.data.frame(all_peaks)

#save as a tsv
write.table(save_GR, file="~/dosage_manuscript/csv/peak_categories_GR_090624.tsv",sep="\t")

categories <- unique(all_peaks$category_string)

#save bed files for each category
for(g in 1:length(categories)){

    subset_peaks <- all_peaks[which(all_peaks$category_string == categories[g]),]

    write.table(subset_peaks, file=paste0("~/dbt_dosage_two_reps/peak_categories/",categories[g],"_peaks_2w_090624.bed"),sep="\t", col.names=F, row.names=F,quote=F)


}


#QC checks
# categories <- unique(all_peaks$category_string) #with zero, 50 categories

# #analyzing peak category distributions

# all_peaks <- read.table("~/dbt_dosage_two_reps/peak_categories_GR_090624.tsv",sep="\t")

# #how many peaks are there total?
# nrow(all_peaks) 

# nrow(unique(all_peaks))
# #still the name number of peaks

# # #how many peaks are in each category?

# n_cat <- all_peaks %>%
#     count(category_string)

