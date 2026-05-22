# finds overlaps between dosage ATAC (RD, SMS) & ChIP peaks (Dbt)
#used for Homer downstream & for Supplemental Figure 12C

rm(list = ls())
graphics.off()

library(tidyverse)
library(GenomicRanges)

#import ATAC peaks files, from find_dosage_specific_peaks.R
atac_peaks <- readRDS("~/RD_SMS_iP3F/ATAC_seq/dosage_specific_peak_categories_plus.rds")

#import ChIP peaks files, from figure 3 analyses
chip_peaks <- readRDS("~/dbt_dosage_two_reps/peak_categories_GR_090624.rds")

# overlap GR objects (assign each ATAC peak a ChIP peak if present in dataset)

#get indices of overlaps in atac_peaks
chip_atac_overlap <- findOverlaps(atac_peaks, chip_peaks)
overlap_indices <- queryHits(findOverlaps(atac_peaks, chip_peaks)) 

atac_peaks$in_chip_peaks <- F

atac_peaks$in_chip_peaks[overlap_indices] <- T

#add chip peak metadata to ATAC peaks
atac_peaks$category_string_chip <- NA
atac_peaks$peak_in_75 <- NA
atac_peaks$peak_in_150 <- NA
atac_peaks$peak_in_250 <- NA
atac_peaks$peak_in_500 <- NA

atac_peaks$category_string_chip[overlap_indices] <- chip_peaks$category_string[subjectHits(chip_atac_overlap)]
atac_peaks$peak_in_75[overlap_indices] <- chip_peaks$peak_in_75[subjectHits(chip_atac_overlap)]
atac_peaks$peak_in_150[overlap_indices] <- chip_peaks$peak_in_150[subjectHits(chip_atac_overlap)]
atac_peaks$peak_in_250[overlap_indices] <- chip_peaks$peak_in_250[subjectHits(chip_atac_overlap)]
atac_peaks$peak_in_500[overlap_indices] <- chip_peaks$peak_in_500[subjectHits(chip_atac_overlap)]

atac_peaks_subset <- atac_peaks[which(atac_peaks$in_chip_peaks)]

categories_SMS <- unique(atac_peaks_subset$category_string_SMS)

#save bed files for each category
for(g in 1:length(categories_SMS)){

    subset_peaks <- atac_peaks_subset[which(atac_peaks_subset$category_string_SMS == categories_SMS[g]),]

    write.table(subset_peaks, file=paste0("~/RD_SMS_iP3F/ATAC_seq/peak_categories_chip/",categories_SMS[g],"_peaks_040926.bed"),sep="\t", col.names=F, row.names=F,quote=F)


}


