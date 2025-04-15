rm(list = ls()) 
graphics.off()

#categorize repliATAC-seq peaks, pulse to use in Fig 5/6

library(GenomicRanges)
library(tidyverse)

sample_prefix <- c("Dox0_pulse","Dox75_pulse","Dox500_pulse") 
file_path <- "/home/gdstantonlab/lab/GSL-RH-3690/ENCODE/results/"

results_list <- vector("list")
#read in bed files to GR objects, store in list
for(i in 1:length(sample_prefix)){
  
  file_name <- paste0(file_path, sample_prefix[i],"/",sample_prefix[i], "_IDR.bed") #IDR peaks from ENCODE pipeline
  temp_table <- read.table(file_name, sep="\t",header=F)
  temp_table$relative_peak_center <- temp_table$V10
  results_list[[i]] <- makeGRangesFromDataFrame(temp_table, seqnames.field="V1",start.field="V2", end.field="V3", keep.extra.columns=T, starts.in.df.are.0based=T)
  
}

#unlist into one GRanges object
all_peaks <- GenomicRanges::reduce(unlist(GRangesList(results_list)), min.gapwidth=0L) #min gap is to combine overlaps that are immediately adjacent, as bedtools appears to do

all_peaks$peak_in_0_pulse <- F
all_peaks$peak_in_75_pulse  <- F
all_peaks$peak_in_500_pulse <- F

for(j in 1:length(sample_prefix)){
  
  #find overlaps with dosage peaks, get the indices of the overlaps in all_peaks
  
  overlap_indices <- queryHits(findOverlaps(all_peaks, results_list[[j]])) 
  
  #increment metadata in all_peaks based on which dosage
  
  if((j == 1) ){
    
    all_peaks$peak_in_0_pulse[overlap_indices] <- T
    
  }
  if((j == 2 )){ 
    all_peaks$peak_in_75_pulse[overlap_indices] <- T
  }
  if((j == 3)){
    all_peaks$peak_in_500_pulse[overlap_indices] <- T
  }
  
  
}

#name peak categories based on which samples they are called peaks in
all_peaks$category_string <- "no category"
for(h in 1:length(all_peaks)){
  
  cat(h, "\n")
  
  build_category_string <- c()
  the_peak <- all_peaks[h]
  
  if(all_peaks[h]$peak_in_0_pulse){
    build_category_string <- c(build_category_string, "0p")
  }
  if(all_peaks[h]$peak_in_75_pulse){
    build_category_string <- c(build_category_string, "75p")
  }
  if(all_peaks[h]$peak_in_500_pulse){
    build_category_string <- c(build_category_string, "500p")
  }
  
  category_string <- paste(build_category_string, collapse="_")
  
  all_peaks[h]$category_string <- category_string
  
  rm(category_string, build_category_string)
  
}

saveRDS(all_peaks,"~/dosage_manuscript/rds/peak_categories_GR_103024_pulse.rds")

save_GR <- as.data.frame(all_peaks)

write.table(save_GR, file="~/dosage_manuscript/csv/peak_categories_GR_103024_pulse.tsv",sep="\t")

categories <- unique(all_peaks$category_string) #with zero, 50 categories


#save bed files for each category, for motif calling
for(g in 1:length(categories)){
  
  subset_peaks <- all_peaks[which(all_peaks$category_string == categories[g]),]
  
  write.table(subset_peaks, file=paste0("~/repliATAC/peak_categories/",categories[g],"_peaks_103024_pulse.bed"),sep="\t", col.names=F, row.names=F,quote=F)
  
  
}



all_peaks <- read.table("~/dosage_manuscript/csv/peak_categories_GR_103024_pulse.tsv",sep="\t")

#total IDR peaks
len_0 <- length(results_list[[1]])
len_75 <- length(results_list[[2]])
len_500 <- length(results_list[[3]])

#unique IDR peaks
select_categories <- c("0p", "75p","500p" , "0p_75p_500p")

combo_two <- all_peaks[which(all_peaks$category_string %in% select_categories),]

n_cat <- combo_two %>%
  count(category_string)

#plot total & unique peaks
dosage <- c("0","75","500","common")
total_peak_count <- c(len_0, len_75, len_500, NA)
unique_peak_count <- c( n_cat[1,2],n_cat[4,2] , n_cat[3,2], n_cat[2,2])

dbt.df <- as.data.frame(cbind(dosage,unique_peak_count,total_peak_count))
dbt.df$unique_peak_count <- as.numeric(unique_peak_count)
dbt.df$total_peak_count <- as.numeric(total_peak_count)

library(viridis)
col.v1 <- viridis(8)
col.v <- c(col.v1[1],col.v1[2], col.v1[5], col.v1[7])
names(col.v) <- dosage


########same y-axis as chase

#unique peaks
#Extended Data Figure 4A, left

p1 <- ggplot(dbt.df)+
  geom_col(aes(x=dosage, y=unique_peak_count, fill=dosage), show.legend = F )+ 
  scale_fill_manual(values = col.v)+ 
  scale_x_discrete(limits =  c("0","75","500","common"), labels=c("unique to 0 ng/mL","unique to 75 ng/mL","unique to 500 ng/mL","common") )+
  labs(y="Accessible sites\nin nascent chromatin",x="", title="Nascent chromatin (10 min)")+
  scale_y_continuous(expand = c(0,0), limits=c(0,9331))+
  theme_classic(base_size = 30)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),plot.title = element_text(hjust = 0.5))

png("~/dosage_manuscript/figure_5/dbt_mycn_ip3f_unique_peak_pulse_newscale.png", width = 9, height = 7, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()

#Figure 5C, left

#total peaks
p2 <- ggplot(dbt.df)+
  geom_col(aes(x=dosage, y=total_peak_count, fill=dosage), show.legend = F )+
  scale_fill_manual(values = col.v[1:3])+
  scale_x_discrete(limits =  c("0","75","500"))+
  labs(y="Accessible sites\nin nascent chromatin",x="Doxycycline dose (ng/mL)",title = "Nascent chromatin (10 min)")+
  scale_y_continuous(expand = c(0,0), limits=c(0,35066))+
  theme_classic(base_size = 30)+
  theme(plot.title = element_text(hjust = 0.5))

png("~/dosage_manuscript/figure_5/dbt_mycn_ip3f_total_peak_pulse_newscale.png", width = 9, height = 7, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off() 