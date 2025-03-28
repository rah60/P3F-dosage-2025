rm(list = ls())
graphics.off()

#plot figures 3B & C, extended data figures 3A, B

library(ggplot2)
library(viridis)
library(dplyr)
library(stringr)

#set up for all plots first

all_peaks <- read.table("~/dosage_manuscript/csv/peak_categories_GR_090624.tsv",sep="\t")

#unique peaks and total peaks dataframe set-up
dosage <- c("0","75","150","250","500","1000","75_150_250_500_1000")

n_cat <- all_peaks[which(all_peaks$category_string %in% dosage),] %>%
    count(category_string, name="unique_peak_count")

n_cat[which(n_cat$category_string == "75_150_250_500_1000"),1] <- "common"

n_cat$total_peak_count <- NA

#determined from peak files for each dosage. number of peaks per dosage total in all_peaks is slightly different due to overlapping peaks that get combined in peak_set_overlaps.R
#I decided to go with the initial number of peaks for each dosage unaffected by information from other dosages. pattern is be the same regardless
#system("wc -l ~/dbt_dosage_two_reps/75_2w_50bp_peaks.bed"), etc used to determine peak number in peak files

n_cat[which(n_cat$category_string == "75"),3] <- 24240

n_cat[which(n_cat$category_string == "150"),3] <- 39703

n_cat[which(n_cat$category_string == "250"),3] <- 35536

n_cat[which(n_cat$category_string == "500"),3] <- 29738

n_cat[which(n_cat$category_string == "1000"),3] <- 6919

n_cat[which(n_cat$category_string == "0"),3] <- 1261

#two dosage overlap peaks dataframe
dosage_two <- c("150_500", "150_250","250_500", "75_150","75_250" , "500_1000" ,"75_500" ,"250_1000" ,"150_1000"  ,"75_0","75_1000")
combo_two <- all_peaks[which(all_peaks$category_string %in% dosage_two),]
n_cat_2 <- combo_two %>%
    count(category_string)

#three dosage overlap peaks dataframe
dosage_three <- c("75_150_500","150_250_500","75_150_250","75_150_1000","250_500_1000" ,"75_250_500" ,"150_500_1000" ,"150_250_1000","75_500_0","75_150_0")
combo_three <- all_peaks[which(all_peaks$category_string %in% dosage_three),]
n_cat_3 <- combo_three %>%
    count(category_string)

col.v1 <- viridis(8)
col.v <- col.v1[1:7]
names(col.v) <- c("0","75","150","250","500","1000","common")

#Figure 3C
dbt.df <- n_cat

#unique peaks
p1 <- ggplot(dbt.df)+
  geom_col(aes(x=category_string, y=unique_peak_count, fill=category_string), show.legend = F )+ 
  scale_fill_manual(values = col.v[1:7])+ 
  scale_x_discrete(limits =  c("0","75","150","250","500","1000","common"), labels=c("unique to 0 ng/mL","unique to 75 ng/mL","unique to 150 ng/mL","unique to 250 ng/mL","unique to 500 ng/mL","unique to 1000 ng/mL","common") )+
  labs(y="PAX3::FOXO1 sites",x="")+
  scale_y_continuous(expand = c(0,0), limits=c(0,max(n_cat$unique_peak_count)))+
  theme_classic(base_size = 30)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

png("~/dosage_manuscript/figure_3/dbt_mycn_ip3f_unique_peak_color_by_dosage.png", width = 7, height = 8, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()

#Figure 3B

#total peaks
p2 <- ggplot(dbt.df)+
  geom_col(aes(x=category_string, y=total_peak_count, fill=category_string), show.legend = F )+
  scale_fill_manual(values = col.v[1:6])+
  scale_x_discrete(limits =  c("0","75","150","250","500","1000"))+
  labs(y="PAX3::FOXO1 sites",x="Doxycycline dose (ng/mL)")+
  scale_y_continuous(expand = c(0,0), limits=c(0,max(n_cat$total_peak_count)))+
  theme_classic(base_size = 30)

png("~/dosage_manuscript/figure_3/dbt_mycn_ip3f_total_peak_color_by_dosage.png", width = 7, height = 5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()



#two sample overlaps

dosage_two_new <- str_replace(dosage_two, "_", " & ")

#Extended Data Figure 3A

p3 <- ggplot(n_cat_2)+
  geom_col(aes(x=reorder(category_string, -n), y=n, fill=category_string), show.legend = F )+ #fill=dumb_color_var
  scale_fill_viridis(discrete = T, end=0.8)+
  scale_x_discrete(limits =  dosage_two, labels=dosage_two_new  )+
  labs(y="PAX3::FOXO1 sites",x="Doxycycline (ng/ml)")+
  scale_y_continuous(expand = c(0,0))+
  theme_classic(base_size = 30)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

png("~/dosage_manuscript/figure_3/dbt_mycn_ip3f_two_dosages_S3.png", width = 8, height = 8, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p3)
dev.off()

#three sample overlaps

#Extended Data Figure 3B

dosage_three_new <- dosage_three %>%
                    str_replace( "_", ", ") %>%
                    str_replace( "_", ", & ")


p4 <- ggplot(n_cat_3)+
  geom_col(aes(x=reorder(category_string, -n), y=n, fill=category_string), show.legend = F )+ #fill=dumb_color_var
  scale_fill_viridis(discrete = T, end=0.8)+
  scale_x_discrete(limits = dosage_three, labels=dosage_three_new  )+
  labs(y="PAX3::FOXO1 sites",x="Doxycycline (ng/ml)")+
  scale_y_continuous(expand = c(0,0))+
  theme_classic(base_size = 30)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

png("~/dosage_manuscript/figure_3/dbt_mycn_ip3f_three_dosages_S3.png", width = 8, height = 8, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p4)
dev.off()
