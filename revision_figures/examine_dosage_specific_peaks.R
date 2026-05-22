# look at dosage peaks

rm(list = ls())
graphics.off()

library(tidyverse)
library(viridis)

## plots total & unique peaks in RD & SMS ATAC-seq, Supplemental Figures 12A & B

all_peaks <- read.table("~/RD_SMS_iP3F/ATAC_seq/dosage_specific_peak_categories_plus.tsv",sep="\t",header=T)

n_count_RD <- all_peaks %>%
            dplyr::count(category_string_RD) %>%
            arrange(-n) %>%
            filter( category_string_RD != "")

n_count_RD$category_string_RD <- factor(n_count_RD$category_string_RD, levels=n_count_RD$category_string_RD)

n_count_SMS <- all_peaks %>%
            dplyr::count(category_string_SMS) %>%
            arrange(-n) %>%
            filter( category_string_SMS != "")

n_count_SMS$category_string_SMS <- factor(n_count_SMS$category_string_SMS, levels=n_count_SMS$category_string_SMS)

RD1 <- c("RD", "RD-0", "RD-75", "RD-250", "RD-500")

temp1 <- subset(n_count_RD, category_string_RD %in% RD1 )
temp1$category_string_RD <- factor(temp1$category_string_RD, levels = RD1)

#unique peaks RD per category
p1 <- ggplot(temp1)+
    geom_col(aes(x=category_string_RD, y=n, fill=category_string_RD), show.legend=F)+
    scale_y_continuous(expand=c(0,0))+
    scale_x_discrete(labels = c("RD, 0", "RD-iP3F, 0", "RD-iP3F, 75","RD-iP3F, 250","RD-iP3F, 500"))+
    labs(x="Cell line &\ndoxycycline (ng/mL)",y="Dosage-specific\naccessible sites",title="RD")+
    theme_classic(base_size=30)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), plot.title=element_text(hjust=0.5))

png("~/dosage_manuscript/revision_figures/RD_unique_peak_count.png", width = 6, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()

SMS1 <- c("SMS", "SMS-0", "SMS-75", "SMS-250", "SMS-500")

temp2 <- subset(n_count_SMS, category_string_SMS %in% SMS1 )
temp2$category_string_SMS <- factor(temp2$category_string_SMS, levels = SMS1)

#unique peaks SMS per category
p2 <- ggplot(temp2)+
    geom_col(aes(x=category_string_SMS, y=n, fill=category_string_SMS), show.legend=F)+
    scale_y_continuous(expand=c(0,0))+
    scale_x_discrete(labels = c("SMS-CTR, 0", "SMS-CTR-iP3F, 0", "SMS-CTR-iP3F, 75","SMS-CTR-iP3F, 250","SMS-CTR-iP3F, 500"))+
    labs(x="Cell line &\ndoxycycline (ng/mL)",y="Dosage-specific\naccessible sites",title="SMS-CTR")+
    theme_classic(base_size=30)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), plot.title=element_text(hjust=0.5))

png("~/dosage_manuscript/revision_figures/SMS_unique_peak_count.png", width = 6, height = 7, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()



######### total peaks

#total peaks by idr files from ENCODE pipeline
#RD 151452
#RD 0 165280
#RD 75 172060
#RD 250 142980
#RD 500 161121
#SMS 155725
#SMS 0 157356
#SMS 75 157876
#SMS 250 154111
#SMS 500 139899

temp1$n_total <- NA
temp2$n_total <- NA

temp1[which(temp1$category_string_RD == "RD"),3] <- 151452
temp1[which(temp1$category_string_RD == "RD-0"),3] <- 165280
temp1[which(temp1$category_string_RD == "RD-75"),3] <- 172060
temp1[which(temp1$category_string_RD == "RD-250"),3] <- 142980
temp1[which(temp1$category_string_RD == "RD-500"),3] <- 161121

temp2[which(temp2$category_string_SMS == "SMS"),3] <- 155725
temp2[which(temp2$category_string_SMS == "SMS-0"),3] <- 157356
temp2[which(temp2$category_string_SMS == "SMS-75"),3] <- 157876
temp2[which(temp2$category_string_SMS == "SMS-250"),3] <- 154111
temp2[which(temp2$category_string_SMS == "SMS-500"),3] <- 139899

p3 <- ggplot(temp1)+
    geom_col(aes(x=category_string_RD, y=n_total, fill=category_string_RD), show.legend=F)+
    scale_y_continuous(expand=c(0,0))+
    scale_x_discrete(labels = c("RD, 0", "RD-iP3F, 0", "RD-iP3F, 75","RD-iP3F, 250","RD-iP3F, 500"))+
    labs(x="Cell line &\ndoxycycline (ng/mL)",y="Accessible sites",title="RD")+
    theme_classic(base_size=30)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),plot.title=element_text(hjust=0.5))

png("~/dosage_manuscript/revision_figures/RD_total_peak_count.png", width = 6, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p3)
dev.off()

p4 <- ggplot(temp2)+
    geom_col(aes(x=category_string_SMS, y=n_total, fill=category_string_SMS), show.legend=F)+
    scale_y_continuous(expand=c(0,0))+
    scale_x_discrete(labels = c("SMS-CTR, 0", "SMS-CTR-iP3F, 0", "SMS-CTR-iP3F, 75","SMS-CTR-iP3F, 250","SMS-CTR-iP3F, 500"))+
    labs(x="Cell line & \ndoxycycline (ng/mL)",y="Accessible sites",title="SMS-CTR")+
    theme_classic(base_size=30)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),plot.title=element_text(hjust=0.5))

png("~/dosage_manuscript/revision_figures/SMS_total_peak_count.png", width = 6, height = 7, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p4)
dev.off()