
rm(list = ls())
graphics.off()

library(ggplot2)
library(readr)
library(viridis)
library(dplyr)
library(ggsignif)

#Figure 5A

col.v1 <- viridis(8)
col.v <- col.v1[1:7] # 0, 75, 150, 250, 500, 1000, common  so 0, 75, 500 is 1, 2, 5
col_3 <- c(col.v[1],col.v[2],col.v[5])

#Dbt/MYCN/iP3F
data <- read_csv("~/dosage_manuscript/csv/edu_2wk_dmp.csv") #reads in data

data$doxycycline <- as.character(data$doxycycline)

#ANOVA significance testing for edu+
data_plus <- data[which(data$edu_status=="plus"),]
res.aov <- aov(percent ~ doxycycline, data=data_plus)
summary(res.aov)

TukeyHSD(res.aov) #only 0 and 500 are significantly different, 0.0279

#significance testing for edu-
data_minus <- data[which(data$edu_status=="minus"),]
res.aov <- aov(percent ~ doxycycline, data=data_minus)
summary(res.aov)

TukeyHSD(res.aov) #again, 0 and 500 are significantly different, 0.0279

#process data for plotting
data_sum <- data %>% #summarize to get mean & std dev & standard error
  group_by(edu_status, doxycycline) %>%
  summarise(sd=sd(percent), mean=mean(percent), se=sd(percent)/sqrt(2))

data_sum <- data_sum[which(data_sum$edu_status == "plus"),]

p1 <- ggplot(data_sum, aes(x=factor(doxycycline, c("0","75","500")), y=mean, fill=doxycycline ))+ 
  geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7)+
  labs(x="PAX3::FOXO1 induction level \n(doxycycline dosage, ng/mL)", y="Percent S phase cells in population", fill="")+ #change cell type label
  geom_errorbar( aes(ymin=mean-se, ymax=mean+se),
                 width=0.3, position = position_dodge(width=0.7) )+
  scale_y_continuous(expand = c(0,0), limits=c(0,40))+
  theme_classic(base_size = 25)+
  scale_fill_manual(values=col_3, limits = c("0","75","500"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "none")

p2 <- p1 + geom_signif(comparisons = list(c("0","500")), map_signif_level = T,
                 annotations = c("*"),
                 y_position = c(35),
                 textsize = 7,
                 size=1)

png("~/dosage_manuscript/figure_5/EdU_DMP.png", width = 9, height = 7, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()


