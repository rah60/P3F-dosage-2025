
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

TukeyHSD(res.aov) #0,500 0.000906; 0,75 0.0283; 75,500 0.0254

#significance testing for edu-
data_minus <- data[which(data$edu_status=="minus"),]
res.aov <- aov(percent ~ doxycycline, data=data_minus)
summary(res.aov)

TukeyHSD(res.aov) #0,500 0.000906; 0,75 0.0283; 75,500 0.0254

data_sum <- data[which(data$edu_status == "plus"),]

p1 <- ggplot(data_sum, aes(x=factor(doxycycline, c("0","75","500")), y=percent, color=doxycycline ))+ 
  geom_boxplot(width=0.7, outliers=F )+
  labs(x="PAX3::FOXO1 induction level\n(doxycycline dosage, ng/mL)", y="Percent S phase\ncells in population", fill="")+ #change cell type label
  geom_jitter(height=0, width=0.1,size=2)+
  scale_y_continuous(expand = c(0,0), limits=c(0,50))+
  theme_classic(base_size = 35)+
  scale_color_manual(values=col_3, limits = c("0","75","500"))+
  theme(plot.title=element_text(hjust=0.5), legend.position = "none")

p2 <- p1 + geom_signif(comparisons = list(c("0","500"),c("0","75"),c("75","500")), map_signif_level = T,
                 annotations = c("***","*","*"),
                 y_position = c(38,43, 46),
                 textsize = 7,
                 size=1,
                 color="black")

png("~/dosage_manuscript/figure_5/EdU_DMP_updated.png", width = 9, height = 7, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()


