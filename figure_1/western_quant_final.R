rm(list = ls())
graphics.off()

#runs ANOVA & plots Western blot quantifications in Dbt/MYCN/iP3F WCE after 8 hr, 24 hr, or 2 week doxycycline inductions, Figure 1

library(tidyverse)
library(ggsignif)
library(patchwork)

all_tmpts <- read_csv("~/dosage_manuscript/csv/western_quant_final.csv")[1:54,]

eight_hr <- all_tmpts[which(all_tmpts$timepoint == "8_hr"),]
eight_hr$dosage <- factor(eight_hr$dosage, levels = c("0","75","150","250","500","1000"))

aov(ratio ~ dosage, data=eight_hr) %>% TukeyHSD()
my_comparisons <- list(c("0","150"),c("0","250") ,c("0","500"), c("0","1000"),c("75","250"),c("75","500"),c("75","1000"),c("150","500"),c("150","1000"))

p1 <- ggplot(eight_hr, aes(y=ratio,x=dosage, group=dosage,color=dosage))+
  geom_boxplot(outliers = F, show.legend = F, lwd=1)+
  geom_jitter(width=0.1,show.legend = F)+
  scale_y_continuous(limits=c(NA,3.2))+
  geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, 
              annotation=c("*","***","****","****","**","***","****","*","**"), y_position=c(1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3), 
              show.legend=F, color = "black", textsize=5, size=1)+
  labs(y="Relative PAX3::FOXO1 expression",x="Doxycycline (ng/mL)",title="8 hour induction")+
  theme_classic(base_size=25)+
  scale_color_viridis_d(option = "C",end=0.8)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size = 20) )


twentyfour_hr <- all_tmpts[which(all_tmpts$timepoint == "24_hr"),]
twentyfour_hr$dosage <- factor(twentyfour_hr$dosage, levels = c("0","75","150","250","500","1000"))

aov(ratio ~ dosage, data=twentyfour_hr) %>% TukeyHSD()
my_comparisons <- list(c("0","1000"))

p2 <- ggplot(twentyfour_hr, aes(y=ratio,x=dosage, group=dosage, color=dosage))+
  geom_boxplot(outliers = F,show.legend = F,lwd=1)+
  scale_y_continuous(limits=c(NA,3))+
  geom_jitter(width=0.1,show.legend = F)+
  geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, 
              annotation=c("*"), y_position=c(2.7), 
              show.legend=F, color = "black", textsize=5, size=1)+
  scale_color_viridis_d(option = "C",end=0.8)+
  labs(y="",x="Doxycycline (ng/mL)",title="24 hour induction")+
  theme_classic(base_size=25)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size = 20))


two_week <- all_tmpts[which(all_tmpts$timepoint == "2_wk"),]
two_week$dosage <- factor(two_week$dosage, levels = c("0","75","150","250","500","1000"))

aov(ratio ~ dosage, data=two_week) %>% TukeyHSD()
my_comparisons <- list(c("0","250"),c("0","500"),c("0","1000"))

p3 <- ggplot(two_week, aes(y=ratio,x=dosage, group=dosage,color=dosage))+
  geom_boxplot(outliers = F,show.legend = F, lwd=1)+
  geom_jitter(width=0.1,show.legend = F)+
  scale_y_continuous(limits=c(NA,2))+
  geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, 
              annotation=c("**","**","**"), y_position=c(1.4,1.6,1.8), 
              show.legend=F, color = "black", textsize=5, size=1)+
  scale_color_viridis_d(option = "C",end=0.8)+
  labs(y="",x="Doxycycline (ng/mL)",title="2 week induction")+
  theme_classic(base_size=25)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size = 20))



png("~/dosage_manuscript/figure_1/western_quant_alltimepts_3.png", width = 14, height = 6, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1 + p2 + p3)
dev.off()
