rm(list = ls())
graphics.off()

###### runs ANOVA and plots quantifications for Western blot data across chromatin fractions in Dbt/MYCN/iP3F, Supplemental Figure 1

library(tidyverse)
library(ggsignif)
library(patchwork)

all_tmpts <- read_csv("~/dosage_manuscript/csv/western_quant_dmp_fractionations_final.csv")
all_tmpts$dosage <- factor(all_tmpts$dosage, levels = c("0","75","150","250","500","1000"))

############# determine significance and set up for plotting

#cytoplasm
subset_fraction <- all_tmpts[which(all_tmpts$fraction== "cyto"),]
temp <- aov(ratio ~ dosage, data=subset_fraction)
summary(temp)
temp %>% TukeyHSD()
cyto_comp <- list(c("0","250"),c("0","500"),c("0","1000"),c("75","1000"))
cyto_sig <- c("**","**","***","*")

#nucleoplasm
subset_fraction <- all_tmpts[which(all_tmpts$fraction== "nucleoplasm"),]
temp <- aov(ratio ~ dosage, data=subset_fraction)
summary(temp)
temp %>% TukeyHSD()
nuc_comp <- list(c("0","1000"))
nuc_sig <- c("*")

#euchromatin
subset_fraction <- all_tmpts[which(all_tmpts$fraction== "euchromatin"),]
temp <- aov(ratio ~ dosage, data=subset_fraction)
summary(temp)
temp %>% TukeyHSD()
euch_comp <- list(c("0","250"),c("0","1000"))
euch_sig <- c("*","*")

#soluble heterochromatin
subset_fraction <- all_tmpts[which(all_tmpts$fraction== "soluble_het"),]
temp <- aov(ratio ~ dosage, data=subset_fraction)
summary(temp)
temp %>% TukeyHSD()
sol_comp <- list(c("150","250"),c("150","500"),c("150","1000"),c("75","250"),c("75","500"),c("75","1000"), c("0","250"),c("0","500"),c("0","1000"))
sol_sig <- c("*","*","*","*","*","*","****","****","****")

#heterochromatin pellet
subset_fraction <- all_tmpts[which(all_tmpts$fraction== "pellet"),]
temp <- aov(ratio ~ dosage, data=subset_fraction)
summary(temp)
temp %>% TukeyHSD()
pellet_comp <- list(c("0","250"),c("0","500"),c("0","1000"))
pellet_sig <- c("*","**","**")

############# plot

fractions.v <- c("cyto","nucleoplasm","euchromatin","soluble_het","pellet")

#cyto    
fraction1 <- all_tmpts[which(all_tmpts$fraction== fractions.v[1]),]

p1 <- ggplot(fraction1, aes(y=ratio,x=dosage, group=dosage,color=dosage))+
    geom_boxplot(outliers = F, show.legend = F)+
    geom_jitter(width=0.1,show.legend = F)+
    scale_y_continuous(limits=c(NA,3.7))+
    geom_signif(comparisons = cyto_comp, map_signif_level = TRUE, 
                annotation=cyto_sig, y_position = c(2.9, 3.1,3.3,3.5) , 
                show.legend=F, color = "black")+
    labs(y="Relative PAX3::FOXO1 expression",x="Doxycycline (ng/mL)", title="Cytoplasm")+
    theme_classic(base_size=16)+
    scale_color_viridis_d(option = "C",end=0.8)+
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))

## nucleoplasm
fraction2 <- all_tmpts[which(all_tmpts$fraction== fractions.v[2]),]

p2 <- ggplot(fraction2, aes(y=ratio,x=dosage, group=dosage,color=dosage))+
    geom_boxplot(outliers = F, show.legend = F)+
    geom_jitter(width=0.1,show.legend = F)+
    geom_signif(comparisons = nuc_comp, map_signif_level = TRUE, 
                annotation=nuc_sig, 
                show.legend=F, color = "black")+
    labs(y="",x="Doxycycline (ng/mL)",title="Nucleoplasm")+
    theme_classic(base_size=16)+
    scale_color_viridis_d(option = "C",end=0.8)+
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))

## euchromatin
fraction3 <- all_tmpts[which(all_tmpts$fraction== fractions.v[3]),]

p3 <- ggplot(fraction3, aes(y=ratio,x=dosage, group=dosage,color=dosage))+
    geom_boxplot(outliers = F, show.legend = F)+
    geom_jitter(width=0.1,show.legend = F)+
    scale_y_continuous(limits=c(NA,3.4))+
    geom_signif(comparisons = euch_comp, map_signif_level = TRUE, 
                annotation=euch_sig, y_position = c(3, 3.2) , 
                show.legend=F, color = "black")+
    labs(y="",x="Doxycycline (ng/mL)", title="Euchromatin")+
    theme_classic(base_size=16)+
    scale_color_viridis_d(option = "C",end=0.8)+
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))

## sol het
fraction4 <- all_tmpts[which(all_tmpts$fraction== fractions.v[4]),]

p4 <- ggplot(fraction4, aes(y=ratio,x=dosage, group=dosage,color=dosage))+
    geom_boxplot(outliers = F, show.legend = F)+
    geom_jitter(width=0.1,show.legend = F)+
    scale_y_continuous(limits=c(NA,3.5))+
    geom_signif(comparisons = sol_comp, map_signif_level = TRUE, 
                annotation=sol_sig, y_position = c(1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3) , 
                show.legend=F, color = "black")+
    labs(y="",x="Doxycycline (ng/mL)",title="Soluble heterochromatin")+
    theme_classic(base_size=16)+
    scale_color_viridis_d(option = "C",end=0.8)+
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))

## pellet
fraction5 <- all_tmpts[which(all_tmpts$fraction== fractions.v[5]),]

p5 <- ggplot(fraction5, aes(y=ratio,x=dosage, group=dosage,color=dosage))+
    geom_boxplot(outliers = F, show.legend = F)+
    geom_jitter(width=0.1,show.legend = F)+
    scale_y_continuous(limits=c(NA,3.8))+
    geom_signif(comparisons = pellet_comp, map_signif_level = TRUE, 
                annotation=pellet_sig, y_position = c(3.2,3.4, 3.6) , 
                show.legend=F, color = "black")+
    labs(y="",x="Doxycycline (ng/mL)",title="Heterochromatin pellet")+
    theme_classic(base_size=16)+
    scale_color_viridis_d(option = "C",end=0.8)+
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))

#plot all fractions together
png(paste0("~/dosage_manuscript/figure_1/", "fractionation_quant_all", ".png"), width = 16, height = 4.5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1 + p2 + p3 + p4 + p5 + plot_layout(ncol=5))
dev.off()
