rm(list = ls())
graphics.off()

#plots colony and focus formation quantification, Figure 1

library(tidyverse)
library(broom)
library(gt)
library(ggpubr)
library(rstatix)
library(cowplot)

#read in quantifications
foci <- read_csv("~/dosage_manuscript/csv/foci_colony_quant.csv")

foci$dox <- factor(foci$dox, levels=c("0","75","150","250","500","1000"))

#run ANOVA on focus formation data
summary(aov(foci ~ dox, data = foci)) #***
tukey_hsd(aov(foci ~ dox, data = foci))

# # A tibble: 15 × 9
#    term  group1 group2 null.value estimate conf.low conf.high    p.adj p.adj.signif
#  * <chr> <chr>  <chr>       <dbl>    <dbl>    <dbl>     <dbl>    <dbl> <chr>       
#  1 dox   0      75              0    1.000  -14.7      16.7   1   e+ 0 ns          
#  2 dox   0      150             0   66.5     50.8      82.2   1.11e- 9 ****        
#  3 dox   0      250             0  113.      97.3     129.    1.94e-13 ****        
#  4 dox   0      500             0   51.0     35.3      66.7   7.64e- 8 ****        
#  5 dox   0      1000            0   15.0     -0.739    30.7   6.69e- 2 ns          
#  6 dox   75     150             0   65.5     49.8      81.2   1.42e- 9 ****        
#  7 dox   75     250             0  112.      96.3     128.    2.19e-13 ****        
#  8 dox   75     500             0   50       34.3      65.7   1.04e- 7 ****        
#  9 dox   75     1000            0   14.0     -1.74     29.7   9.82e- 2 ns          
# 10 dox   150    250             0   46.5     30.8      62.2   3.12e- 7 ****        
# 11 dox   150    500             0  -15.5    -31.2       0.239 5.5 e- 2 ns          
# 12 dox   150    1000            0  -51.5    -67.2     -35.8   6.57e- 8 ****        
# 13 dox   250    500             0  -62.0    -77.7     -46.3   3.46e- 9 ****        
# 14 dox   250    1000            0  -98.0   -114.      -82.3   1.5 e-12 ****        
# 15 dox   500    1000            0  -36      -51.7     -20.3   1.21e- 5 ****   

#set up for plotting the statistical comparisons
foci_comp <- list(c("75","150"),c("0","150"),c("500","1000"),c("150","250"),c("250","500"), c("0","250"), c("0","500"), c("75","250"),c("75","500"),
             c("150","1000"), c("250","1000"))
stars_foci <- rep("****",times=11)

#plot focus formation data
p1 <- ggplot(foci, aes(x=dox,y=foci,color=dox))+
    geom_boxplot(show.legend=F,outliers=F)+
    geom_jitter(show.legend=F, height=0)+
    geom_signif(comparisons = foci_comp, map_signif_level = TRUE, annotation=stars_foci, 
    y_position=c(80, 87, 65, 125, 132, 139,146, 153, 160, 167,174), show.legend=F, color = "black",textsize=4,vjust = 0.5)+
    labs(y="Number of foci", x="Doxcycline (ng/mL)",title="Focus formation")+
    scale_y_continuous(limits=c(-1,190))+
    theme_classic(base_size=25)+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_viridis_d(option = "C",end=0.8)

#run ANOVA on colony formation data
summary(aov(colonies ~ dox, data = foci)) #***
tukey_hsd(aov(colonies ~ dox, data = foci))

# # A tibble: 15 × 9
#    term  group1 group2 null.value estimate conf.low conf.high         p.adj p.adj.signif
#  * <chr> <chr>  <chr>       <dbl>    <dbl>    <dbl>     <dbl>         <dbl> <chr>       
#  1 dox   0      75              0     39.2   -0.508    79.0   0.0542        ns          
#  2 dox   0      150             0    -53.3  -93.0     -13.5   0.00536       **          
#  3 dox   0      250             0    -52.3  -92.0     -12.5   0.00635       **          
#  4 dox   0      500             0    -77.5 -117.      -37.7   0.0000959     ****        
#  5 dox   0      1000            0   -118.  -158.      -78.2   0.000000291   ****        
#  6 dox   75     150             0    -92.5 -132.      -52.7   0.00000963    ****        
#  7 dox   75     250             0    -91.5 -131.      -51.7   0.0000112     ****        
#  8 dox   75     500             0   -117.  -157.      -77.0   0.000000342   ****        
#  9 dox   75     1000            0   -157.  -197.     -117.    0.00000000324 ****        
# 10 dox   150    250             0      1    -38.8      40.8   1             ns          
# 11 dox   150    500             0    -24.3  -64.0      15.5   0.412         ns          
# 12 dox   150    1000            0    -64.8 -105.      -25.0   0.00077       ***         
# 13 dox   250    500             0    -25.3  -65.0      14.5   0.37          ns          
# 14 dox   250    1000            0    -65.8 -106.      -26.0   0.000652      ***         
# 15 dox   500    1000            0    -40.5  -80.3      -0.742 0.0444        *    

#set up statistical comparisons to plot
colony_comp <- list(c("0","150"),c("0","250"),c("0","500"),c("0","1000"),c("75","150"),c("75","250"),
                c("75","500"), c("75","1000"), c("500","1000"), c("250","1000"),c("150","1000") )

stars_colony <- c("**","**","****","****","****","****","****","****","*","***","***")

#plots colony formation data
p2 <- ggplot(foci, aes(x=dox,y=colonies,color=dox))+
    geom_boxplot(show.legend=F,outliers=F)+
    geom_jitter(show.legend=F, height=0)+
    geom_signif(comparisons = colony_comp, map_signif_level = TRUE, annotation=stars_colony, 
    y_position=c(202, 212, 222, 232, 242, 252,262, 272, 80, 100,110), show.legend=F, color = "black", textsize=4,vjust = 0.5)+
    labs(y="Number of colonies", x="Doxcycline (ng/mL)",title="Colony formation")+
    scale_y_continuous(limits=c(0,290))+
    theme_classic(base_size=25)+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_viridis_d(option = "C",end=0.8)

#save focus and colony formation data plots
png("~/dosage_manuscript/figure_1/dbt_foci_colony_quant_2.png", width = 11, height = 4.75, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(plot_grid(p2, p1))
dev.off()
