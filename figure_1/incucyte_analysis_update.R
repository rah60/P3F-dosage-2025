
rm(list = ls())
graphics.off()

#Plots Incucyte data figures across doxycycline dosages in 2 cell lines (Dbt/MYCN/iP3F, Dbt/MYCN/iEV), Figure 1
#Calculates & plots doubling time, Supplemental Figure 1
#needs csv files exported from Incucyte as input
#note that statistical analysis was carried out in figure_1/incucyte_test.R

library(tidyverse)
library(broom)
library(gt)
library(ggpubr)
library(rstatix)

###Dbt/MYCN/iP3F

#confluence plot, Figure 1D
data <- read_csv("~/dosage_manuscript/csv/dbt_032024_confluence_2.csv")[,c(1:7,11)] # data from Incucyte
data_2 <- pivot_longer(data, cols=c("0_hr","12_hr","24_hr","36_hr","48_hr","60_hr"), names_to = c("time"), values_to = c("confluence") ) #select timepoints and rearrange df
data_3 <- data_2[ which(data_2[,2] == "medium"),] #choose the cell density (seeded 3, low-high)
data_3$confluence <- data_3$confluence/100 #make initial confluence 1
data_3$time <- as.numeric(sapply(str_split(data_3$time,"_"), "[[",1)) #make timepoints into numerics
data_3$dox <- factor(data_3$dox, levels=c("0","75","150","250","500","1000")) #factor dox dosages

data_sum <- data_3 %>% #summarize to get mean & std dev
  group_by(time, dox) %>%
  summarise(sd=sd(confluence), mean=mean(confluence), se=sd(confluence)/sqrt(2),.groups="drop" )

#plot
p1 <-  ggplot()+ 
  geom_line(data_sum, mapping=aes(x=time,y=mean,color=dox ), linewidth=1.5)+
  geom_errorbar( data_sum, mapping=aes(x=time,y=mean, ymin=mean-se, ymax=mean+se), width=1 )+
  geom_text(mapping=aes(label = c("*","**","****","****","****"), x=c(12,24,36,48,60), y=c(1.6,2.1,3,4,5.1)), size=6 )+ #see incucyte_test.R for stats testing
  labs(y="Confluence relative to 0 hrs", x="Time (hrs)",
        color="doxycycline\n(ng/mL)")+
  theme_classic(base_size=20)+
  scale_x_continuous(expand = c(0,0), limits=c(0, 63) )+ 
  scale_color_viridis_d(option = "C",end=0.8)

png("~/dosage_manuscript/figure_1/dbt_mycn_ip3f_confluence_revision.png", width = 7, height = 5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()


###Dbt/MYCN/iEV

#confluence plot, Figure 1E
data <- read_csv("~/dosage_manuscript/csv/dbt_mycn_iev_041624.csv")[,c(1:7,12)]  # data from Incucyte
data_2 <- pivot_longer(data, cols=c("0_hr","12_hr","24_hr","36_hr","48_hr","60_hr"), names_to = c("time"), values_to = c("confluence") )  #select timepoints and rearrange df
data_4 <- data_2[ which(data_2[,2] == "medium"),] #choose the cell density (seeded 3, low-high)
data_4$confluence <- data_4$confluence/100  #make initial confluence 1
data_4$time <- as.numeric(sapply(str_split(data_4$time,"_"), "[[",1))  #make timepoints into numerics
data_4$dox <- factor(data_4$dox, levels=c("0","75","150","250","500","1000")) #factor dox dosages

data_sum <- data_4 %>% #summarize to get mean & std dev
  group_by(time, dox) %>%
  summarise(sd=sd(confluence), mean=mean(confluence), se=sd(confluence)/sqrt(2),.groups="drop" )

#plot
p1 <-  ggplot(data_sum, aes(x=time,y=mean ))+ 
  geom_line(linewidth=1.5, aes(color=dox))+
  geom_errorbar( aes(ymin=mean-se, ymax=mean+se),
                 width=1 )+
  labs(y="Confluence relative to 0 hrs", x="Time (hrs)",
        color="doxycycline\n(ng/mL)")+
  theme_classic(base_size=20)+
  scale_x_continuous(expand = c(0,0) )+ #, breaks=c(75,150,250,500,1000)
  scale_color_viridis_d(option = "C",end=0.8)


png("~/dosage_manuscript/figure_1/dbt_mycn_iEV_confluence_v2.png", width = 7, height = 5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()


##############################
# determine cell cycle time 

## Dbt/MYCN/iP3F

#set up vector with replicate labels
new_string <- c()
for(l in 1:(nrow(data_3)/6/3)){
    new_string <- c(new_string, rep(1, times=6),rep(2, times=6),rep(3, times=6))
}
data_3$well <- new_string
data_4$well <- new_string

#Dbt/MYCN/iP3F confluence data
ip3f <- data_3 #full time series seems give best fit

#set df to store results
doubling_time.df <- data.frame(matrix(ncol=11, nrow=18))
colnames(doubling_time.df) <- c("dosage","doubling_time","intercept","slope","std_error_intercept","std_error_slope","R_squared","df","F_stat","p_val","replicate")
num_rep <- 3
doubling_time.df$dosage <- factor(c(rep(0,times=num_rep), rep(75,times=num_rep), rep(150,times=num_rep), 
                        rep(250,times=num_rep), rep(500,times=num_rep), rep(1000,times=num_rep)))

#calculate doubling time for each dosage using linear regression & store results
row_n <- 0
for(i in 1:length(unique(ip3f$dox))){  

  subset_1 <- subset(ip3f, dox==unique(ip3f$dox)[i]) #subset dosage

  for(f in 1:num_rep){

    model1 <- lm(log(confluence) ~ time, subset(subset_1, well == f) )
    test <- summary(model1 ) #linear regression for that dosage

    #calculate doubling time
    dbt <- log(2)/test$coefficients[2]
    row_n <- row_n + 1
    doubling_time.df[ row_n, 2] <- dbt
    doubling_time.df[ row_n, 3] <- test$coefficients[1] #intercept
    doubling_time.df[ row_n, 4] <- test$coefficients[2] #slope
    doubling_time.df[ row_n, 5] <- test$coefficients[3] #std error intercept
    doubling_time.df[ row_n, 6] <- test$coefficients[4] #std error slope
    doubling_time.df[ row_n, 7] <- test$r.squared #r-squared
    doubling_time.df[ row_n,8] <- test$df[2] #degrees of freedom
    doubling_time.df[ row_n,9] <- test$fstatistic[1] #f statistic
    doubling_time.df[ row_n, 10] <- glance(test)[,"p.value"] #pulls overall p-value from summary output
    doubling_time.df[ row_n, 11] <- f

  }
}

##### repeat these calculations for Dbt/MYCN/iEV
iev <- data_4 

#store in df
doubling_time2.df <- data.frame(matrix(ncol=11, nrow=18))
colnames(doubling_time2.df) <- c("dosage","doubling_time","intercept","slope","std_error_intercept","std_error_slope","R_squared","df","F_stat","p_val","replicate")
num_rep <- 3
doubling_time2.df$dosage <- factor(c(rep(0,times=num_rep), rep(75,times=num_rep), rep(150,times=num_rep), 
                        rep(250,times=num_rep), rep(500,times=num_rep), rep(1000,times=num_rep)))


row_n <- 0
for(i in 1:length(unique(iev$dox))){  

  subset_1 <- subset(iev, dox==unique(iev$dox)[i]) #subset dosage

  for(f in 1:num_rep){

    model1 <- lm(log(confluence) ~ time, subset(subset_1, well == f) )
    test <- summary(model1 ) #linear regression for that dosage

    #calculate doubling time
    dbt <- log(2)/test$coefficients[2]
    row_n <- row_n + 1
    doubling_time2.df[ row_n, 2] <- dbt
    doubling_time2.df[ row_n, 3] <- test$coefficients[1] #intercept
    doubling_time2.df[ row_n, 4] <- test$coefficients[2] #slope
    doubling_time2.df[ row_n, 5] <- test$coefficients[3] #std error intercept
    doubling_time2.df[ row_n, 6] <- test$coefficients[4] #std error slope
    doubling_time2.df[ row_n, 7] <- test$r.squared #r-squared
    doubling_time2.df[ row_n,8] <- test$df[2] #degrees of freedom
    doubling_time2.df[ row_n,9] <- test$fstatistic[1] #f statistic
    doubling_time2.df[ row_n, 10] <- glance(test)[,"p.value"] #pulls overall p-value from summary output
    doubling_time2.df[ row_n, 11] <- f

  }
}

#run ANOVA on Dbt/MYCN/iP3F doubling times
aov_dmp <- aov(doubling_time ~ dosage, data = doubling_time.df) %>% tukey_hsd()

comp_dmp <- list( c("75","150"), c("0","150"),c("75","250"),c("0","250"),c("150","500"),
              c("75","500"),c("0","500"), c("150","1000"),  c("75","1000"), c("0", "1000"))

signif_dmp <- c("*","**","**","****","*",
                "****","****","**","****","****")

#plot doubling time for Dbt/MYCN/iP3F
p3 <- ggplot(doubling_time.df)+
  geom_boxplot(aes(x=dosage, y=doubling_time, fill=dosage))+
  scale_fill_viridis_d(option = "C",end=0.8)+
  scale_y_continuous(limits=c(20,72),expand=c(0,0))+
  theme_classic(base_size=30)+
  geom_signif(aes(x=dosage, y=doubling_time), comparisons = comp_dmp, map_signif_level = TRUE, annotation=signif_dmp, 
  y_position=c(45, 47, 52, 54, 57,59, 61,64,66,68), show.legend=F, color = "black")+
  theme(legend.position = "none")+
  labs(y="Doubling time (hrs)", x="Doxycycline\n(ng/mL)")

png("~/dosage_manuscript/figure_1/dbt_mycn_iP3F_dt.png", width = 7, height = 7, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p3)
dev.off()


#run ANOVA on Dbt/MYCN/iEV doubling times
aov_iev <- aov(doubling_time ~ dosage, data = doubling_time2.df) #all ns

#plot doubling time for Dbt/MYCN/iEV
p4 <- ggplot(doubling_time2.df)+
  geom_boxplot(aes(x=dosage, y=doubling_time, fill=dosage))+
  scale_fill_viridis_d(option = "C",end=0.8)+
  theme_classic(base_size=30)+
  scale_y_continuous(limits=c(15,25),expand=c(0,0))+
  theme(legend.position = "none")+
  labs(y="Doubling time (hrs)", x="Doxycycline\n(ng/mL)")

png("~/dosage_manuscript/figure_1/dbt_mycn_iEV_dt.png", width = 7, height = 7, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p4)
dev.off()