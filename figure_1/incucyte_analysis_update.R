
rm(list = ls())
graphics.off()

#Plots Incucyte data figures, both confluence & cytotoxicity, across doxycycline dosages in 2 cell lines (Dbt/MYCN/iP3F, Dbt/MYCN/iEV)
#needs csv files exported from Incucyte as input

library(tidyverse)

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
p1 <-  ggplot(data_sum, aes(x=time,y=mean ))+ 
  geom_line(linewidth=1.5, aes(color=dox))+
  geom_errorbar( aes(ymin=mean-se, ymax=mean+se),
                 width=1 )+
  labs(y="Confluence relative to 0 hrs", x="Time (hrs)",
        color="doxycycline\n(ng/mL)")+
  theme_classic(base_size=20)+
  scale_x_continuous(expand = c(0,0) )+ 
  scale_color_viridis_d(option = "C",end=0.8)+

png("~/dosage_manuscript/figure_1/dbt_mycn_ip3f_confluence_v2.png", width = 7, height = 5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()



#cytotoxicity plot, Extended Data figure 1A
data1 <- read_csv("~/dosage_manuscript/csv/dbt_0302024_high_red_percent.csv")[,c(1:7,11)] # data from Incucyte

data1_2 <- pivot_longer(data1, cols=c("0_hr","12_hr","24_hr","36_hr","48_hr","60_hr"), names_to = c("time"), values_to = c("Cytotoxicity") )  #select timepoints and rearrange df
data1_3 <- data1_2[ which(data1_2[,2] == "medium"),] #choose the cell density (seeded 3, low-high)
data1_3$time <- as.numeric(sapply(str_split(data1_3$time,"_"), "[[",1)) #make timepoints into numerics
data1_3$dox <- factor(data1_3$dox, levels=c("0","75","150","250","500","1000")) #factor dox dosages

data1_sum <- data1_3 %>% #summarize to get mean & std dev
  group_by(time, dox) %>%
  summarise(sd=sd(Cytotoxicity), mean=mean(Cytotoxicity), se=sd(Cytotoxicity)/sqrt(2),.groups="drop" )

#plot
p2 <- ggplot(data1_sum, aes(x=time,y=mean ))+
  geom_line(linewidth=1.5, aes(color=dox))+
  geom_errorbar( aes(ymin=mean-se, ymax=mean+se),
                 width=1)+
  labs(y="% Cytotoxicity", x="Time (hrs)",
       color="doxycycline\n(ng/mL)")+
  theme_classic(base_size=20)+
  scale_x_continuous(expand = c(0,0))+ 
  scale_y_continuous(expand = c(0,0), limits=c(0,35))+
  scale_color_viridis_d(option = "C",end=0.8)


png("~/dosage_manuscript/figure_1/dbt_mycn_ip3f_cytotoxicity_v2.png", width = 7, height = 5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()

###Dbt/MYCN/iEV

rm(list = ls())

#confluence plot, Figure 1E
data <- read_csv("~/dosage_manuscript/csv/dbt_mycn_iev_041624.csv")[,c(1:7,12)]  # data from Incucyte
data_2 <- pivot_longer(data, cols=c("0_hr","12_hr","24_hr","36_hr","48_hr","60_hr"), names_to = c("time"), values_to = c("confluence") )  #select timepoints and rearrange df
data_3 <- data_2[ which(data_2[,2] == "medium"),] #choose the cell density (seeded 3, low-high)
data_3$confluence <- data_3$confluence/100  #make initial confluence 1
data_3$time <- as.numeric(sapply(str_split(data_3$time,"_"), "[[",1))  #make timepoints into numerics
data_3$dox <- factor(data_3$dox, levels=c("0","75","150","250","500","1000")) #factor dox dosages

data_sum <- data_3 %>% #summarize to get mean & std dev
  group_by(time, dox) %>%
  summarise(sd=sd(confluence), mean=mean(confluence), se=sd(confluence)/sqrt(2),.groups="drop" )

#plot
p1 <-  ggplot(data_sum, aes(x=time,y=mean ))+ 
  geom_line(linewidth=1.5, aes(color=dox))+
  geom_errorbar( aes(ymin=mean-se, ymax=mean+se),
                 width=1 )+
  labs(y="Confluence relative to 0 hrs", x="time (hrs)",
        color="doxycycline\n(ng/mL)")+
  theme_classic(base_size=20)+
  scale_x_continuous(expand = c(0,0) )+ #, breaks=c(75,150,250,500,1000)
  scale_color_viridis_d(option = "C",end=0.8)+


png("~/dosage_manuscript/figure_1/dbt_mycn_iEV_confluence_v2.png", width = 7, height = 5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p1)
dev.off()

#cytotoxicity plot, Extended Data figure 1B
data1 <- read_csv("~/dosage_manuscript/csv/dbt_mycn_iev_cytox_050124.csv")[,c(1:7,13)] # data from Incucyte

data1_2 <- pivot_longer(data1, cols=c("0_hr","12_hr","24_hr","36_hr","48_hr","60_hr"), names_to = c("time"), values_to = c("Cytotoxicity") )  #select timepoints and rearrange df
data1_3 <- data1_2[ which(data1_2[,2] == "medium"),] #choose the cell density (seeded 3, low-high)
data1_3$time <- as.numeric(sapply(str_split(data1_3$time,"_"), "[[",1)) #make timepoints into numerics
data1_3$dox <- factor(data1_3$dox, levels=c("0","75","150","250","500","1000")) #factor dox dosages

data1_sum <- data1_3 %>% #summarize to get mean & std dev
  group_by(time, dox) %>%
  summarise(sd=sd(Cytotoxicity), mean=mean(Cytotoxicity), se=sd(Cytotoxicity)/sqrt(2),.groups="drop" )

#plot
p2 <- ggplot(data1_sum, aes(x=time,y=mean ))+ 
  geom_line(linewidth=1.5, aes(color=dox))+
  geom_errorbar( aes(ymin=mean-se, ymax=mean+se),
                 width=1)+
  labs(y="% Cytotoxicity", x="Time (hrs)",
       color="doxycycline\n(ng/mL)")+
  theme_classic(base_size=20)+
  scale_x_continuous(expand = c(0,0))+ 
  scale_y_continuous(expand = c(0,0), limits=c(0,35))+
  scale_color_viridis_d(option = "C",end=0.8)

png("~/dosage_manuscript/figure_1/dbt_mycn_iEV_cytox_v2.png", width = 7, height = 5, units = "in", res = 200, bg = "transparent", type = "cairo-png")
print(p2)
dev.off()