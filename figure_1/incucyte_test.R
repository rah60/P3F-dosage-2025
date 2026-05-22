rm(list = ls())
graphics.off()

#Runs stats on Dbt/MYCN/iP3F and Dbt/MYCN/iEV confluence data across dox dosages, Figure 1
#outputs pairwise comparisons included in Supplemental Table 4
#needs csv files exported from Incucyte as input

library(tidyverse)
library(broom)
library(gt)
library(ggpubr)
library(lme4)
library(rstatix)

#iP3F
data <- read_csv("~/dosage_manuscript/csv/dbt_032024_confluence_2.csv")[,c(1:7,11)] # data from Incucyte
data_2 <- pivot_longer(data, cols=c("0_hr","12_hr","24_hr","36_hr","48_hr","60_hr"), names_to = c("time"), values_to = c("confluence") ) #select timepoints and rearrange df
data_3 <- data_2[ which(data_2[,2] == "medium"),] #choose the cell density (seeded 3, low-high)
data_3$confluence <- data_3$confluence/100 #make initial confluence 1
data_3$time <- as.numeric(sapply(str_split(data_3$time,"_"), "[[",1)) #make timepoints into numerics
data_3$dox <- factor(data_3$dox, levels=c("0","75","150","250","500","1000")) #factor dox dosages

#make vector for replicate labels
new_string <- c()
for(l in 1:(nrow(data_3)/6)){
    new_string <- c(new_string, rep(l, times=6))
}

#need to give ids to the individual wells over time
data_3$well <- new_string
data_3$well <- factor(data_3$well)
data_3$time <- factor(data_3$time)

##### run tests to ensure I can run ANOVA
data_3 %>% group_by(dox, time) %>% identify_outliers(confluence)
#none

#note, have to exclude time 0 because all the same, normalized to 1
data_3[-which(data_3$time ==0),] %>% group_by(dox, time) %>% shapiro_test(confluence)
#norm. dist.

ggqqplot(data_3, "confluence", ggtheme = theme_bw()) +
  facet_grid(dox ~ time, labeller = "label_both")
#fine

#repeated measures ANOVA
res.aov <- anova_test(
  data = data_3[-which(data_3$time ==0),], dv = confluence, wid=well,
  within =  time,
  between = dox,
  )
get_anova_table(res.aov)

# ANOVA Table (type II tests)

#     Effect  DFn   DFd       F        p p<.05   ges
# 1      dox 5.00 12.00  42.300 3.27e-07     * 0.921
# 2     time 1.26 15.08 827.032 4.04e-15     * 0.959
# 3 dox:time 6.28 15.08  46.061 6.47e-09     * 0.868

#post-hoc tests

#effect of dox 
data_3[-which(data_3$time ==0),] %>%
  group_by(time) %>%
  anova_test(dv = confluence, wid = well, between = dox) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")

#effect of dox is significant at each time point
# # A tibble: 5 × 9
#   time  Effect   DFn   DFd     F            p `p<.05`   ges       p.adj
#   <fct> <chr>  <dbl> <dbl> <dbl>        <dbl> <chr>   <dbl>       <dbl>
# 1 12    dox        5    12  5.54 0.007        *       0.698 0.035       *
# 2 24    dox        5    12  9.77 0.000655     *       0.803 0.00328    **
# 3 36    dox        5    12 28.9  0.0000027    *       0.923 0.0000135  ****
# 4 48    dox        5    12 62.3  0.0000000364 *       0.963 0.000000182 ****
# 5 60    dox        5    12 46.1  0.000000201  *       0.951 0.00000101 ****

pairwise_comp <- data_3[-which(data_3$time ==0),] %>%
  group_by(time) %>%
  pairwise_t_test(
    confluence ~ dox, paired = F,
    p.adjust.method = "bonferroni"
    )

#save for supplemental table
write_csv(pairwise_comp[,c(1,3:ncol(pairwise_comp))], "~/dosage_manuscript/figure_1/incucyte_pairwise_comp.csv")

#iEV
data <- read_csv("~/dosage_manuscript/csv/dbt_mycn_iev_041624.csv")[,c(1:7,12)]  # data from Incucyte
data_2 <- pivot_longer(data, cols=c("0_hr","12_hr","24_hr","36_hr","48_hr","60_hr"), names_to = c("time"), values_to = c("confluence") )  #select timepoints and rearrange df
data_4 <- data_2[ which(data_2[,2] == "medium"),] #choose the cell density (seeded 3, low-high)
data_4$confluence <- data_4$confluence/100  #make initial confluence 1
data_4$time <- as.numeric(sapply(str_split(data_4$time,"_"), "[[",1))  #make timepoints into numerics
data_4$dox <- factor(data_4$dox, levels=c("0","75","150","250","500","1000")) #factor dox dosages

#make vector for replicate labels
new_string <- c()
for(l in 1:(nrow(data_4)/6)){
    new_string <- c(new_string, rep(l, times=6))
}

#need to give ids to the individual wells over time
data_4$well <- new_string
data_4$well <- factor(data_4$well)
data_4$time <- factor(data_4$time)

data_4 %>% group_by(dox, time) %>% identify_outliers(confluence)
#none

#note, have to exclude time 0 because all the same, normalized to 1
data_4[-which(data_4$time ==0),] %>% group_by(dox, time) %>% shapiro_test(confluence)
#norm. dist.

ggqqplot(data_4, "confluence", ggtheme = theme_bw()) +
  facet_grid(dox ~ time, labeller = "label_both")
#fine

#repeated measures ANOVA
res.aov <- anova_test(
  data = data_4[-which(data_4$time ==0),], dv = confluence, wid=well,
  within =  time,
  between = dox,
  )
get_anova_table(res.aov)

#only time is significant, not dox
# ANOVA Table (type II tests)

#     Effect  DFn   DFd        F        p p<.05   ges
# 1      dox 5.00 12.00    1.348 3.10e-01       0.241
# 2     time 1.37 16.47 5816.845 2.93e-23     * 0.995
# 3 dox:time 6.86 16.47    1.321 3.01e-01       0.193

#one-way to look at dox at each timept confirms this
data_4[-which(data_4$time ==0),] %>%
  group_by(time) %>%
  anova_test(dv = confluence, wid = well, between = dox) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")

