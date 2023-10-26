
#if("dplyr" %in% (.packages())){
#  detach("package:dplyr", unload=TRUE) 
#  detach("package:plyr", unload=TRUE) 
#} 
#library(plyr)
library(dplyr)

########################################################################################
# Read sequence meta data assigned to each reference strain

data_assigned_NE_16 = read.csv("../data_clade_assigned/NE_clade_assigned_season_2016.csv")
data_assigned_NE_17 = read.csv("../data_clade_assigned/NE_clade_assigned_season_2017.csv")

data_assigned_US_16 = read.csv("../data_clade_assigned/US_clade_assigned_season_2016.csv")
data_assigned_US_17 = read.csv("../data_clade_assigned/US_clade_assigned_season_2017.csv")

data_assigned_NA_16 = read.csv("../data_clade_assigned/NA_clade_assigned_season_2016.csv")
data_assigned_NA_17 = read.csv("../data_clade_assigned/NA_clade_assigned_season_2017.csv")

make_frequency_table = function(data_assigned){
  data_assigned = data_assigned[data_assigned$clade != "other", ]
  data_assigned$clade = factor(data_assigned$clade,
                      levels = c("c3a", "c2a", "A1", "A1_2", "A1_3", "A2", "A2_2", "A3"))
  freq_tb = data_assigned %>% 
            group_by(clade) %>%
            summarise(count = length(clade))
  freq_tb$freq = freq_tb$count / sum(freq_tb$count)
  
  return(freq_tb)
}

freq_NE_16 = make_frequency_table(data_assigned_NE_16)
freq_NE_17 = make_frequency_table(data_assigned_NE_17)

freq_US_16 = make_frequency_table(data_assigned_US_16)
freq_US_17 = make_frequency_table(data_assigned_US_17)

freq_NA_16 = make_frequency_table(data_assigned_NA_16)
freq_NA_17 = make_frequency_table(data_assigned_NA_17)


#Assign season

freq_NE_16$season = "2016-17"
freq_US_16$season = "2016-17"
freq_NA_16$season = "2016-17"

freq_NE_17$season = "2017-18"
freq_US_17$season = "2017-18"
freq_NA_17$season = "2017-18"

#Assign region

freq_NE_16$region = "Northeastern US"
freq_US_16$region = "United States"
freq_NA_16$region = "North America"

freq_NE_17$region = "Northeastern US"
freq_US_17$region = "United States"
freq_NA_17$region = "North America"

#bind all together
tb = rbind(freq_NE_16, freq_NE_17,
           freq_US_16, freq_US_17,
           freq_NA_16, freq_NA_17)

tb$region = factor(tb$region, levels = c("North America", "United States", "Northeastern US"))


################################################################
# read inferred relative susceptibility using cutoff 40

rname = paste0("../../forecasting/result/df_susc2_0.rds")
df_susc = readRDS(rname)

df_susc$rank = factor(df_susc$rank, levels = c(1,2,3,4,5,6,7,8))



####################################################
# Correlation between 2017-18 frequencies and inferred relative susceptibility

#correlation test using North America
freq = freq_NA_17$freq
freq = c(freq[1], 0, 0, freq[2], freq[3], 0, freq[4], freq[5])
susc = df_susc$susceptibility
cor.test(freq, susc, method="spearman")
#test without a2-1
cor.test(freq[-6], susc[-6], method = "spearman")



#correlation test using United States
freq = freq_US_17$freq
freq = c(freq[1], 0, 0, freq[2], freq[3], 0, freq[4], freq[5])
susc = df_susc$susceptibility
cor.test(freq, susc, method="spearman")
#test without a2-1
cor.test(freq[-6], susc[-6], method = "spearman")


#correlation test using United States
freq = freq_NE_17$freq
freq = c(freq[1], 0, 0, freq[2], freq[3], 0, freq[4], freq[5])
susc = df_susc$susceptibility
cor.test(freq, susc, method="spearman")
#test without a2-1
cor.test(freq[-6], susc[-6], method = "spearman")

######################################


make_frequency_table_by_ageGroup_17 = function (data_assigned_17) {
  
  data_assigned_17 = data_assigned_17[data_assigned_17$clade != "other",]
  data_assigned_17 = data_assigned_17[ !is.na(data_assigned_17$age),]
  data_assigned_17 = data_assigned_17[ data_assigned_17$age <= 90, ]
  data_assigned_17$age = ifelse(data_assigned_17$age == 0, 1, data_assigned_17$age)
  
  data_assigned_17 = transform(data_assigned_17, age_group = cut(age, breaks = c(0, 4, 17, 44, 64, 90)))
  
  
  f17 = data_assigned_17 %>%
    group_by(age_group, clade) %>%
    summarise(n = n()) %>%
    mutate(freq = n/sum(n))
  
  f17$label = f17$age_group
  f17$label = gsub("\\(0,4\\]", "1-4 years", f17$label)
  f17$label = gsub("\\(4,17\\]", "5-17 years", f17$label)
  f17$label = gsub( "\\(17,44\\]", "18-44 years", f17$label)
  f17$label = gsub( "\\(44,64\\]", "45-64 years", f17$label)
  f17$label = gsub("\\(64,90\\]", "65-90 years", f17$label)
  
  f17 = f17[,c(2,3,4,5)]
  t17 = tibble(clade = rep(c("c2a", "A1", "A2"),5 ),
               n = rep(0, 15),
               freq = rep(0, 15),
               label = c("1-4 years", "1-4 years", "1-4 years",
                         "5-17 years", "5-17 years","5-17 years",
                         "18-44 years","18-44 years","18-44 years",
                         "45-64 years","45-64 years","45-64 years",
                         "65-90 years","65-90 years","65-90 years"))
  ft17 = rbind(f17, t17)
  
  ft17$clade = as.factor(ft17$clade)
  ft17$clade = factor(ft17$clade, levels =c("c3a", "c2a", "A1", "A1_2", "A1_3", "A2", "A2_2", "A3", "other"))
  ft17$label = as.factor(ft17$label)
  ft17$label = factor(ft17$label, levels = c("1-4 years", "5-17 years", 
                                             "18-44 years", "45-64 years", "65-90 years"))
  
  return (ft17)
  
}

freq_NA_17_by_AG = make_frequency_table_by_ageGroup_17(data_assigned_NA_17)
freq_US_17_by_AG = make_frequency_table_by_ageGroup_17(data_assigned_US_17)
freq_NE_17_by_AG = make_frequency_table_by_ageGroup_17(data_assigned_NE_17)

freq_NA_17_by_AG$region = "North America"
freq_US_17_by_AG$region = "United States"
freq_NE_17_by_AG$region = "Northeastern US"


freq_17_by_AG = rbind(freq_NA_17_by_AG, freq_US_17_by_AG, freq_NE_17_by_AG)
freq_17_by_AG$region = factor(freq_17_by_AG$region, levels = c("North America", "United States", "Northeastern US"))

###############################################################
# Chi-square test to test if proportion of A2-2 is lower in children


#make data frame of 2*2 proportion

ag_data = data.frame(freq_17_by_AG)

#1-4 vs adults
ag_region = ag_data[ag_data$region == "North America", ]
#ag_region = ag_data[ag_data$region == "United States", ]
#ag_region = ag_data[ag_data$region == "Northeastern US", ]


young_a2 = ag_region[ag_region$clade == "A2_2" & ag_region$label == "1-4 years", ]$n
young_others = ag_region[ag_region$clade != "A2_2" & ag_region$label == "1-4 years", ]$n
young_others = sum(young_others)

adults_a2 = ag_region[ag_region$clade == "A2_2" & 
                       (ag_region$label != "1-4 years" & ag_region$label != "5-17 years"), ]$n
adults_a2 = sum(adults_a2)
adults_others = ag_region[ag_region$clade != "A2_2" & 
                                          (ag_region$label != "1-4 years" & ag_region$label != "5-17 years"), ]$n
adults_others = sum(adults_others)

young = c(young_a2, young_others)
adults = c(adults_a2, adults_others)
ct = data.frame(young, adults)
row.names(ct) = c("A2", "Others")

chisq_group1 = chisq.test(ct)
chisq_group1$statistic
chisq_group1$p.vadlue
df_chisq_group1 = data.frame("1-4 years vs adults", chisq_group1$statistic, chisq_group1$p.value)
colnames(df_chisq_group1) = c("Age_group", "stat", "p")

################################################################

#5-17 vs adults
ag_region = ag_data[ag_data$region == "North America", ]
#ag_region = ag_data[ag_data$region == "United States", ]
#ag_region = ag_data[ag_data$region == "Northeastern US", ]


young_a2 = ag_region[ag_region$clade == "A2_2" & ag_region$label == "5-17 years", ]$n
young_others = ag_region[ag_region$clade != "A2_2" & ag_region$label == "5-17 years", ]$n
young_others = sum(young_others)

adults_a2 = ag_region[ag_region$clade == "A2_2" & 
                        (ag_region$label != "1-4 years" & ag_region$label != "5-17 years"), ]$n
adults_a2 = sum(adults_a2)
adults_others = ag_region[ag_region$clade != "A2_2" & 
                            (ag_region$label != "1-4 years" & ag_region$label != "5-17 years"), ]$n
adults_others = sum(adults_others)

young = c(young_a2, young_others)
adults = c(adults_a2, adults_others)
ct = data.frame(young, adults)
row.names(ct) = c("A2", "Others")

chisq.test(ct)
chisq_group2 = chisq.test(ct)
chisq_group2$statistic
chisq_group2$p.vadlue
df_chisq_group2 = data.frame("5-17 years vs adults", chisq_group2$statistic, chisq_group2$p.value)
colnames(df_chisq_group2) = c("Age_group", "stat", "p")

df_chisq = rbind(df_chisq_group1, df_chisq_group2)
write.csv(df_chisq, "../stat_result/chi_square_a2_2_of_children.csv")

