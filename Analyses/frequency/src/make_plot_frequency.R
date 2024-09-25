#if("dplyr" %in% (.packages())){
#  detach("package:dplyr", unload=TRUE) 
#  detach("package:plyr", unload=TRUE) 
#} 
#library(plyr)
library(dplyr)
library(stringr)

########################################################################################
data_assigned_NE_16 = read.csv("../data_clade_assigned/NE_clade_assigned_season_2016.csv")
data_assigned_NE_17 = read.csv("../data_clade_assigned/NE_clade_assigned_season_2017.csv")

data_assigned_US_16 = read.csv("../data_clade_assigned/US_clade_assigned_season_2016.csv")
data_assigned_US_17 = read.csv("../data_clade_assigned/US_clade_assigned_season_2017.csv")

data_assigned_NA_16 = read.csv("../data_clade_assigned/NA_clade_assigned_season_2016.csv")
data_assigned_NA_17 = read.csv("../data_clade_assigned/NA_clade_assigned_season_2017.csv")

# MV: Later inclusion of data from 18-19 for North American (for frequency plot at finer timescale)
data_assigned_NA_18 = read.csv("../data_clade_assigned/NA_clade_assigned_season_2018.csv")



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

# Export to result folder
write.csv(tb, "../result/clade_frequencies.csv", row.names = F)

##################################################
# plot frequency and population susceptibility

library(ggplot2)
susc_x_label = c("3C.3A", "3C.2A", "3C.2A1-1", "3C.2A1-2", "3C.2A1-3",
                 "3C.2A2-1", "3C.2A2-2", "3C.2A3")

p_freq = ggplot(tb, aes(x=clade, y=freq, fill=clade, linetype = factor(region))) +
  geom_bar(color="black", width = 0.7, position=position_dodge(width=0.9), stat="identity") +
  facet_wrap(~season, nrow=2) +
  ylab("Proportion of sequences") +
  xlab("") +
  scale_x_discrete(labels=susc_x_label,
                   breaks = c("c3a", "c2a", "A1", "A1_2", "A1_3", "A2", "A2_2", "A3"),
                   limits = c("c3a", "c2a", "A1", "A1_2", "A1_3", "A2", "A2_2", "A3")) +
  scale_fill_discrete(guide="none") +
  scale_linetype_manual(name = "Region", values = c("solid", "dotted", "dashed")) +
  guides(linetype=guide_legend(override.aes=list(fill=c(NA, NA, NA))) ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=310, hjust=0, vjust=1),
      text=element_text(size=8),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.text=element_text(size=8),
        legend.position=c(0.8, 0.85),
        legend.key.size=unit(0.7, "lines"))


p_freq        

################################################################

rname = paste0("../../forecasting/result/df_susc2_0.rds")
df_susc = readRDS(rname)

df_susc$rank = factor(df_susc$rank, levels = c(1,2,3,4,5,6,7,8))


p_susc = ggplot(df_susc) +
  geom_point(aes(x=Clade, y=susceptibility, color=Clade), size=3) +
  geom_errorbar(aes(x=Clade, ymin=ci1, ymax=ci2, color=Clade), size=1.5) +
  ylab("Inferred relative susceptibility") +
  scale_x_discrete(labels = susc_x_label) +
  scale_color_discrete(guide="none") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=310, hjust=0, vjust=1),
        text=element_text(size=8),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"))


#p_susc = p_susc + theme(text = element_text(size=11),
#                        legend.position="none",
#                        axis.text.x = element_text(angle=0, size=10, hjust=0, vjust=0))

library(ggpubr)
ggarrange(p_freq, p_susc,
          nrow=2, heights= c(2,1.1),
          labels=c("A", "B"),
          label.y=c(1, 1.07),
          vjust = c(1.9, 0.1)) 

ggsave("../fig/frequency.png", height=6, width=4)
  
####################################################
#correlation test

freq = freq_NA_17$freq
freq = c(freq[1], 0, 0, freq[2], freq[3], 0, freq[4], freq[5])
susc = df_susc$susceptibility
cor.test(freq, susc, method="spearman")

#test without a2-1
cor.test(freq[-6], susc[-6], method = "spearman")

######################################
#plot frequency by age group


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


ggplot(freq_17_by_AG) +
  geom_bar(aes(x=clade, y=freq, fill=clade, linetype = factor(region)), 
           stat="identity", position=position_dodge(width=0.95), 
           width=0.65, color="black") +
  facet_wrap(~label, nrow=3) +
  scale_x_discrete(labels = susc_x_label) +
  ylab("Proportion of sequences") +
  xlab("") +
  scale_fill_discrete(guide="none") + 
  scale_linetype_manual(name = "Region", 
                        values = c("solid", "dotted", "dashed")) +
  guides(linetype=guide_legend(override.aes=list(fill=c(NA, NA, NA))) ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=300, hjust=0, vjust=1),
        text=element_text(size=8),
        axis.title = element_text(size=9),
        axis.text.y = element_text(size=8),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.text=element_text(size=8),
        legend.position=c(0.8, 0.16),
        legend.key.size=unit(0.7, "lines"))
  

ggsave("../fig/frequency_by_agegroup_17.png", height=6, width=4)

################################################################



make_frequency_table_by_ageGroup_16 = function (data_assigned_16) {
  
  data_assigned_16 = data_assigned_16[data_assigned_16$clade != "other",]
  data_assigned_16 = data_assigned_16[ !is.na(data_assigned_16$age),]
  data_assigned_16 = data_assigned_16[ data_assigned_16$age <= 90, ]
  data_assigned_16$age = ifelse(data_assigned_16$age == 0, 1, data_assigned_16$age)
  
  data_assigned_16 = transform(data_assigned_16, age_group = cut(age, breaks = c(0, 4, 16, 44, 64, 90)))
  
  
  f16 = data_assigned_16 %>%
    group_by(age_group, clade) %>%
    summarise(n = n()) %>%
    mutate(freq = n/sum(n))
  
  f16$label = f16$age_group
  f16$label = gsub("\\(0,4\\]", "1-4 years", f16$label)
  f16$label = gsub("\\(4,16\\]", "5-16 years", f16$label)
  f16$label = gsub( "\\(16,44\\]", "18-44 years", f16$label)
  f16$label = gsub( "\\(44,64\\]", "45-64 years", f16$label)
  f16$label = gsub("\\(64,90\\]", "65-90 years", f16$label)
  
  f16 = f16[,c(2,3,4,5)]

  ft16 = f16
  
  ft16$clade = as.factor(ft16$clade)
  ft16$clade = factor(ft16$clade, levels =c("c3a", "c2a", "A1", "A1_2", "A1_3", "A2", "A2_2", "A3", "other"))
  ft16$label = as.factor(ft16$label)
  ft16$label = factor(ft16$label, levels = c("1-4 years", "5-16 years", 
                                             "18-44 years", "45-64 years", "65-90 years"))
  
  return (ft16)
  
}



freq_NA_16_by_AG = make_frequency_table_by_ageGroup_16(data_assigned_NA_16)
freq_US_16_by_AG = make_frequency_table_by_ageGroup_16(data_assigned_US_16)
freq_NE_16_by_AG = make_frequency_table_by_ageGroup_16(data_assigned_NE_16)

freq_NA_16_by_AG$region = "North America"
freq_US_16_by_AG$region = "United States"
freq_NE_16_by_AG$region = "Northeastern US"


freq_16_by_AG = rbind(freq_NA_16_by_AG, freq_US_16_by_AG, freq_NE_16_by_AG)
freq_16_by_AG$region = factor(freq_16_by_AG$region, levels = c("North America", "United States", "Northeastern US"))

ggplot(freq_16_by_AG) +
  geom_bar(aes(x=clade, y=freq, fill=clade, linetype = factor(region)), 
           stat="identity", position=position_dodge(width=0.95), 
           width=0.65, color="black") +
  facet_wrap(~label, nrow=3) +
  scale_x_discrete(labels = susc_x_label) +
  ylab("Proportion of sequences") +
  xlab("") +
  scale_fill_discrete(guide="none") + 
  scale_linetype_manual(name = "Region", 
                        values = c("solid", "dotted", "dashed")) +
  guides(linetype=guide_legend(override.aes=list(fill=c(NA, NA, NA))) ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=300, hjust=0, vjust=1),
        text=element_text(size=8),
        axis.title = element_text(size=9),
        axis.text.y = element_text(size=8),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.text=element_text(size=8),
        legend.position=c(0.8, 0.16),
        legend.key.size=unit(0.7, "lines"))


#ggsave("../fig/frequency_by_agegroup_16.png", height=6, width=4)

### ======== MV: plot with frequency over time at a finer resolution ===========

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

clade_colors <- tibble(clade = c("other", "3C.3a", "3C.2a", "3C.2a1-1", "3C.2a1-2","3C.2a1-3",
                                 "3C.2a2-1", "3C.2a2-2", "3C.2a3"),
                       color = c("grey85","#F8766D", "#CD9600",
                                 "#7CAE00", "#00BE67", "#00BFC4", 
                                 "#6666FF", "#C77CFF", "#FF61CC"))

north_america_monthly_freqs <- bind_rows(data_assigned_NA_16,
                                         data_assigned_NA_17) %>%
  mutate(id = as.character(id)) %>%
  bind_rows(data_assigned_NA_18) %>%
  as_tibble() %>%
  # Keep only correctly formatted dates
  filter(str_count(collection, "-") == 2) %>%
  mutate(year = lubridate::year(collection),
         month = lubridate::month(collection),
         # Reference date for each year/month (for plotting)
         year_month_ref_date = lubridate::ymd(paste(year, month, 1, sep = '-'))) %>%
  group_by(year, month, year_month_ref_date, clade) %>%
  count() %>%
  group_by(year, month) %>%
  mutate(clade_freq = n/sum(n)) %>%
  ungroup() 

north_america_monthly_freqs <- north_america_monthly_freqs %>%
  mutate(
    clade = case_match(clade,
                       "A1" ~ "3C.2a1-1",
                       "A1_2" ~ "3C.2a1-2",
                       "A1_3" ~ "3C.2a1-3",
                       "A2" ~ "3C.2a2-1",
                       "A2_2" ~ "3C.2a2-2",
                       "A3" ~ "3C.2a3",
                       "c2a" ~ "3C.2a",
                       "c3a" ~ "3C.3a",
                       "other" ~ "other")
  )%>%
  mutate(clade = factor(clade, levels = clade_colors$clade))


# Quick test, frequencies should sum to 1
stopifnot(
  all(north_america_monthly_freqs %>% group_by(month, year) %>%
        summarise(S = sum(clade_freq)) %>%
        pull(S) %>% unique() %>%
        round(10) == 1)
)

extended_freq_plot <- north_america_monthly_freqs %>%
  # Using discrete instead of date x axis to make spacing between columns even
  ggplot(aes(x = as.character(year_month_ref_date), y = clade_freq)) +
  scale_fill_manual(values = clade_colors$color, name = 'Clade') +
  geom_col(aes(fill = clade),
           width = 0.95, just = 0) +
  scale_x_discrete(breaks = lubridate::ymd(
                          c("2016-11-01",
                            paste("2017", seq(1,11,2), "01", sep = '-'),
                            paste("2018", seq(1,11,2), "01", sep = '-'),
                            paste("2019", seq(1,9,2), "01", sep = '-')
                            )) %>%
                     as.character(),
               labels = function(x){
                 
                 case_when(
                   lubridate::month(x) == 1 & lubridate::day(x) == 1 ~
                     paste0(lubridate::month(x, label = T, abbr = T),"\n",lubridate::year(x)),
                   T ~ lubridate::month(x, label = T, abbr = T)
                 )

               }) +
  # Hacky way to draw a line between Sep and Oct 2017 given discrete scale of x axis
  geom_col(data = tibble(year_month_ref_date = c("2017-10-01","2018-10-01"),
                         clade_freq = 1.1),
           color = 'black', width = 0.05, just = 0, linetype = 1, linewidth = 1.1) +
  geom_text(data = tibble(
    year_month_ref_date = c("2017-04-01", "2018-04-01", "2019-04-01"),
    clade_freq = 1.1,
    label = c("2016-2017 season", "2017-2018 season", "2018-2019 season")),
    aes(label = label), 
    size = 2.5) +
  theme(legend.position = "right",
        axis.text = element_text(lineheight = 1.1)) +
  scale_y_continuous(breaks = seq(0,1,0.25)) +
  xlab("") +
  ylab("Frequency")

save(extended_freq_plot, file =  "../result/extended_freq_plot.RData")







