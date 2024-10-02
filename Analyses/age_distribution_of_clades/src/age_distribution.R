
#if("dplyr" %in% (.packages())){
#  detach("package:dplyr", unload=TRUE) 
#  detach("package:plyr", unload=TRUE) 
#} 
#library(plyr)
library(dplyr)
library(cowplot)
library(ggplot2)
theme_set(theme_cowplot())

########################################################################################

data_assigned_NE_16 = read.csv("../../frequency/data_clade_assigned/NE_clade_assigned_season_2016.csv")
data_assigned_NE_17 = read.csv("../../frequency/data_clade_assigned/NE_clade_assigned_season_2017.csv")

data_assigned_US_16 = read.csv("../../frequency/data_clade_assigned/US_clade_assigned_season_2016.csv")
data_assigned_US_17 = read.csv("../../frequency/data_clade_assigned/US_clade_assigned_season_2017.csv")

data_assigned_NA_16 = read.csv("../../frequency/data_clade_assigned/NA_clade_assigned_season_2016.csv")
data_assigned_NA_17 = read.csv("../../frequency/data_clade_assigned/NA_clade_assigned_season_2017.csv")


process_data = function(data_assigned){
  data_assigned = data_assigned[data_assigned$clade != "other",]
  data_assigned = data_assigned[ !is.na(data_assigned$age),]
  data_assigned = data_assigned[ data_assigned$age <= 90, ]
  data_assigned$age = ifelse(data_assigned$age == 0, 1, data_assigned$age)
  data_assigned$A2 = ifelse(data_assigned$clade == "A2_2", "3C.2a2", "Not 3C.2a2")
  return (data_assigned)
}

data_assigned_NA_17 = process_data(data_assigned_NA_17) %>%
  mutate(region = "North America")
data_assigned_US_17 = process_data(data_assigned_US_17) %>%
  mutate(region = 'United States')
data_assigned_NE_17 = process_data(data_assigned_NE_17) %>%
  mutate(region = "Northeastern US")


combined_data <- bind_rows(data_assigned_NA_17, data_assigned_US_17, data_assigned_NE_17) %>%
  as_tibble() %>%
  mutate(region = factor(region,levels = c("North America", "United States", "Northeastern US")))


binw = 3



# Plot with US alone (for main text)
density_histogram_1718_all_regions <- combined_data %>% 
  filter(region == 'United States') %>%
  ggplot(aes(x = age)) +
  geom_histogram(aes(y = ..count.. * 0.0001, fill = A2), alpha = 0.8, 
                 color = 'black', binwidth = binw, position="identity") +
  scale_y_continuous(name = 'Number of GISAID isolates',
                     labels = function(x){x/0.0001},
                     sec.axis = sec_axis(~ ., name = 'Density')) +
  geom_density(aes(color = A2), linewidth = 1.5, adjust = 1/2)  + 
  theme(legend.position = c(0.35,0.95),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 11),
        legend.direction = "horizontal") +
  scale_fill_manual(name = '', values = c("#C77CFF","gray")) +
  scale_color_manual(name = '', values = c("purple","black")) +
  guides(color = NULL) +
  xlab('Age')

save_plot("../figure/density_histogram_1718_US.pdf",
          density_histogram_1718_all_regions,
          base_height = 3, base_width = 5)


# Plot with all 3 regions (for the supplement)
density_histogram_1718_all_regions <- combined_data %>% 
  ggplot(aes(x = age)) +
  geom_histogram(aes(y = ..count.. * 0.0001, fill = A2), alpha = 0.8, 
                 color = 'black', binwidth = binw, position="identity") +
  scale_y_continuous(name = 'Number of GISAID isolates',
                     labels = function(x){x/0.0001},
                     sec.axis = sec_axis(~ ., name = 'Density')) +
  geom_density(aes(color = A2), linewidth = 2, adjust = 1/2)  + 
  facet_wrap("region", ncol = 1, scales = 'free_y') +
  theme(legend.position = 'top') +
  scale_fill_manual(name = '', values = c("#C77CFF","gray")) +
  scale_color_manual(name = '', values = c("purple","black")) +
  guides(color = NULL) +
  xlab('Age')

save_plot("../figure/density_histogram_1718_all_regions.pdf",
          density_histogram_1718_all_regions,
          base_height = 9, base_width = 5)





