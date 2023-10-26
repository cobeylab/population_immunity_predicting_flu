
#if("dplyr" %in% (.packages())){
#  detach("package:dplyr", unload=TRUE) 
#  detach("package:plyr", unload=TRUE) 
#} 
#library(plyr)
library(dplyr)

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
  data_assigned$A2 = ifelse(data_assigned$clade == "A2_2", "3C.2A2", "Not 3C.2A2")
  return (data_assigned)
}

data_assigned_NA_17 = process_data(data_assigned_NA_17)
data_assigned_US_17 = process_data(data_assigned_US_17)
data_assigned_NE_17 = process_data(data_assigned_NE_17)

######################################################################################

library(ggpubr)
library(cowplot)

# NA 17

phist = gghistogram( data_assigned_NA_17,
                      x="age", fill="A2", bins=20,
                     palette = c("#C77CFF", "Black")) +
  theme(legend.title = element_blank(),
        legend.position = "top")


pdensity = ggdensity( data_assigned_NA_17,
                      x="age", color="A2", alpha=0,
                      palette = c("#C77CFF", "Black"), size=1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), position = "right") +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis") +
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") +
  ggtitle("North America") +
  coord_cartesian(clip="off")


aligned_plots = align_plots(phist, pdensity, align="hv", axis="tblr")
p_NA_17 = ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p_NA_17
#ggsave("../figure/density_histogram_NA_17.png", height=3, width=3)


# US 17

phist = gghistogram( data_assigned_US_17,
                     x="age", fill="A2", bins=20,
                     palette = c("#C77CFF", "Black")) +
  theme(legend.position="none")


pdensity = ggdensity( data_assigned_US_17,
                      x="age", color="A2", alpha=0,
                      palette = c("#C77CFF", "Black"), size=1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), position = "right") +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis") +
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") +
  ggtitle("United States")+
  coord_cartesian(clip="off")

aligned_plots = align_plots(phist, pdensity, align="hv", axis="tblr")
p_US_17 = ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p_US_17
#ggsave("../figure/density_histogram_US_17.png", height=3, width=3)


# NE 17


phist = gghistogram( data_assigned_NE_17,
                     x="age", fill="A2", bins=20,
                     palette = c("#C77CFF", "Black")) +
  theme(legend.position="none")


pdensity = ggdensity( data_assigned_NE_17,
                      x="age", color="A2", alpha=0,
                      palette = c("#C77CFF", "Black"), size=1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), position = "right") +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis") +
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") + 
  ggtitle("Northeastern US")+
  coord_cartesian(clip="off")

aligned_plots = align_plots(phist, pdensity, align="hv", axis="tblr")
p_NE_17 = ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p_NE_17
#ggsave("../figure/density_histogram_NE_17.png", height=3, width=3)


ggarrange(p_NA_17, p_US_17, p_NE_17,
          nrow=3, heights=c(7,5,5),
          labels=c("A", "B", "C")) 

ggsave("../figure/density_histogram_1718.png", height=6, width=4)

