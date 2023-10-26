library(ggplot2)
library(ggridges)



ag_brk = c(0,4,17,44,64,90)
sq = seq(0, 9, 1)
lb = 2^(sq) * 10


##################################################################################
source("../../../Data/script/load_data/load_data.R")

#df = transform(df, Age_group = cut(Age, breaks = ag_brk))

levels(df$Test_virus) = c("3C.3A", "3C.2A", "3C.2A1-1", "3C.2A1-2", "3C.2A1-3",
                          "3C.2A2-1", "3C.2A2-2", "3C.2A3")
levels(df$Age_group) = c("1-4 years", "5-17 years", "18-44 years",
                         "45-64 years", "65-90 years")

ggplot(df, aes(x=Titer, y=Age_group, height=stat(density), fill=Test_virus)) +
  geom_density_ridges(stat="binline", bins=10, scale=0.95, draw_baseline=F) +
  facet_grid(~Test_virus) +
  scale_x_continuous(limits = c(-1, 9), breaks=c(0,2,4,6,8), 
                     labels=c(10, 40, 160, 640, 2560)) +
  scale_y_discrete(name = "") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=300, vjust=0.7, hjust=0, size=7.5)) 
  
ggsave("../fig/titer_distribution_ha.png", height=5, width=6.5)
ggsave("../fig/titer_distribution_ha.pdf", height=5, width=6.5)
ggsave("../fig/titer_distribution_ha.tiff", height=5, width=6.5)


#####################################################################################

source("../../../Data/script/load_data/load_na.R")



#levels(df_na$Test_virus) = c("3C.3A", "3C.2A_ha", "A1-1", "A1-2", "A1-3",
#                          "A2-1", "A2-2_ha", "A3", "A2-2", "3C.2A")
#levels(df_na$Age_group) = c("1-4 years", "5-17 years", "18-44 years",
#                         "45-64 years", "65-90 years")

#df_na$Test_virus = factor(df_na$Test_virus, levels = c("3C.2A", "A2-2"))

#df$Test_virus = ifelse(df$Test_virus == "A2", "3C.2A2-2 (NA)", "3C.2A (NA)")

levels(df_na$Test_virus) = c("3C.2A (NA)", "3C.2A2-2 (NA)")
levels(df_na$Age_group) = c("1-4 years", "5-17 years", "18-44 years",
                         "45-64 years", "65-90 years")


ggplot(df_na, aes(x=Titer, y=Age_group, height=stat(density), fill=Test_virus)) +
  geom_density_ridges(stat="binline", bins=10, scale=0.95, draw_baseline=F) +
  facet_grid(~Test_virus) +
  scale_x_continuous(limits = c(-1, 9), breaks=c(0,2,4,6,8), 
                     labels=c(10, 40, 160, 640, 2560)) +
  scale_y_discrete(name = "") +
  scale_fill_manual(values=c("#CD9600", "#C77CFF")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=300, vjust=0.7, hjust=0, size=7.5)) 


ggsave("../fig/titer_distribution_na.png", height=5, width=3)
ggsave("../fig/titer_distribution_na.pdf", height=5, width=3)
ggsave("../fig/titer_distribution_na.tiff", height=5, width=3)












