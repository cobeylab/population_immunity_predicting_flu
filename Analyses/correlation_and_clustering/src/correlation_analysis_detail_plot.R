library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
source("analysis_util.R")

ag_brk = c(0, 4, 17, 44, 64, 90)

##############################################################
# Data: ag_sera (age group + HA titers)
################################################################


source("../../../Data/script/load_data/load_data.R")

# levels for viruses used for correlation analysis

vlevels = c("X3C3.A", "X3C2.A", 
            "N171K",  "N121K_N171K", "N121K_T135K_N171K",
            "T131K_R142K", "T131K_R142K_R261Q",
            "N121K_S144K")


sera0 = sera

sera0 = transform(sera0, Age_group = cut(Age, breaks = c(0, 4, 17, 44, 64, 90)))

sera0$label = sera0$Age_group
sera0$label = gsub("\\(0,4\\]", "1-4", sera0$label)
sera0$label = gsub("\\(4,17\\]", "5-17", sera0$label)
sera0$label = gsub( "\\(17,44\\]", "18-44", sera0$label)
sera0$label = gsub( "\\(44,64\\]", "45-64", sera0$label)
sera0$label = gsub("\\(64,90\\]", "65-90", sera0$label)

sera0$label = factor(sera0$label, levels = c("1-4", "5-17", "18-44", "45-64", "65-90"))

# make data frame for analysis
# age group + titers

ag_sera = sera0[, c("Age_group", vlevels)]



################################################################

#
rs = rowSums(ag_sera[ , c(2:ncol(ag_sera))], na.rm = T)
ag_sera = ag_sera[rs != 0, ]
ag_sera$ID = row.names(ag_sera)

#colnames(ag_sera) = c("Age_group", "X3C.3A", "X3C.2A", 
#                      "A1_1", "A1_2", "A1_3",
#                      "A2_1", "A2_2", 
#                      "A3", "ID")
##################################################################



cor_theme = theme_bw() + theme(
  axis.text.x = element_text(angle=290, vjust=0.7, hjust=0, size=7.5),
  axis.text.y = element_text(size=8),
  axis.title = element_text(size=8),
  strip.text = element_text(size=8),
  legend.title=element_text(size=8),
  legend.text = element_text(size=8),
  legend.key.size = unit(0.35, "cm"),
  legend.margin = unit(-0.2, "cm"),
  #panel.grid = element_blank(),
  strip.background = element_rect(fill = "white", colour = "white"))


############################################################################
#In all age group, 3C.3A has low correlation with other titers
#This can be driven by 3C.3A vs 3C.2A epitope targeting, 
#especially people who recognize 3C.3A. 

melt_3c = melt(ag_sera, id = c("Age_group", "X3C3.A", "ID"))
melt_3c$variable =  gsub("X3C2.A", "3C.2A", melt_3c$variable)
melt_3c$variable =  gsub("N121K_T135K_N171K", "3C.2A1-3", melt_3c$variable)
melt_3c$variable =  gsub("N121K_N171K", "3C.2A1-2", melt_3c$variable)
melt_3c$variable =  gsub("N171K", "3C.2A1-1", melt_3c$variable)
melt_3c$variable =  gsub("T131K_R142K_R261Q", "3C.2A2-2", melt_3c$variable)
melt_3c$variable =  gsub("T131K_R142K", "3C.2A2-1", melt_3c$variable)
melt_3c$variable =  gsub("N121K_S144K", "3C.2A3", melt_3c$variable)


ggplot(melt_3c, aes(x=X3C3.A, y=value)) +
  geom_jitter(width=0.1, height=0.1, alpha=0.25) +
  facet_wrap(~variable) +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0,10),
                     labels=c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0,10),
                     labels=c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)) +
  xlab("3C.3A titer") + ylab("Titer") +
  cor_theme +
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y=10, label.x=5, aes(label = ..eq.label..), size=2) +
  stat_regline_equation(label.y=9, label.x=5, aes(label = ..rr.label..), size=2)

ggsave("../figure/3c3a_titers_x_other_titers_all_age_groups.png", height=5.5, width=5)
ggsave("../figure/3c3a_titers_x_other_titers_all_age_groups.pdf", height=5.5, width=5)
ggsave("../figure/3c3a_titers_x_other_titers_all_age_groups.tiff", height=5.5, width=5)


#In 45-90, A1-1, A1-2, and A3 have high correlation to each other
#However, their titers are not necessarily higher than titers to other viruses
# Some people target open epitope B at 3C3A (high 3C.3A titers), 
# some people target disrupted epitope by A2 mutations (131T, 142R, 261R, low A2 titers), 
# some people target 135T (high A1-3 titiers)
# But most of these people target epitope conserved among A1-1, A1-2, A3, 

####### A1-1
melt_a11 = melt(ag_sera, id=c("Age_group", "N171K", "ID"))
melt_a11_4590 = melt_a11[melt_a11$Age_group == "(44,64]" |  melt_a11$Age_group == "(64,90]", ]
melt_a11_4590$variable =  gsub("X3C2.A", "3C.2A", melt_a11_4590$variable)
melt_a11_4590$variable =  gsub("N121K_T135K_N171K", "3C.2A1-3", melt_a11_4590$variable)
melt_a11_4590$variable =  gsub("N121K_N171K", "3C.2A1-2", melt_a11_4590$variable)
melt_a11_4590$variable =  gsub("X3C3.A", "3C.3A", melt_a11_4590$variable)
melt_a11_4590$variable =  gsub("T131K_R142K_R261Q", "3C.2A2-2", melt_a11_4590$variable)
melt_a11_4590$variable =  gsub("T131K_R142K", "3C.2A2-1", melt_a11_4590$variable)
melt_a11_4590$variable =  gsub("N121K_S144K", "3C.2A3", melt_a11_4590$variable)
melt_a11_4590$variable = factor(melt_a11_4590$variable, 
                                levels = c("3C.3A", "3C.2A", "3C.2A1-1", "3C.2A1-2", "3C.2A1-3",
                                           "3C.2A2-1", "3C.2A2-2", "3C.2A3"))

ggplot(melt_a11_4590, aes(x=N171K, y=value)) +
  geom_jitter(width=0.1, height=0.1, alpha=0.25) +
  facet_wrap(~variable)+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0,10),
                     labels=c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0,10),
                     labels=c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)) +
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y=10, label.x=5, aes(label = ..eq.label..), size=2) +
  stat_regline_equation(label.y=9, label.x=5, aes(label = ..rr.label..), size=2) +
  xlab("3C.2A1-1 titer") + ylab("Titer") +
  cor_theme


ggsave("../figure/A1-1_titers_x_other_titers_4590.png", height=5.5, width=5)
ggsave("../figure/A1-1_titers_x_other_titers_4590.pdf", height=5.5, width=5)
ggsave("../figure/A1-1_titers_x_other_titers_4590.tiff", height=5.5, width=5)


####### A1-2
melt_a12 = melt(ag_sera, id=c("Age_group", "N121K_N171K", "ID"))
melt_a12_4590 = melt_a12[melt_a12$Age_group == "(44,64]" | melt_a12$Age_group == "(64,90]", ]

melt_a12_4590$variable =  gsub("X3C2.A", "3C.2A", melt_a12_4590$variable)
melt_a12_4590$variable =  gsub("N121K_T135K_N171K", "3C.2A1-3", melt_a12_4590$variable)
melt_a12_4590$variable =  gsub("N171K", "3C.2A1-1", melt_a12_4590$variable)
melt_a12_4590$variable =  gsub("X3C3.A", "3C.3A", melt_a12_4590$variable)
melt_a12_4590$variable =  gsub("T131K_R142K_R261Q", "3C.2A2-2", melt_a12_4590$variable)
melt_a12_4590$variable =  gsub("T131K_R142K", "3C.2A2-1", melt_a12_4590$variable)
melt_a12_4590$variable =  gsub("N121K_S144K", "3C.2A3", melt_a12_4590$variable)
melt_a12_4590$variable = factor(melt_a12_4590$variable, 
                                levels = c("3C.3A", "3C.2A", "3C.2A1-1", "3C.2A1-2", "3C.2A1-3",
                                           "3C.2A2-1", "3C.2A2-2", "3C.2A3"))


ggplot(melt_a12_4590, aes(x=N121K_N171K, y=value)) +
  geom_jitter(width=0.1, height=0.1, alpha=0.25) +
  facet_wrap(~variable)+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0,10),
                     labels=c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0,10),
                     labels=c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)) +
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y=10, label.x=5, aes(label = ..eq.label..), size=2) +
  stat_regline_equation(label.y=9, label.x=5, aes(label = ..rr.label..), size=2) +
  
  xlab("3C.2A1-2 titer") + ylab("Titer") +
  cor_theme


ggsave("../figure/A1-2_titers_x_other_titers_4590.png", height=5.5, width=5)
ggsave("../figure/A1-2_titers_x_other_titers_4590.pdf", height=5.5, width=5)
ggsave("../figure/A1-2_titers_x_other_titers_4590.tiff", height=5.5, width=5)

##### A3

melt_a3 = melt(ag_sera, id=c("Age_group", "N121K_S144K", "ID"))
melt_a3_4590 = melt_a3[melt_a3$Age_group == "(44,64]" | melt_a3$Age_group == "(64,90]", ]

melt_a3_4590$variable =  gsub("X3C3.A", "3C.3A", melt_a3_4590$variable)
melt_a3_4590$variable =  gsub("X3C2.A", "3C.2A", melt_a3_4590$variable)
melt_a3_4590$variable =  gsub("N121K_T135K_N171K", "3C.2A1-3", melt_a3_4590$variable)
melt_a3_4590$variable =  gsub("N121K_N171K", "3C.2A1-2", melt_a3_4590$variable)
melt_a3_4590$variable =  gsub("N171K", "3C.2A1-1", melt_a3_4590$variable)
melt_a3_4590$variable =  gsub("T131K_R142K_R261Q", "3C.2A2-2", melt_a3_4590$variable)
melt_a3_4590$variable =  gsub("T131K_R142K", "3C.2A2-1", melt_a3_4590$variable)
melt_a3_4590$variable =  gsub("N121K_S144K", "3C.2A3", melt_a3_4590$variable)
melt_a3_4590$variable = factor(melt_a3_4590$variable, 
                               levels = c("3C.3A", "3C.2A", "3C.2A1-1", "3C.2A1-2", "3C.2A1-3",
                                          "3C.2A2-1", "3C.2A2-2", "3C.2A3"))


ggplot(melt_a3_4590, aes(x=N121K_S144K, y=value)) +
  geom_jitter(width=0.1, height=0.1, alpha=0.25) +
  facet_wrap(~variable)+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0,10),
                     labels=c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0,10),
                     labels=c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)) +
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y=10, label.x=5, aes(label = ..eq.label..), size=2) +
  stat_regline_equation(label.y=9, label.x=5, aes(label = ..rr.label..), size=2) +
  
  xlab("3C.2A3 titer") + ylab("Titer") +
  cor_theme


ggsave("../figure/A3_titers_x_other_titers_4590.png", height=5.5, width=5)
ggsave("../figure/A3_titers_x_other_titers_4590.pdf", height=5.5, width=5)
ggsave("../figure/A3_titers_x_other_titers_4590.tiff", height=5.5, width=5)


#In 18-90, A1-3 and A2-2 have low correlation to other titers
# A1-3 shows high titers

melt_a13 = melt(ag_sera, id=c("Age_group", "N121K_T135K_N171K", "ID"))
melt_a13_1890 = melt_a13[melt_a13$Age_group == "(17,44]" |
                          melt_a13$Age_group == "(44,64]" | 
                           melt_a13$Age_group == "(64,90]", ]
melt_a13_1890$variable =  gsub("X3C2.A", "3C.2A", melt_a13_1890$variable)
melt_a13_1890$variable =  gsub("N121K_T135K_N171K", "3C.2A1-3", melt_a13_1890$variable)
melt_a13_1890$variable =  gsub("N121K_N171K", "3C.2A1-2", melt_a13_1890$variable)
melt_a13_1890$variable =  gsub("N171K", "3C.2A1-1", melt_a13_1890$variable)
melt_a13_1890$variable =  gsub("X3C3.A", "3C.3A", melt_a13_1890$variable)
melt_a13_1890$variable =  gsub("T131K_R142K_R261Q", "3C.2A2-2", melt_a13_1890$variable)
melt_a13_1890$variable =  gsub("T131K_R142K", "3C.2A2-1", melt_a13_1890$variable)
melt_a13_1890$variable =  gsub("N121K_S144K", "3C.2A3", melt_a13_1890$variable)
melt_a13_1890$variable = factor(melt_a13_1890$variable, 
                                levels = c("3C.3A", "3C.2A", "3C.2A1-1", "3C.2A1-2", "3C.2A1-3",
                                           "3C.2A2-1", "3C.2A2-2", "3C.2A3"))


ggplot(melt_a13_1890, aes(x=N121K_T135K_N171K, y=value)) +
  geom_jitter(width=0.1, height=0.1, alpha=0.25) +
  facet_wrap(~variable)+
 
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0,10),
                     labels=c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0,10),
                     labels=c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)) +
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y=10, label.x=5, aes(label = ..eq.label..), size=2) +
  stat_regline_equation(label.y=9, label.x=5, aes(label = ..rr.label..), size=2) +
  xlab("3C.2A1-3 titer") + ylab("Titer") +
  cor_theme


ggsave("../figure/A1-3_titers_x_other_titers_1890.png", height=5.5, width=5)
ggsave("../figure/A1-3_titers_x_other_titers_1890.pdf", height=5.5, width=5)
ggsave("../figure/A1-3_titers_x_other_titers_1890.tiff", height=5.5, width=5)


###### A2-2
melt_a22 = melt(ag_sera, id=c("Age_group", "T131K_R142K_R261Q", "ID"))
melt_a22_4590 = melt_a22[melt_a22$Age_group == "(17,44]"|
  melt_a22$Age_group == "(44,64]" | melt_a22$Age_group == "(64,90]", ]

melt_a22_4590$variable =  gsub("X3C2.A", "3C.2A", melt_a22_4590$variable)
melt_a22_4590$variable =  gsub("N121K_T135K_N171K", "3C.2A1-3", melt_a22_4590$variable)
melt_a22_4590$variable =  gsub("N121K_N171K", "3C.2A1-2", melt_a22_4590$variable)
melt_a22_4590$variable =  gsub("N171K", "3C.2A1-1", melt_a22_4590$variable)
melt_a22_4590$variable =  gsub("X3C3.A", "3C.3A", melt_a22_4590$variable)
melt_a22_4590$variable =  gsub("T131K_R142K_R261Q", "3C.2A2-2", melt_a22_4590$variable)
melt_a22_4590$variable =  gsub("T131K_R142K", "3C.2A2-1", melt_a22_4590$variable)
melt_a22_4590$variable =  gsub("N121K_S144K", "3C.2A3", melt_a22_4590$variable)
melt_a22_4590$variable = factor(melt_a22_4590$variable, 
                                levels = c("3C.3A", "3C.2A", "3C.2A1-1", "3C.2A1-2", "3C.2A1-3",
                                           "3C.2A2-1", "3C.2A2-2", "3C.2A3"))


ggplot(melt_a22_4590, aes(x=T131K_R142K_R261Q, y=value)) +
  geom_jitter(width=0.1, height=0.1, alpha=0.25) +
  facet_wrap(~variable)+

  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0,10),
                     labels=c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0,10),
                     labels=c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)) +
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(label.y=10, label.x=5, aes(label = ..eq.label..), size=2) +
  stat_regline_equation(label.y=9, label.x=5, aes(label = ..rr.label..), size=2) +
  labs(x="3C.2A2-2 titer", y="Titer") +
  cor_theme


ggsave("../figure/A2-2_titers_x_other_titers_1890.png", height=5.5, width=5)
ggsave("../figure/A2-2_titers_x_other_titers_1890.pdf", height=5.5, width=5)
ggsave("../figure/A2-2_titers_x_other_titers_1890.tiff", height=5.5, width=5)
