library(ggplot2)
library(dplyr)
library(reshape2)
source("analysis_util.R")

ag_brk = c(0, 15, 30, 45, 60, 75, 90)

##############################################################
# Make data frame "ag_sera" (age group + HA titers)
################################################################

# Read titer data
source("../../../Data/script/load_data/load_data.R")


sera0 = sera

# transform function add 'X' in front of column name if it starts with number
sera0 = transform(sera0, Age_group = cut(Age, breaks = ag_brk))


# levels for viruses used for correlation analysis
vlevels = c("X3C3.A", "X3C2.A", "N121K_S144K",
            "N171K",  "N121K_N171K", "N121K_T135K_N171K",
            "T131K_R142K", "T131K_R142K_R261Q")


# Make data frame for analysis
# ag_sera: age group + titers

ag_sera = sera0[, c("Age_group", vlevels)]


##################################################################
#remove all undetectable
rs = rowSums(ag_sera[ , c(2:ncol(ag_sera))], na.rm = T)
ag_sera = ag_sera[rs != 0, ]


###################################################################
source("correlation_functions.R")

#1. no scale


ag_cormat = get_ag_cormat_rp(ag_sera, "spearman")
ag_cormat$label = ag_cormat$Age_group 

ag_cormat$label = gsub("\\(0,15\\]", "1-15 years", ag_cormat$label)
ag_cormat$label = gsub("\\(15,30\\]", "16-30 years", ag_cormat$label)
ag_cormat$label = gsub("\\(30,45\\]", "31-45 years", ag_cormat$label)
ag_cormat$label = gsub("\\(45,60\\]", "46-60 years", ag_cormat$label)
ag_cormat$label = gsub("\\(60,75\\]", "61-75 years", ag_cormat$label)
ag_cormat$label = gsub("\\(75,90\\]", "76-90 years", ag_cormat$label)


ag_cormat$label = factor(ag_cormat$label, levels = c("1-15 years", "16-30 years", "31-45 years",
                                                     "46-60 years", "61-75 years", "76-90 years"))


#########################################################################################


plot_label = c("3C.3A", "3C.2A", "3C.2A3","3C.2A1-1", "3C.2A1-2", "3C.2A1-3",
               "3C.2A2-1", "3C.2A2-2")


theme_cor_equ =   theme_bw() + theme(axis.text.x = element_text(angle=290, vjust=0.75, hjust=0),
                                 text = element_text(size=8),
                                 panel.spacing = unit(0.5, "lines"), 
                                 legend.margin=margin(1,1,1,1), 
                                 legend.box.margin=margin(1,1,1,1),
                                 #legend.position = c(1, 0.1),
                                 legend.key.size = unit(0.4, "cm"),
                                 
                                 #plot.margin = c(1,10,1,1),
                                 panel.grid = element_blank(),
                                 strip.background = element_rect(fill = "white", colour = "white"))


#########################################################################################

ggplot(ag_cormat, aes(Virus1, Virus2, fill=r)) +
  geom_tile(col="white") +
  scale_fill_gradient2(low="white", mid = "light blue", high="blue",  
                       midpoint = 0.75, 
                       name="Spearman rho") +
  coord_fixed() + 
  facet_wrap(~label) +
  ylab("") + xlab("") +
  scale_x_discrete(labels = plot_label) +
  scale_y_discrete(labels = rev( plot_label) ) +
  theme_cor_equ

ggsave(paste0("../figure/Correlation_heatmap_spearman_real_equal_bin_15_rmv_undetectable.png"),
       height=3.5, width = 5.2)
ggsave(paste0("../figure/Correlation_heatmap_spearman_real_equal_bin_15_rmv_undetectable.pdf"),
       height=3.5, width = 5.2)
ggsave(paste0("../figure/Correlation_heatmap_spearman_real_equal_bin_15_rmv_undetectable.tiff"),
       height=3.5, width = 5.2)

