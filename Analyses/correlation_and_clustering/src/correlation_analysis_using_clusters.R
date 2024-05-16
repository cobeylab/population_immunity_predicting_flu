library(ggplot2)
library(dplyr)
library(reshape2)
source("analysis_util.R")

get_cormat = function(sera, cols, method){
  #sera: matrix of titers, each virus at each column
  cormat = cor(sera[,cols], method = method, use="pairwise.complete.obs")
  cormat[upper.tri(cormat)] = NA
  m_cormat = melt(cormat, na.rm=T)
  m_cormat
}

ida_titer = readRDS("../result/ida_titer_by_cluster.rds")


vlevels = c("3C3.A", "3C2.A", "N121K_S144K",
            "N171K",  "N121K_N171K", "N121K_T135K_N171K",
            "T131K_R142K", "T131K_R142K_R261Q")
ag_sera = ida_titer[, c("Sample_ID", "Age", vlevels, "cluster")]


#calculate correlation of each virus pairs for each cluster group

ag_cormat = ag_sera %>%
  group_by(cluster) %>%
  do(get_cormat(., c(3:10), "spearman"))

ag_cormat$label = ag_cormat$cluster

ag_cormat = ag_cormat[, c("cluster", "Var2", "Var1", "value", "label")]
colnames(ag_cormat) = c("Group", "Virus1", "Virus2", "r", "label")

saveRDS(ag_cormat, "../result/ag_cormat_by_cluster.rds")


#####################################################################################

plot_label = c("3C.3A", "3C.2A", "3C.2A3","3C.2A1-1", "3C.2A1-2", "3C.2A1-3",
               "3C.2A2-1", "3C.2A2-2")

theme_cor =   theme_bw() + theme(axis.text.x = element_text(angle=290, vjust=0.75, hjust=0),
                                 text = element_text(size=9),
                                 panel.spacing = unit(0.5, "lines"), 
                                 legend.margin=margin(0,0,0,0), 
                                 legend.box.margin=margin(1,1,1,1),
                                 #legend.position = c(0.9, 0.13),
                                 legend.key.size = unit(0.4, "cm"),
                                 
                                 panel.grid = element_blank(),
                                 strip.background = element_rect(fill = "white", colour = "white"))

#1. Correlation heatmap


ggplot(ag_cormat, aes(Virus1, Virus2, fill=r)) +
  geom_tile(col="white") +
  scale_fill_gradient2(low="white", mid = "light blue", high="blue",  
                       midpoint = 0.55, 
                       name="Spearman rho") +
  coord_fixed() + 
  facet_wrap(~label) +
  ylab("") + xlab("") +
  scale_x_discrete(labels = plot_label) +
  scale_y_discrete(labels = rev(plot_label), limits = rev ) +
  theme_cor

ggsave(paste0("../figure/Correlation_heatmap_spearman_real_by_cluster.png"),
       height=2.5, width = 5)



write.csv(ag_cormat, "../result/ag_cormat_by_cluster.csv")



