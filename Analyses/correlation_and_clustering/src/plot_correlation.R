library(ggplot2)

ag_cormat = readRDS("../result/ag_cormat_real_rmv_undetectable.rds")
  

ag_cormat$label = ag_cormat$Age_group
ag_cormat$label = gsub("\\(0,4\\]", "1-4 years", ag_cormat$label)
ag_cormat$label = gsub("\\(4,17\\]", "5-17 years", ag_cormat$label)
ag_cormat$label = gsub( "\\(17,44\\]", "18-44 years", ag_cormat$label)
ag_cormat$label = gsub( "\\(44,64\\]", "45-64 years", ag_cormat$label)
ag_cormat$label = gsub("\\(64,90\\]", "65-90 years", ag_cormat$label)

ag_cormat$label = factor(ag_cormat$label, levels = c("1-4 years", "5-17 years", 
                                                     "18-44 years", "45-64 years", "65-90 years"))


##############################################################



plot_label = c("3C.3a", "3C.2a", "3C.2a3","3C.2a1-1", "3C.2a1-2", "3C.2a1-3",
               "3C.2a2-1", "3C.2a2-2")

theme_cor =   theme_bw() + theme(axis.text.x = element_text(angle=290, vjust=0.75, hjust=0),
                    text = element_text(size=8),
                    panel.spacing = unit(0.5, "lines"), 
                    legend.margin=margin(0,0,0,0), 
                    legend.box.margin=margin(5,5,5,12),
                    legend.position = c(0.9, 0.13),
                    legend.key.size = unit(0.4, "cm"),
                    
                    panel.grid = element_blank(),
                    strip.background = element_rect(fill = "white", colour = "white"))

#1. Correlation heatmap



ggplot(ag_cormat, aes(Virus1, Virus2, fill=r)) +
  geom_tile(col="white") +
  scale_fill_gradient2(low="white", mid = "light blue", high="blue",  
                       midpoint = 0.75, 
                       name="Spearman rho") +
  coord_fixed() + 
  facet_wrap(~label) +
  ylab("") + xlab("") +
  scale_x_discrete(labels = plot_label) +
  scale_y_discrete(labels = rev( plot_label) )  +
  theme_cor

ggsave(paste0("../figure/Correlation_heatmap_spearman_real_rmv_undetectable.png"),
              height=3.5, width = 5)
ggsave(paste0("../figure/Correlation_heatmap_spearman_real_rmv_undetectable.pdf"),
       height=3.5, width = 5)
ggsave(paste0("../figure/Correlation_heatmap_spearman_real_rmv_undetectable.tiff"),
       height=3.5, width = 5)


# 2.1. The age group greater than 1-4 


# ggplot(ag_cormat, aes(Virus1, Virus2, fill=r_ccg)) +
#   geom_tile(col="white") +
#   scale_fill_gradient2(low="white", mid = "light blue", high="blue",  
#                        midpoint = 0.75, 
#                        name="Spearman rho") +
#   coord_fixed() + 
#   facet_wrap(~label) +
#   ylab("") + xlab("") +
#   scale_x_discrete(labels = plot_label) +
#   scale_y_discrete(labels = rev(plot_label) ) +
#   theme_cor

# ggsave(paste0("../figure/Correlation_heatmap_spearman_greater_than_baseline_real_rmv_undetectable.png"),
#        height=3.5, width = 5)
# ggsave(paste0("../figure/Correlation_heatmap_spearman_greater_than_baseline_real_rmv_undetectable.pdf"),
#        height=3.5, width = 5)
# ggsave(paste0("../figure/Correlation_heatmap_spearman_greater_than_baseline_real_rmv_undetectable.tiff"),
#        height=3.5, width = 5)


# 2.2. The age group weaker than 1-4


ggplot(ag_cormat, aes(Virus1, Virus2, fill=r_ccl)) +
  geom_tile(col="white") +
  scale_fill_gradient2(low="white", mid = "light blue", high="blue",  
                       midpoint = 0.75, 
                       name="Spearman rho") +
  coord_fixed() + 
  facet_wrap(~label) +
  ylab("") + xlab("") +
  scale_x_discrete(labels = plot_label) +
  scale_y_discrete(labels = rev(plot_label) ) +
  theme_cor

ggsave(paste0("../figure/Correlation_heatmap_spearman_weaker_than_baseline_real_rmv_undetectable.png"),
       height=3.5, width = 5)
ggsave(paste0("../figure/Correlation_heatmap_spearman_weaker_than_baseline_real_rmv_undetectable.pdf"),
       height=3.5, width = 5)
ggsave(paste0("../figure/Correlation_heatmap_spearman_weaker_than_baseline_real_rmv_undetectable.tiff"),
       height=3.5, width = 5)

# Within an age group, rank of correlation

ggplot(ag_cormat, aes(Virus1, Virus2, fill=rank)) +
  geom_tile(col="white") +
  scale_fill_gradient2(low="white",  high="red",  
                       name=" Number of other pairs \n with stronger correlations") +
  coord_fixed() + 
  facet_wrap(~label) +
  ylab("") + xlab("") +
  scale_x_discrete(labels = plot_label) +
  scale_y_discrete(labels = rev(plot_label) ) +
  theme_cor

ggsave(paste0("../figure/Correlation_heatmap_spearman_rank_real_rmv_undetectable.png"),
       height=3.5, width = 5)
ggsave(paste0("../figure/Correlation_heatmap_spearman_rank_real_rmv_undetectable.pdf"),
       height=3.5, width = 5)
ggsave(paste0("../figure/Correlation_heatmap_spearman_rank_real_rmv_undetectable.tiff"),
       height=3.5, width = 5)




##################################################################

# #ggplot(ag_sera) +
#   geom_jitter(aes(x=X3C2.A, y=N121K_N171K), width=0.1, height=0.1, alpha=0.3) +
#   geom_jitter(aes(x=X3C2.A, y=X3C3.A), width=0.1, height=0.11, alpha=0.3, col="green") +
#   facet_wrap(~Age_group) +
#   ggtitle("3c3A")
# #ggsave("cor_3C3A.png")
# 
# ggplot(ag_sera) +
#   geom_jitter(aes(x=X3C2.A, y=N121K_N171K), width=0.1, height=0.1, alpha=0.3) +
#   geom_jitter(aes(x=X3C2.A, y=N121K_T135K_N171K), width=0.1, height=0.1, alpha=0.3, col="green") +
#   facet_wrap(~Age_group) +
#   ggtitle("121 135 171") 
# #ggsave("cor_121_135_171.png")
# 
# ggplot(ag_sera) +
#   geom_jitter(aes(x=X3C2.A, y=N121K_N171K), width=0.1, height=0.1, alpha=0.3) +
#   geom_jitter(aes(x=X3C2.A, y=T131K_R142K), width=0.1, height=0.1, alpha=0.3, col="green") +
#   facet_wrap(~Age_group) +
#   ggtitle("131 142")
# #ggsave("cor_131_142.png")
# 
# ggplot(ag_sera) +
#   geom_jitter(aes(x=X3C2.A, y=N121K_N171K), width=0.1, height=0.1, alpha=0.3) +
#   geom_jitter(aes(x=X3C2.A, y=T131K_R142K_R261Q), width=0.1, height=0.1, alpha=0.3, col="green") +
#   facet_wrap(~Age_group) +
#   ggtitle("131 142 261")
# #ggsave("cor_131_142_261.png")
