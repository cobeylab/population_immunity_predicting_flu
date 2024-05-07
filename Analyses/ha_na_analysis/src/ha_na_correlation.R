library(ggplot2)
library(ggpubr)


ag_brk = c(0,4,17,44,64,90)
sq = seq(0, 10, 1)
lb = 2^(sq) * 10


#####################################################################################

source("../../../Data/script/load_data/load_data.R")
source("../../../Data/script/load_data/load_na.R")

df_ha = df[df$Test_virus == "3C2.A" | df$Test_virus == "T131K_R142K_R261Q", ]
df_na = df_na[df_na$Test_virus == "3c2.A" | df_na$Test_virus == "A2", ]

df_ha$winner = ifelse(df_ha$Test_virus == "3C2.A", "F", "T")
df_na$winner = ifelse(df_na$Test_virus == "3c2.A", "F", "T")

df = merge(df_ha, df_na, by= c("Sample_ID", "Age_group", "Age"))

df$label = df$Age_group
df$label = gsub("\\(0,4\\]", "1-4 years", df$label)
df$label = gsub("\\(4,17\\]", "5-17 years", df$label)
df$label = gsub( "\\(17,44\\]", "18-44 years", df$label)
df$label = gsub( "\\(44,64\\]", "45-64 years", df$label)
df$label = gsub("\\(64,90\\]", "65-90 years", df$label)

df$label = factor(df$label, levels = c("1-4 years", "5-17 years", "18-44 years", 
                                                               "45-64 years", "65-90 years"))


#######################################################################################


cor_theme = theme_bw() + theme(
  axis.text.x = element_text(angle=300, vjust=0.7, hjust=0, size=7.5),
  axis.text.y = element_text(size=8),
  axis.title = element_text(size=8),
  strip.text = element_text(size=8),
  legend.title=element_text(size=8),
  legend.text = element_text(size=8),
  legend.key.size = unit(0.35, "cm"),
  legend.margin = unit(-0.2, "cm"),
  #panel.grid = element_blank(),
  strip.background = element_rect(fill = "white", colour = "white"))

#######################################################################################
# remove individuals who have undetectable titers to both HA and NA

df_rmv = df[ (df$Titer.x != 0) | (df$Titer.y != 0) , ]
cor.test(df_rmv$Titer.x, df_rmv$Titer.y)


p_rmv_wt = ggplot(df_rmv[df_rmv$winner.x == "F", ], aes(x=Titer.x, y=Titer.y)) +
  geom_jitter(alpha=0.2, width=0.15, height=0.15) +
  facet_wrap(~label) +
  labs(x="HA Titer", y="NA Titer") +
  scale_y_continuous(limits = c(0, 10), breaks = sq, labels = lb) +
  scale_x_continuous(limits = c(0, 10), breaks = sq, labels = lb) +
  #theme_bw() +
  cor_theme +
  theme(legend.position="none") +
  stat_cor(aes(label = ..r.label..), method = 'spearman', cor.coef.name = 'rho')


p_rmv_wn = ggplot(df_rmv[df_rmv$winner.x == "T", ], aes(x=Titer.x, y=Titer.y)) +
  geom_jitter(alpha=0.2, width=0.15, height=0.15) +
  facet_wrap(~label) +
  labs(x="HA Titer", y="NA Titer") +
  scale_y_continuous(limits = c(0, 10), breaks = sq, labels = lb) +
  scale_x_continuous(limits = c(0, 10), breaks = sq, labels = lb) +
  #theme_bw() +
  cor_theme +
  theme(legend.position="none") +
  stat_cor(aes(label = ..r.label..), method = 'spearman', cor.coef.name = 'rho')


ggarrange( p_rmv_wt, p_rmv_wn,
           nrow=2,
           labels = c("A", "B")
)

#ggsave(paste0("../figure/ha_na_correlation_rmv_all0s.png"), height=8, width=5.2)


# correlation coefficient and p value 
coeffs = c()
pvals = c()
ag = levels(df_rmv$Age_group)
for (i in 1:length(ag)) {

    df_ag = df_rmv[df_rmv$Age_group == ag[i] & df_rmv$winner.x == "F",]
    ct = cor.test(df_ag$Titer.x, df_ag$Titer.y)
    coeff = ct[[4]]
    coeff = round(coeff, 2)
    pval = ct[[3]]
    pval = ifelse(pval<0.001, "<0.001", round(pval,4) )
      
    coeffs = rbind(coeffs, coeff)
    pvals = rbind(pvals, pval)
}
rmv_wt_r2 = cbind(coeffs, pvals)
colnames(rmv_wt_r2) = c("coefficient", "p value")
rownames(rmv_wt_r2) = c("1-4 years", "5-17 years", 
                        "18-44 years", "45-64 years", "65-90 years")

coeffs = c()
pvals = c()
ag = levels(df_rmv$Age_group)
for (i in 1:length(ag)) {
  
  df_ag = df_rmv[df_rmv$Age_group == ag[i] & df_rmv$winner.x == "T",]
  ct = cor.test(df_ag$Titer.x, df_ag$Titer.y)
  coeff = ct[[4]]
  coeff = round(coeff, 2)
  pval = ct[[3]]
  pval = ifelse(pval<0.001, "<0.001", round(pval,4) )
  
  coeffs = rbind(coeffs, coeff)
  pvals = rbind(pvals, pval)
}
rmv_wn_r2 = cbind(coeffs, pvals)
colnames(rmv_wn_r2) = c("coefficient", "p value")
rownames(rmv_wn_r2) = c("1-4 years", "5-17 years", 
                        "18-44 years", "45-64 years", "65-90 years")

rmv_r2 = rbind(rmv_wt_r2, rmv_wn_r2)
#write.csv(rmv_r2, "../stat_test/rmv_r2.csv")

##########################################################################



cor.test(df$Titer.x, df$Titer.y)

p_wt = ggplot(df[df$winner.x == "F", ], aes(x=Titer.x, y=Titer.y)) +
  geom_jitter(alpha = 0.2, width=0.15, height=0.15) +
  facet_wrap(~label) +
  labs(x="HA Titer", y="NA Titer") +
  scale_y_continuous(limits = c(0, 10), breaks = sq, labels = lb) +
  scale_x_continuous(limits = c(0, 10), breaks = sq, labels = lb) +
  #theme_bw() +
  cor_theme +
  theme(legend.position="none") +
  stat_cor(aes(label = ..r.label..), method = 'spearman', cor.coef.name = 'rho')

p_wn = ggplot(df[df$winner.x == "T", ], aes(x=Titer.x, y=Titer.y)) +
  geom_jitter(alpha = 0.2, width=0.15, height=0.15) +
  facet_wrap(~label) +
  labs(x="HA Titer", y="NA Titer") +
  scale_y_continuous(limits = c(0, 10), breaks = sq, labels = lb) +
  scale_x_continuous(limits = c(0, 10), breaks = sq, labels = lb) +
  #theme_bw() +
  cor_theme +
  theme(legend.position="none") +
  stat_cor(aes(label = ..r.label..), method = 'spearman', cor.coef.name = 'rho')



ggarrange( p_wt, p_wn,
           nrow=2,
           labels = c("A", "B")
)

ggsave(paste0("../figure/ha_na_correlation.png"), height=8, width=5.2)


# correlation coefficient and p value
coeffs = c()
pvals = c()
ag = levels(df$Age_group)
for (i in 1:length(ag)) {
  
  df_ag = df[df$Age_group == ag[i] & df$winner.x == "F",]
  ct = cor.test(df_ag$Titer.x, df_ag$Titer.y)
  coeff = ct[[4]]
  coeff = round(coeff, 2)
  pval = ct[[3]]
  pval = ifelse(pval<0.001, "<0.001", round(pval,4) )
  
  coeffs = rbind(coeffs, coeff)
  pvals = rbind(pvals, pval)
}
wt_r2 = cbind(coeffs, pvals)
colnames(wt_r2) = c("coefficient", "p value")
rownames(wt_r2) = c("1-4 years", "5-17 years", 
                        "18-44 years", "45-64 years", "65-90 years")


coeffs = c()
pvals = c()
ag = levels(df$Age_group)
for (i in 1:length(ag)) {
  
  df_ag = df[df$Age_group == ag[i] & df$winner.x == "T",]
  ct = cor.test(df_ag$Titer.x, df_ag$Titer.y)
  coeff = ct[[4]]
  coeff = round(coeff, 2)
  pval = ct[[3]]
  pval = ifelse(pval<0.001, "<0.001", round(pval,4) )
  
  coeffs = rbind(coeffs, coeff)
  pvals = rbind(pvals, pval)
}
wn_r2 = cbind(coeffs, pvals)
colnames(wn_r2) = c("coefficient", "p value")
rownames(wn_r2) = c("1-4 years", "5-17 years", 
                    "18-44 years", "45-64 years", "65-90 years")

r2 = rbind(wt_r2, wn_r2)
#write.csv(r2, "../stat_test/r2.csv")
