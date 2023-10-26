source("analysis_util.R")
source("clustering_functions.R")


ag_brk = c(0, 4, 17, 44, 64, 90)

##################################################################
# Data: 
#tsera (HA titers, each column each test virus)
# ida_titer (sample id, age, ha titers)  
####################################################################

source("../../../Data/script/load_data/load_data.R")
  
# make data for analysis 

ida_titer = sera[,c("Sample_ID", "Age", "3C3.A", "3C2.A", "N171K", "N121K_S144K", "N121K_N171K", "T131K_R142K", 
                "N121K_T135K_N171K", "T131K_R142K_R261Q")]
ida_titer = na.omit(ida_titer)

tsera = (ida_titer[,3:ncol(ida_titer)])
  
###########################################################################
# Clustering using normalized titers (by individual)
#############################################################################

analysis = "real_try1"

# scale

tz = data.frame(matrix(nrow=0, ncol=8))
colnames(tz) = colnames(tsera)
included_ids = c()
included_ages = c()

for(i in 1:nrow(tsera)){
  row1 = tsera[i,]
  if (sum(is.na(row1)) != 0){ #if there is na in the row, skip
    next
  }
  if (sum(row1) == 0) {# if all titers are undetectable, skip
    next
  }
  
  #normalize (l2) to use cosine similarity for kmeans clustering
  row1 = row1/sqrt(sum(row1^2))
  tz = rbind(tz, row1)
  included_ids = c(included_ids, ida_titer[i,"Sample_ID"])
  included_ages = c(included_ages, ida_titer[i, "Age"])
}

rownames(tz) = as.character(included_ids)


# k-means clustering
cl_fit = cascadeKM(tz, 1, 7, ite=1500)

#assign cluster
ida_titer = assign_cluster(ida_titer, cl_fit, fig_dir, analysis)
cluster_original = ida_titer$cluster

# result support 2 groups
unique(ida_titer$cluster)

####################################################################3
# The group with high 3C.3A is shown as group 1 in the main text and this group has smaller number of people
# Set the group with smaller number of people as "group 1"

num_group1 = sum(ida_titer$cluster == "Group 1")
num_group2 = sum(ida_titer$cluster == "Group 2")
if (num_group1 > num_group2){
  ida_titer$cluster = ifelse(cluster_original == "Group 1", "Group 2", "Group 1")
}
ida_titer$cluster = factor(ida_titer$cluster)


# Test if there is age difference by cluster
sm = summary( lm(Age ~ cluster, data=ida_titer) )
sm$coefficients

# save clustering result and statistical test result

saveRDS(ida_titer, "../result/ida_titer_by_cluster.rds")

write.csv(sm$coefficients, "../stat_test/age_diff_real_cosine_similarity.csv")
  

############################################################################
# plot clusters

br = seq(0, 11)
lb = (2^(br))*10
fig_dir = "../figure/"


plot_x_label = c("3C.3A", "3C.2A", "3C.2A1-1", "3C.2A1-2", "3C.2A1-3",
                 "3C.2A2-1", "3C.2A2-2", "3C.2A3")

m_ida_titer = melt(ida_titer, id.vars= c("Sample_ID", "Age", "cluster"))
colnames(m_ida_titer)[c(4,5)] = c("Virus", "Titer")
m_ida_titer$Virus = factor(m_ida_titer$Virus, 
                           levels = c("3C3.A", "3C2.A", "N171K", "N121K_N171K", "N121K_T135K_N171K",
                                      "T131K_R142K", "T131K_R142K_R261Q", "N121K_S144K"))


theme_cl = theme_bw() +  
      theme(axis.text = element_text(size=8),
      axis.title = element_text(size=9),
      
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size=7),
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-7,-5,-7,-5),
      
      panel.grid = element_blank(),
      
      strip.text = element_text(size=8),
      strip.background = element_rect(fill = "white", colour = "white"))

############################################################################

p_titer = ggplot(m_ida_titer, aes(x=Age, y=Titer)) +
  geom_jitter(aes(col=Virus), width=0.2, height=0.2, alpha=0.2) +
  geom_smooth(aes(col=Virus), alpha=0, size=1.5) +
  facet_wrap(~cluster) +
  scale_color_manual(labels = plot_x_label, name="",
                     values = c("#F8766D", "#CD9600",
                                "#7CAE00", "#00BE67", "#00BFC4", 
                                "#6666FF", "#C77CFF", "#FF61CC")) +
  scale_y_continuous(limits = c(0, 9), breaks = br, labels = lb) +
  xlab("Age (years)") +
  theme_cl 


p_dist = ggplot(ida_titer) +
  geom_histogram(aes(x=Age)) +
  facet_wrap(~cluster) +
  ylab("Count") +
  xlab("Age (years)") +
  theme_cl +
  theme(plot.margin = unit(c(0.2, 3.0, 0.2, 0.5), "lines"))



library(ggpubr)

ggarrange(p_titer, p_dist,
          nrow=2, ncol=1,
          labels=c("A", "B"))


ggsave(paste0(fig_dir, "Clustering_scale_by_indiv_real_cosine_similarity.png"),
       width=4, height=4)
ggsave(paste0(fig_dir, "Clustering_scale_by_indiv_real_cosine_similarity.pdf"),
       width=4, height=4)
ggsave(paste0(fig_dir, "Clustering_scale_by_indiv_real_cosine_similarity.tiff"),
       width=4, height=4)


