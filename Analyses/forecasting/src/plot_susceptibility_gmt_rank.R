######################################################################
#The user needs to set options as below.
######################################################################

# main
# main = 0: for figures going to supplement
# main = 1: for figures incorporating plots of observed titers in addition to 
#           plot of relative susceptibility and going to main text 
# opt_na
# opt_na = 0: HA titers
# opt_na = 1: NA titers
#
# threshold
# threshold = 0: using GMT to calculate relative susceptibility
# threshold = 1: using cutoff of 1:20 to calculate relative susceptibility
# threshold = 2: using cutoff of 1:40 
# threshold = 3: using cutoff of 1:80 
# threshold = 4: using cutoff of 1:160 
#
###########################################################################

library(ggplot2)
fig_dir = "../fig/"

main = 0
opt_na=0
threshold = 0


if(opt_na == 0){
  opt = "df_susc"
  out = "real"
}else if (opt_na == 1){
  opt = "df_susc_na_"
  out = "real_na"
} 
  


#######################################################################################

susc_theme = theme_bw() + theme(
  axis.text.x = element_text(angle=290, vjust=0.7, hjust=0, size=7),
  axis.text.y = element_text(size=8),
  axis.title = element_text(size=8),
  strip.text = element_text(size=8),
  legend.title=element_text(size=8),
  legend.text = element_text(size=8),
  legend.key.size = unit(0.35, "cm"),
  legend.margin = unit(-0.2, "cm"),
  panel.grid = element_blank(),
  strip.background = element_rect(fill = "white", colour = "white"))

titer_theme = theme_bw() + theme(
  axis.text.x = element_text(size=8),
  axis.text.y = element_text(size=8),
  axis.title = element_text(size=10),
  strip.text = element_text(size=7),
  legend.title=element_text(size=8),
  legend.text = element_text(size=8),
  legend.key.size = unit(0.35, "cm"),
  panel.grid = element_blank(),
  strip.background = element_rect(fill = "white", colour = "white"))



if(opt_na==0){
  
  susc_x_label = c("3C.3A", "3C.2A", "3C.2A1-1", "3C.2A1-2", "3C.2A1-3",
                   "3C.2A2-1", "3C.2A2-2", "3C.2A3")
  

  susc_color_manual =  scale_color_manual(values = c("#FED976",  
                                                     "#FEB24C", "#FC4E2A", "#E31A1C", "#B10026",
                                                     "#660000", "#160000", "grey70"),
                                          breaks= c(1,2,3,4,5,6,7,8),
                                          name = "Rank") 
  
}else if (opt_na==1) {
  
  susc_x_label = c("3C.2A (NA)", "3C.2A2-2 (NA)")
  
  susc_color_manual = scale_color_manual(values = c("#FED976", "#FEB24C"),
                                         breaks= c(1,2),
                                         name = "Rank") 
  
}

#######################################################################################

susc_brk = c(0,-1,-2,-3,-4,-5)
susc_lb = lb = 2^(-susc_brk) * 10
if (opt_na == 0){
  susc_lmt = c(-3, 0)
}else if (opt_na == 1){
  susc_lmt = c(-5.5, 0)
}

rname = paste0("../result/", opt, threshold, "_0.rds")
df_susc = readRDS(rname)
df_susc$rank = as.integer( rank (-df_susc$num_sig, ties.method="min" ))
df_susc$rank = as.factor(df_susc$rank)



p_all = ggplot(df_susc) +
  geom_point(aes(x=Clade, y=susceptibility, col=rank), size=2) +
  geom_errorbar(aes(x=Clade, ymin=ci1, ymax=ci2, col=rank), size=1) +
  ylab("GMT (inverted scale)") +
  xlab("") +
  scale_x_discrete(labels = susc_x_label) +
  scale_y_continuous( breaks = susc_brk, labels = susc_lb, limit=susc_lmt) +
  #coord_cartesian(ylim = c(-2,-0.7)) +
  susc_theme +
  susc_color_manual


#######################################################################################

rname1 = paste0("../result/", opt, threshold, "_1.rds")
ag1 = readRDS(rname1)
ag1$Age_group = "1-4 years"
ag1$rank = as.integer( rank (-ag1$num_sig, ties.method="min" ))
ag1$rank = as.factor(ag1$rank)


rname2 = paste0("../result/", opt, threshold, "_2.rds")
ag2 = readRDS(rname2)
ag2$Age_group = "5-17 years"
ag2$rank = as.integer( rank (-ag2$num_sig, ties.method="min" ))
ag2$rank = as.factor(ag2$rank)


rname3 = paste0("../result/", opt, threshold, "_3.rds")
ag3 = readRDS(rname3)
ag3$Age_group = "18-44 years"
ag3$rank = as.integer( rank (-ag3$num_sig, ties.method="min" ))
ag3$rank = as.factor(ag3$rank)


rname4 = paste0("../result/", opt, threshold, "_4.rds")
ag4 = readRDS(rname4)
ag4$Age_group = "45-64 years"
ag4$rank = as.integer( rank (-ag4$num_sig, ties.method="min" ))
ag4$rank = as.factor(ag4$rank)


rname5 = paste0("../result/", opt, threshold, "_5.rds")
ag5 = readRDS(rname5)
ag5$Age_group = "65-90 years"
ag5$rank = as.integer( rank (-ag5$num_sig, ties.method="min" ))
ag5$rank = as.factor(ag5$rank)



df_susc_ag = rbind(ag1, ag2, ag3, ag4, ag5)


df_susc_ag$Age_group = factor(df_susc_ag$Age_group, levels = c("1-4 years", "5-17 years", "18-44 years", 
                                                               "45-64 years", "65-90 years"))

if(opt_na == 0){
  df_susc_ag$rank = factor(df_susc_ag$rank, levels = c(1,2,3,4,5,6,7,8))
  
}else if (opt_na == 1){
  df_susc_ag$rank = factor(df_susc_ag$rank, levels = c(1,2))
  
}

p_ag = ggplot(df_susc_ag) +
  geom_point(aes(x=Clade, y=susceptibility, col=rank), size=1.2) +
  geom_errorbar(aes(x=Clade, ymin=ci1, ymax=ci2, col=rank), size=0.6) +
  facet_wrap(~Age_group) +
  ylab("GMT (inverted scale)") +
  xlab("") +
  scale_y_continuous(breaks = susc_brk, labels = susc_lb, limits = susc_lmt) +
  scale_x_discrete(labels = susc_x_label) +
  susc_theme +
  susc_color_manual 

#####################################################################

ag_brk = c(0,4,17,44,64,90)
sq = seq(0, 9, 1)
lb = 2^(sq) * 10
if(opt_na == 0) {
  
  source("../../../Data/script/load_data/load_data.R")
  
  p_obs_titer = ggplot(df) +
    geom_jitter(aes(x=Age, y=Titer, col=Test_virus), alpha=0.1, height=0.1) +
    geom_smooth(aes(x=Age, y=Titer, col=Test_virus), size=1.3, alpha=0) +
    scale_color_manual(name="", labels = susc_x_label,
                       values = c("#F8766D", "#CD9600",
                                  "#7CAE00", "#00BE67", "#00BFC4",
                                  "#6666FF", "#C77CFF", "#FF61CC")) +
    scale_y_continuous(limits = c(0, 9), breaks = sq, labels = lb) +
    xlab("Age (years)") +
    titer_theme
  
}else if (opt_na == 1) {
  source("../../../Data/script/load_data/load_na.R")
  
  df = df_na
  df$Test_virus = ifelse(df$Test_virus == "A2", "3C.2A2-2 (NA)", "3C.2A (NA)")
  
  p_obs_titer = ggplot(df) +
    geom_jitter(aes(x=Age, y=Titer, col=Test_virus), alpha=0.1, height=0.1) +
    geom_smooth(aes(x=Age, y=Titer, col=Test_virus), size=1.3, alpha=0) +
    scale_color_manual(values = c("#CD9600", "#C77CFF"), name = "", labels = susc_x_label) +
    scale_y_continuous(limits = c(0, 9), breaks = sq, labels = lb) +
    xlab("Age (years)")+
    titer_theme
  
  
}



#####################################################################
library(ggpubr)



if(main == 1){
  p_obs_titer = p_obs_titer + theme(plot.margin = unit(c(0.2, 4, 1, 4), "lines"))
  
  p_all = p_all + theme(plot.margin = unit(c(0.1, 1.0, 0.1, 0.1), "lines"),
                        legend.position = "none")
  
  p_ag = p_ag + theme(plot.margin = unit(c(0.1 ,0.1, 0.1, 0.1), "lines"))
  
 ggarrange( p_obs_titer,
            ggarrange( p_all, p_ag,
                       ncol = 2, nrow = 1, widths=c(17,30), heights=c(10,8),
                       labels = c("B"), label.y = 1.15),
            nrow=2,
            labels = "A"
 )

 ggsave(paste0("../fig/rank_susceptibility_", out, "_threshold", threshold, ".png"), height=4.6, width=5.2)
 ggsave(paste0("../fig/rank_susceptibility_", out, "_threshold", threshold, ".tiff"), height=4.6, width=5.2)
 ggsave(paste0("../fig/rank_susceptibility_", out, "_threshold", threshold, ".pdf"), height=4.6, width=5.2)
 
  
}else{

  p_all = p_all + theme(plot.margin = unit(c(0.1, 1.0, 0.1, 0.1), "lines"),
                        legend.position = "none")
  
  p_ag = p_ag + theme(plot.margin = unit(c(0.1 ,0.1, 0.1, 0.1), "lines"))
  
  ggarrange(p_all, p_ag, 
            ncol=2, nrow=1, widths=c(17,30),
            labels = c("A", "B"))
  
  ggsave(paste0("../fig/rank_susceptibility_", out, "_threshold", threshold, ".png"), height=2.5, width=5.2)
  ggsave(paste0("../fig/rank_susceptibility_", out, "_threshold", threshold, ".tiff"), height=2.5, width=5.2)
  ggsave(paste0("../fig/rank_susceptibility_", out, "_threshold", threshold, ".pdf"), height=2.5, width=5.2)
  
  
   
}




