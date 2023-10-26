library(ggplot2)
library(ggtree)
library(tidytree)
library(ape)
library(treeio)

fp = file.path("na_sample_for_genealogy_1217_mcc.nex")
info = read.csv("na_sample_for_genealogy_1217_clade.csv", stringsAsFactors = F)

tree <- read.beast(fp) 
tt = as_tibble(tree)

y = full_join(tt, as_tibble(info), by="label")

#find mrca of each clade
clades = c("c3a", "c2a", "A1", "A1_2", "A1_3", "A2", "A2_2", "A3")

clades_mrca = c()
for (c in 1:length(clades)) {
  print (c)
  c_nodes = c(y$clade == clades[c])
  c_nodes = ifelse(c_nodes == F, NA, c_nodes)
  
  if (sum(c_nodes, na.rm=T) == 0){
    c_mrca = NA
    clades_mrca = c(clades_mrca, c_mrca)
    next
  }
  
  c_mrca = MRCA(y, which(y$clade == clades[c]))
  c_mrca = c_mrca$node
  clades_mrca = c(clades_mrca, c_mrca)  
}

clades_mrca = clades_mrca[!(is.na(clades_mrca))]

mrca_3c3a = clades_mrca[1]
mrca_3c2a = clades_mrca[2]
mrca_a1 = MRCA(y, clades_mrca[3], clades_mrca[4], clades_mrca[5])$node
mrca_a2 = MRCA(y, clades_mrca[6], clades_mrca[7])$node
mrca_a3 = clades_mrca[8]
clades_mrca = c(mrca_3c3a, mrca_3c2a, mrca_a1, mrca_a2, mrca_a3)
clades_mrca_name = c("3C3.A", "3C2.A", "3C2.A1", "3C2.A2", "3C2.A3")

#######################

#group clades based on mrca nodes

y2 = groupClade(y, clades_mrca) #this adds "group" variable to y and make y2
y2$clade = factor(y2$clade, levels = c("c3a", "c2a", 
                                          "A1", "A1_2", "A1_3", 
                                          "A2", "A2_2", "A3"))

#for 1718 season North America nodes, identify clade
ss = c(y2$season == "1718") * c(y2$location == "North America ")
y2$na_clade_1718 = ifelse(ss == 1, as.character(y2$clade), NA)
y2$na_clade_1718 = factor(y2$na_clade_1718, levels = c("c3a", "c2a", 
                                       "A1", "A1_2", "A1_3", 
                                       "A2", "A2_2", "A3"))

ss1617 = c(y2$season == "1617") * c(y2$location == "North America ")
y2$na_clade_1617 = ifelse(ss1617 == 1, as.character(y2$clade), NA)
y2$na_clade_1617 = factor(y2$na_clade_1617, levels = c("c3a", "c2a", 
                                                       "A1", "A1_2", "A1_3", 
                                                       "A2", "A2_2", "A3"))

y2$show1617 = ifelse(y2$height > 1.0, T, F )

#for text of clade name on tree figure
y2$clade_label = y2$node %in% clades_mrca
y2$clade_label = ifelse(y2$clade_label == F, NA, y2$node)
for (i in 1:length(clades_mrca_name)) {
  y2$clade_label[clades_mrca[i]] = clades_mrca_name[i]
  
}

################################################################

atree = as.treedata(y2)

btree=atree
##############################################################

#MRCA(btree, 28, 29)
#child(btree, 300)
#y2$season

#child(y2, 3)
p1617 = ggtree(btree, aes(color=show1617), mrsd = '2018-09-29') +
  geom_tippoint(aes(fill=na_clade_1617), shape=21, size=3, color="transparent") +
  scale_fill_manual(values = c("#F8766D", "#CD9600",
                               "#7CAE00", "#00BE67", "#00BFC4", 
                               "#6666FF", "#C77CFF", "#FF61CC", "white"),
                    labels = c("3C.3A", "3C.2A", "3C.2A1-1", "3C.2A1-2", "3C.2A1-3",
                               "3C.2A2-1", "3C.2A2-2", "3C.2A3", ""),
                    name = "",
                    na.value=NA) +
  scale_color_manual(values = c( "#dddddd", "#444444"), guide = "none") +  
  theme_tree2(panel.grid.major.x=element_line(size=0.5, linetype="solid", color="light gray")) +
  scale_x_continuous(limits = c(2010.6, 2019), breaks=seq(2011, 2018), labels = seq(2011, 2018)) 
p1617
saveRDS(p1617, "na_tree_1617.rds")
ggsave("na_tree_1617.png", height=6, width=5)
ggsave("na_tree_1617.pdf", height=6, width=5)


p1617 + scale_fill_manual(values = c("#F8766D", "#CD9600", 
                                     "#7CAE00", "#00BE67", "#00BFC4", 
                                     "#6666FF", "#C77CFF", "#FF61CC"),
                          labels = c("3C.3A", "3C.2A", "3C.2A1-1", "3C.2A1-2", "3C.2A1-3",
                                     "3C.2A2-1", "3C.2A2-2", "3C.2A3", ""),
                          name = "",
                          na.value=NA,
                          guide="none") 
ggsave("na_tree_1617_nolegend.png", height=6, width=5)
ggsave("na_tree_1617_nolegend.pdf", height=6, width=5)




p1718 = ggtree(btree, color="#444444", mrsd = '2018-09-29') +

  geom_tippoint(aes(fill=na_clade_1718), shape=21, size=3, color="transparent") +
  scale_fill_manual(values = c("#F8766D",  
                               "#00BE67", "#00BFC4", 
                               "#C77CFF", "#FF61CC"),
                    labels = c("3C.3A",  "3C.2A1-2", "3C.2A1-3",
                                "3C.2A2-2", "3C.2A3"),
                    name = "",
                    na.value=NA) +
  theme_tree2(panel.grid.major.x=element_line(size=0.5, linetype="solid", color="light gray")) +
  scale_x_continuous(limits = c(2010.6, 2019), breaks=seq(2011, 2019), labels = seq(2011, 2019)) +
  coord_cartesian(clip = "off")

p1718

saveRDS(p1718, "na_tree_1718.rds")
ggsave("na_tree_1718.png", height=6, width=5)
ggsave("na_tree_1718.pdf", height=6, width=5)


p1718 + scale_fill_manual(values = c("#F8766D",  
                                     "#00BE67", "#00BFC4", 
                                     "#C77CFF", "#FF61CC"),
                          labels = c("3C3.A",  "A1-2", "A1-3",
                                     "A2-2", "A3"),
                          name = "",
                          na.value=NA,
                          guide="none")
ggsave("na_tree_1718_nolegend.png", height=6, width=5)
ggsave("na_tree_1718_nolegend.pdf", height=6, width=5)


#library(scales)
#show_col(hue_pal()(8))
