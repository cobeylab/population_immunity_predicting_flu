library(ggplot2)
library(ggtree)
library(tidytree)
library(ape)
library(treeio)
library(tidytree)


fp = file.path("sample_for_genealogy_1217_mcc.nex")
info = read.csv("sample_for_genealogy_1217_clade.csv", stringsAsFactors = F)

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
mrca_a1 = MRCA(y, clades_mrca[3], clades_mrca[4])$node #Force MRCA to be 401, because MRCA function gives wrong mrca
mrca_a2 = MRCA(y, clades_mrca[6], clades_mrca[7])$node
mrca_a3 = clades_mrca[8]
clades_mrca = c(mrca_3c3a, mrca_3c2a, mrca_a1, mrca_a2, mrca_a3)
clades_mrca_name = c("3C.3a", "3C.2a", "3C.2a1", "3C.2a2", "3C.2a3")


#######################

y2 = groupClade(y, clades_mrca)
y2$clade = factor(y2$clade, levels = c("c3a", "c2a", 
                                       "A1", "A1_2", "A1_3", 
                                       "A2", "A2_2", "A3"))

ss = c(y2$season == "1718") * c(y2$location == "North America ")
y2$ha_clade_1718 = ifelse(ss == 1, as.character(y2$clade), NA)
y2$ha_clade_1718 = factor(y2$ha_clade_1718, levels = c("c3a", "c2a", 
                                                       "A1", "A1_2", "A1_3", 
                                                       "A2", "A2_2", "A3"))

ss1617 = c(y2$season == "1617") * c(y2$location == "North America ")
y2$ha_clade_1617 = ifelse(ss1617 == 1, as.character(y2$clade), NA)
y2$ha_clade_1617 = factor(y2$ha_clade_1617, levels = c("c3a", "c2a", 
                                                       "A1", "A1_2", "A1_3", 
                                                       "A2", "A2_2", "A3"))



#ss = c(y2$season == "1718") * c(y2$location == "North America ")
#y2$na_clade = ifelse(ss == 1, as.character(y2$clade), NA)
#y2$na_clade = factor(y2$na_clade, levels = c("c3a", "c2a", 
#                                             "A1", "A1_2", "A1_3", 
#                                             "A2", "A2_2", "A3"))


y2$clade_label = y2$node %in% clades_mrca
y2$clade_label = ifelse(y2$clade_label == F, NA, y2$node)
for (i in 1:length(clades_mrca_name)) {
  y2$clade_label[clades_mrca[i]] = clades_mrca_name[i]
  
}

levels(y2$group)
y2$group_1617 = ifelse(y2$height < 1, "7", as.character(y2$group))
######################

atree = as.treedata(y2)

#MRCA(atree, clades_mrca[1], clades_mrca[2])
#parent(y2, 285)
#parent(y2, 283)
#child(y2, 283)
#child(y2, 560)
#child(y2, 561)
#child(y2, 562)
#child(y2, 563)
#child(y2, 56)


to_drop = c(130, 132, 140, 142, 143)
btree = drop.tip(atree, to_drop)
###########################


p1617 = ggtree(btree, aes(color=group_1617), mrsd = '2018-09-29') +
  scale_color_manual(values = c("#444444", 
                                "#F8766D", "#CD9600",
                                "#00BFC4", "#C77CFF",
                                "#FF61CC", "#dddddd"
  ),
  guide="none")  +
  geom_tippoint(aes(fill=ha_clade_1617), shape=21,  size=4, color="transparent") +
  scale_fill_manual(values = c("#F8766D", "#CD9600",
                               "#7CAE00", "#00BE67", "#00BFC4", 
                               "#6666FF", "#C77CFF", "#FF61CC", "white"),
                    labels = c("3C.3a", "3C.2a", "3C.2a1-1", "3C.2a1-2", "3C.2a1-3",
                               "3C.2a2-1", "3C.2a2-2", "3C.2a3", ""),
                    name = "",
                    na.value=NA) +
  theme_tree2(panel.grid.major.x=element_line(size=0.5, linetype="solid", color="light gray"),
              axis.text.x = element_text(size = 14)) +
  scale_x_continuous(limits = c(2011, 2019), breaks=seq(2011, 2018), labels = seq(2011, 2018)) +
  geom_text(aes(label=clade_label), hjust=1, vjust=-0.3, fontface="bold", size = 6) +
  coord_cartesian(clip = "off")
p1617 = flip(p1617, 293, 381)
p1617

saveRDS(p1617, "ha_tree_1617.rds")
ggsave("ha_tree_1617.png", height=5, width=5)
ggsave("ha_tree_1617.pdf", height=5, width=5)

p1617 = p1617 + scale_fill_manual(values = c("#F8766D", "#CD9600",
                                                     "#7CAE00", "#00BE67", "#00BFC4", 
                                                     "#6666FF", "#C77CFF", "#FF61CC", "white"),
                                  labels = c("3C3.A", "3C2.A", "3C.2a1-1", "3C.2a1-2", "3C.2a1-3",
                                             "3C.2a2-1", "3C.2a2-2", "3C.2a3", ""),
                  name = "",
                  na.value=NA,
                  guide="none")
ggsave("ha_tree_1617_nolegend.png", height=5, width=5)
ggsave("ha_tree_1617_nolegend.pdf", height=5, width=5)



p_tree = ggtree(btree, aes(color=group), mrsd = '2018-09-29') +
  scale_color_manual(values = c("#444444", 
                                "#F8766D", "#CD9600",
                                "#00BFC4", "#C77CFF",
                                 "#FF61CC"
                                ),
                     guide="none")  +
  geom_tippoint(aes(fill=ha_clade_1718), shape=21,  size=4, color="transparent") +
  scale_fill_manual(values = c("#F8766D",  
                                "#00BE67", "#00BFC4", 
                               "#C77CFF", "#FF61CC"),
                    labels = c("3C.3a",  "3C.2a1-2", "3C.2a1-3",
                                "3C.2a2-2", "3C.2a3"),
                    name = "",
                    na.value=NA) +
  theme_tree2(panel.grid.major.x=element_line(size=0.5, linetype="solid", color="light gray"),
              axis.text.x = element_text(size = 14)) +
  scale_x_continuous(limits = c(2011, 2019), breaks=seq(2011, 2019), labels = seq(2011, 2019)) +
  geom_text(aes(label=clade_label), hjust=1, vjust=-0.3, fontface="bold", size = 6)

#check MRCA of 3C2.A1 and 3c2.A2. 
#Node number of 3C2.A1 and 3C2.A2 clade root has to be checked based on clade label
#MRCA(btree, 293, 383)
#child(btree, 292)

p_tree = flip(p_tree, 293, 381)


#p_tree = p_tree + 
#  geom_text(aes(label=clade_label), hjust=1, vjust=-0.3) 

p_tree

saveRDS(p_tree, "ha_tree_1718.rds")
ggsave("ha_tree_1718.png", height=5, width=5)
ggsave("ha_tree_1718.pdf", height=5, width=5)

p_tree + scale_fill_manual(values = c("#F8766D",  
                                      "#00BE67", "#00BFC4", 
                                      "#C77CFF", "#FF61CC"),
                           labels = c("3C.3a",  "3c.2A1-2", "3c.2A1-3",
                                      "3c.2A2-2", "3c.2A3"),
                           name = "",
                           na.value=NA,
                           guide="none")

ggsave("ha_tree_1718_nolegend.png", height=5, width=5)
ggsave("ha_tree_1718_nolegend.pdf", height=5, width=5)

