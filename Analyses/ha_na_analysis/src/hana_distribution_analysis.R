library(dplyr)

############################################################################
# Load data
############################################################################

source("../../../Data/script/load_data/load_data.R")
source("../../../Data/script/load_data/load_na.R")

ha = sera
head(ha)
colnames(ha)[3:ncol(ha)] = paste0( colnames(ha)[3:ncol(ha)] , "_ha"  )

na = ella
head(na)
colnames(na)[3:ncol(na)] = paste0( colnames(na)[3:ncol(na)] , "_na"  )

hana_indiv = merge(ha, na)

################################################################
# data frame of fraction detectable per age 
################################################################

#detach(package:plyr)

# For each individual, add indicators whether they have HA, NA, have only HA, have only NA,
# have non, and have both
# hana_indiv after adding these indicators will be used for logistic regression

detection_threshold = 1 # By default, consider detectable all titers coded as >=1

hana_indiv <- as_tibble(hana_indiv) %>%
  mutate(haveHA = if_any(.cols = matches('_ha'), function(x){ifelse(is.na(x), F, x >= detection_threshold)}),
         haveNA = if_any(.cols = matches('_na'), function(x){ifelse(is.na(x), F, x >= detection_threshold)}),
         onlyHA = haveHA & !haveNA,
         onlyNA = haveNA & !haveHA,
         none = !haveHA & !haveNA,
         both = haveHA & haveNA)

# This is for making hana_ag_m,
# to make plot of fraction of individuals with detectable antibody by age group for main figure

hana_indiv = transform(hana_indiv, Age_group = cut(Age, breaks = c(0, 10, 20,30,40,50,60,70,90)))
hana_indiv = transform(hana_indiv, Age_group2 = cut(Age, breaks = c(0, 20, 40, 65, 90)))
hana_indiv = transform(hana_indiv, Age_group3 = cut(Age, breaks = as.vector(seq(0, 90, 5))) )


hana_ag = hana_indiv %>%
  group_by(Age_group) %>%
  summarise(frac_none = sum(none)/length(none),
            frac_onlyHA = sum(onlyHA)/length(onlyHA),
            frac_onlyNA = sum(onlyNA)/length(onlyNA),
            frac_both = sum(both)/length(both))
hana_ag = data.frame(hana_ag)


hana_ag$Age_group = gsub("\\(", "", hana_ag$Age_group)
hana_ag$Age_group = gsub("\\]", "", hana_ag$Age_group)

cc = strsplit(hana_ag$Age_group, ",")
part1 = unlist(cc)[2*(1:nrow(hana_ag))-1]
part2 = unlist(cc)[2*(1:nrow(hana_ag)) ]
part1 = as.numeric(part1)+1
hana_ag$Age_group = paste0(part1, "-", part2)


hana_ag_count = hana_indiv %>%
  group_by(Age_group) %>%
  summarise(none = sum(none),
            onlyHA = sum(onlyHA),
            onlyNA = sum(onlyNA),
            both = sum(both))
hana_ag_count = data.frame(hana_ag_count)

hana_ag_count$Age_group = hana_ag$Age_group

#hana_ag_m: for making plot of fraction of individuals with detectable antibody by age group
hana_ag_m = melt(hana_ag, id.vars = c("Age_group"))
colnames(hana_ag_m) = c("Age_group", "Group", "Fraction")

#############################################################
# Andibody distribution by age group 

library(ggplot2) 

theme_detect = theme_bw() +
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title = element_text(size=9),
        
        legend.text = element_text(size=7),
        legend.key.size =  unit(0.3, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-3,-3,-3,-3),
        
        plot.margin = unit(c(1.2, 0.3, 0.3, 0.3), "line"),
        
        text = element_text(size=9),
        panel.grid = element_blank())


ggplot(hana_ag_m, aes(x=Age_group)) +
  geom_bar(aes(y=Fraction, fill=Group), stat="identity") +
  scale_fill_manual( labels = c("No detectable antibody", "HA only", "NA only", "HA and NA"),
                     name = "",
                     values = c("black",  "#56B4E9", "#E69F00", "#009E73")) +
  xlab("Age group (years)") +
  theme_bw() +
  theme(text = element_text(size=9),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title = element_text(size=9),
        
        legend.text = element_text(size=7),
        legend.key.size =  unit(0.4, "cm"),
        legend.margin=margin(0,0,0,0),
        
        panel.grid = element_blank())



ggsave("../figure/hana_ageDist_real_allHA.png", width=4.7, height=2.7)



########################################################################
# Logistic regression using hana_indiv

regression_results = data.frame()

g1 = glm(haveHA ~ Age, data = hana_indiv, family="binomial")
g1 = summary(g1)
regression_results = rbind(regression_results, data.frame(g1$coefficients) )

g2 = glm(haveNA ~ Age, data = hana_indiv, family="binomial")
g2 = summary(g2)
regression_results = rbind(regression_results, data.frame(g2$coefficients) )

g3 = glm(onlyHA ~ Age, data = hana_indiv, family="binomial")
g3 = summary(g3)
regression_results = rbind(regression_results, data.frame(g3$coefficients) )

g4 = glm(onlyNA ~ Age, data = hana_indiv, family="binomial")
g4 = summary(g4)
regression_results = rbind(regression_results, data.frame(g4$coefficients) )


write.csv(regression_results, "../stat_test/reg_age_real_all.csv")

#################################################################################
# looking at all points would help?

p_ha = ggplot(hana_indiv, aes(x=Age, y=haveHA)) +
  geom_jitter(height=0.01, alpha=0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  scale_y_continuous(breaks = c(0, 1)) +
  ylab("Have HA titer") +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90))

p_na = ggplot(hana_indiv, aes(x=Age, y=haveNA)) +
  geom_jitter(height=0.01, alpha=0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  scale_y_continuous(breaks = c(0, 1)) +
  ylab("Have NA titer") +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90))


p_ha_only = ggplot(hana_indiv, aes(x=Age, y=onlyHA)) +
  geom_jitter(height=0.01, alpha=0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  scale_y_continuous(breaks = c(0, 1)) +
  ylab("Have HA titer only") +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90))


p_na_only = ggplot(hana_indiv, aes(x=Age, y=onlyNA)) +
  geom_jitter(height=0.01, alpha=0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  scale_y_continuous(breaks = c(0, 1)) +
  ylab("Have NA titer only") +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90))


ggarrange(p_ha, p_ha_only, p_na, p_na_only,
          nrow=2, ncol=2, 
          labels = c("A", "B", "C", "D"),
          label.y = 1.02)

#ggsave("../figure/fraction_ha_na_real.png", height=4.2, width=5.2)


###############################################################################

#See if results change after removing 1-4 years old
#Results hold: haveHA increase by age, NAonly decrease by age 

# hana_indiv_exp: remove 1-4 years old, because many of them are unexposed

hana_indiv_exp = hana_indiv[hana_indiv$Age > 4,]
hana_indiv_exp = hana_indiv_exp[order(hana_indiv_exp$Age),]


regression_results_exp = data.frame()

g1 = glm(haveHA ~ Age, data = hana_indiv_exp, family="binomial")
g1 = summary(g1)
regression_results_exp = rbind(regression_results_exp, data.frame(g1$coefficients) )

g2 = glm(haveNA ~ Age, data = hana_indiv_exp, family="binomial")
g2 = summary(g2)
regression_results_exp = rbind(regression_results_exp, data.frame(g2$coefficients) )

g3 = glm(onlyHA ~ Age, data = hana_indiv_exp, family="binomial")
g3 = summary(g3)
regression_results_exp = rbind(regression_results_exp, data.frame(g3$coefficients) )

g4 = glm(onlyNA ~ Age, data = hana_indiv_exp, family="binomial")
g4 = summary(g4)
regression_results_exp = rbind(regression_results_exp, data.frame(g4$coefficients) )


#write.csv(regression_results_exp, "../stat_test/reg_age_real_all_exp.csv")


p_ha = ggplot(hana_indiv_exp, aes(x=Age, y=haveHA)) +
  geom_jitter(height=0.01, alpha=0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  scale_y_continuous(breaks = c(0, 1)) +
  ylab("Have HA titer") +
  scale_x_continuous(limits = c(0, 90), breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90))

p_na = ggplot(hana_indiv_exp, aes(x=Age, y=haveNA)) +
  geom_jitter(height=0.01, alpha=0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  scale_y_continuous(breaks = c(0, 1)) +
  ylab("Have NA titer") +
  scale_x_continuous(limits = c(0, 90), breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90))


p_ha_only = ggplot(hana_indiv_exp, aes(x=Age, y=onlyHA)) +
  geom_jitter(height=0.01, alpha=0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  scale_y_continuous(breaks = c(0, 1)) +
  ylab("Have HA titer only") +
  scale_x_continuous(limits = c(0, 90), breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90))


p_na_only = ggplot(hana_indiv_exp, aes(x=Age, y=onlyNA)) +
  geom_jitter(height=0.01, alpha=0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  scale_y_continuous(breaks = c(0, 1)) +
  ylab("Have NA titer only") +
  scale_x_continuous(limits = c(0, 90), breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90))



ggarrange(p_ha, p_ha_only, p_na, p_na_only,
          nrow=2, ncol=2, 
          labels = c("A", "B", "C", "D"),
          label.y = 1.02)

ggsave("../figure/fraction_ha_na_real_exp.png", height=4.2, width=5.2)



####################################################################
# 5. shifting from ha and na to na only


ag = unique(hana_ag$Age_group)
hana_ag$others = hana_ag$frac_none + hana_ag$frac_onlyHA
hana_ag_count$others = hana_ag_count$none + hana_ag_count$onlyHA

cts = c()
for (i in 1:length(ag)){
  for(j in 1:length(ag) ){
    print (paste(i, j))
    ct = chisq.test( rbind( as.numeric( hana_ag_count[j,c(4:5)] ), 
                     as.numeric( hana_ag_count[i,c(4:5)] ) ) )
    print (ct$p.value )
    cts = rbind(cts, data.frame(i, j, ct$p.value))
  }
}


write.csv(cts, "../stat_test/chisqr_real_all.csv")
  


##########################################################################3



hana_age = hana_indiv %>%
  group_by(Age) %>%
  summarise(frac_haveHA = sum(haveHA)/length(haveHA),
            frac_haveNA = sum(haveNA)/length(haveNA),
            frac_none = sum(none)/length(none),
            frac_onlyHA = sum(onlyHA)/length(onlyHA),
            frac_onlyNA = sum(onlyNA)/length(onlyNA),
            frac_haveBoth = sum(both)/length(both),
            num_sample = length(haveHA))

p_ha = ggplot(hana_age, aes(x=Age)) +
  geom_bar(aes(y=frac_haveHA),stat="identity", alpha=0.7) +
  scale_color_manual(labels = c("Similarity", "Imprinting to H3", "Vaccination rate"),
                     values = c("red", "blue", "magenta"),
                     name = "") +
  ylab("Fraction of individuals with \n detectable H3 HA titer") +
  xlab("Age (years)") +
  theme_detect +
  theme(legend.position = "none")

print (p_ha)


# 2. fraction NA by age

p_na =  ggplot(hana_age, aes(x=Age)) +
  geom_bar(aes(y=frac_haveNA),stat="identity", alpha=0.7) +
  scale_color_manual(labels = c("Similarity", "Imprinting to N2", "Vaccination rate"),
                     values = c("red", "blue", "magenta"),
                     name = "") +
  ylab("Fraction of individuals with \n detectable N2 NA titer") +
  xlab("Age (years)") +
  theme_detect +
  theme(legend.position = "none")

print (p_na)


# 3. fraction HA only by age

p_ha_only = ggplot(hana_age, aes(x=Age)) +
  geom_bar(aes(y=frac_onlyHA),stat="identity", alpha=0.7) +
  scale_color_manual(labels = c("Similarity", "Imprinting to H3", "Vaccination rate"),
                     values = c("red", "blue", "magenta"),
                     name = "") +
  ylab("Fraction of individuals with \n detectable H3 HA only") +
  xlab("Age (years)") +
  theme_detect 
print (p_ha_only)



# 4. fraction NA only by age
p_na_only =  ggplot(hana_age, aes(x=Age)) +
  geom_bar(aes(y=frac_onlyNA),stat="identity", alpha=0.7) +
  scale_color_manual(labels = c("Similarity", "Imprinting to N2", "Vaccination rate"),
                     values = c("red", "blue", "magenta"),
                     name = "") +
  ylab("Fraction of individuals with \n detectable N2 NA only")+
  xlab("Age (years)") +
  theme_detect

print (p_na_only)


library(ggpubr)
ggarrange(p_ha, p_ha_only, p_na, p_na_only,
          nrow=2, ncol=2, 
          labels = c("A", "B", "C", "D"),
          label.y = 1.02)

ggsave("../figure/fraction_ha_na_real_byAge.png", height=4.2, width=5.2)
