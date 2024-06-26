library(ggplot2)
library(cowplot)
library(dplyr)
library(reshape2)
source("analysis_util.R")

ag_brk = c(0, 4, 17, 44, 64, 90)

##############################################################
# Make data frame "ag_sera" (age group + HA titers)
################################################################

# Read titer data
source("../../../Data/script/load_data/load_data.R")

sera0 = sera

# transform function add 'X' in front of column name if it starts with number
sera0 = transform(sera0, Age_group = cut(Age, breaks = c(0, 4, 17, 44, 64, 90)))


# levels for viruses used for correlation analysis
vlevels = c("X3C3.A", "X3C2.A", "N121K_S144K",
            "N171K",  "N121K_N171K", "N121K_T135K_N171K",
            "T131K_R142K", "T131K_R142K_R261Q")


sera0$label = sera0$Age_group
sera0$label = gsub("\\(0,4\\]", "1-4", sera0$label)
sera0$label = gsub("\\(4,17\\]", "5-17", sera0$label)
sera0$label = gsub( "\\(17,44\\]", "18-44", sera0$label)
sera0$label = gsub( "\\(44,64\\]", "45-64", sera0$label)
sera0$label = gsub("\\(64,90\\]", "65-90", sera0$label)

sera0$label = factor(sera0$label, levels = c("1-4", "5-17", "18-44", "45-64", "65-90"))

# This is a convenience spot to make an age distribution figure for the supplement:

age_dist <- sera0 %>% 
  ggplot(aes(x = Age)) +
  geom_histogram(binwidth = 1) +
  xlab('Age (years)') +
  ylab('Number of serum sample donors') +
  theme_cowplot() +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 11))

age_dist_binned <- sera0 %>% group_by(label) %>% count() %>%
  ggplot(aes(x = label, y = n)) +
  geom_col() +
  geom_text(aes(label = n), position = position_nudge(x = 0, y = 5)) + 
  ylab('') +
  xlab('Age group (years)') +
  theme_cowplot() +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 11))


save_plot("../../sample_age_dist.pdf",
          plot_grid(age_dist, age_dist_binned, nrow = 1),
          base_height = 3.5, base_width = 8)


# Make data frame for analysis
# ag_sera: age group + titers

ag_sera = sera0[, c("Age_group", vlevels)]

##################################################################
#remove all undetectable
  
rs = rowSums(ag_sera[ , c(2:ncol(ag_sera))], na.rm = T)
ag_sera = ag_sera[rs != 0, ]

##################################################################
#1. Calculate correlation coefficient for each pairs within each age group

source("correlation_functions.R")

ag_cormat <- get_ag_cormat(ag_sera, "spearman", n_imputations = 1000)  

ag_cormat$label = ag_cormat$Age_group 
ag_cormat$label = gsub("\\(0,4\\]", "1-4", ag_cormat$label)
ag_cormat$label = gsub("\\(4,17\\]", "5-17", ag_cormat$label)
ag_cormat$label = gsub( "\\(17,44\\]", "18-44", ag_cormat$label)
ag_cormat$label = gsub( "\\(44,64\\]", "45-64", ag_cormat$label)
ag_cormat$label = gsub("\\(64,90\\]", "65-90", ag_cormat$label)
ag_cormat$label = factor(ag_cormat$label, levels = c("1-4", "5-17", "18-44", "45-64", "65-90"))


# 2. Test correlation difference using bootstrapping

source("correlation_function_bootstrap.R")

ag = unique(ag_sera$Age_group)
ag = sort(ag)


ag_cormat_repl = list()

n_repl=1000
# Here we can do a single random imputation each time titers are reshuffled
for(r in 1:n_repl){
  
  ag_sera_resample = get_resamples(ag_sera, ag)
  ag_cormat_1 = get_ag_cormat(ag_sera_resample, "spearman", n_imputations = 1)
  ag_cormat_repl[[r]] = data.frame(ag_cormat_1)  
}



n_virus_pairs <-  ag_cormat %>% filter(Age_group == '(0,4]') %>% count() %>% pull(n)

r_ccg = rep(NA, n_virus_pairs)
r_ccl = rep(NA, n_virus_pairs)

for(a in 2:length(ag)){

  for(i in 1:n_virus_pairs){
    v1 <- ag_cormat %>% filter(Age_group == ag[a]) %>% slice(i) %>% pull(Virus1)
    v2 <- ag_cormat %>% filter(Age_group == ag[a]) %>% slice(i) %>% pull(Virus2)
    
    if(v1 == v2){
      r_ccg = c(r_ccg, NA)
      r_ccl = c(r_ccl, NA) 
    }else{
      null_dist_greater = c()
      null_dist_less = c()
      
      for(r in 1:n_repl){

        age_group_1 = ag_cormat_repl[[r]][ag_cormat_repl[[r]]$Age_group == ag[1],]
        age_group_2 = ag_cormat_repl[[r]][ag_cormat_repl[[r]]$Age_group == ag[a],]
          
        group_1 = age_group_1 %>% filter(Virus1 == v1, Virus2 == v2)
        group_2 = age_group_2 %>% filter(Virus1 == v1, Virus2 == v2)
          
        diff_greater = group_2[,"r"] - group_1[,"r"]
        diff_less = group_1[,"r"] - group_2[,"r"]
        
        null_dist_greater = c(null_dist_greater, diff_greater)
        null_dist_less = c(null_dist_less, diff_less)
      }
      
      ag1 = data.frame( ag_cormat[ag_cormat$Age_group == ag[1], ] )
      ag2 = data.frame( ag_cormat[ag_cormat$Age_group == ag[a], ] )
      
      r1 = ag1 %>% filter(Virus1 == v1, Virus2 == v2) %>% pull(r)
      
      r2 = ag2 %>% filter(Virus1 == v1, Virus2 == v2) %>% pull(r)
      
      null_dist_greater = null_dist_greater - (r2-r1)
      null_dist_less = null_dist_less - (r1-r2)
      
      null_dist_greater = sort(null_dist_greater)
      null_dist_less = sort(null_dist_less)
      
      if(null_dist_greater[as.integer( n_repl*0.95 )] < (r2-r1) ){
        r_ccg = c(r_ccg, r2) 
      }else {
        r_ccg = c(r_ccg, NA) 
      }
      
      if(null_dist_less[as.integer( n_repl*0.95 )] < (r1-r2) ){
        r_ccl = c(r_ccl, r2)
      }else{
        r_ccl = c(r_ccl, NA)
      }
      
    }
  }
}

# 2.1. The age group greater than 1-4 
ag_cormat$r_ccg = r_ccg


# 2.2. The age group greater than 1-4
ag_cormat$r_ccl = r_ccl


##############################################################################

# 3. Within an age group, rank of correlation

ag_cormat <- ag_cormat %>% mutate(rank = NA)
#rank <- c()

for(a in 1:length(ag)){
  
  #For each age group
  
  for (i in 1:length(vlevels)){
    
    for (j in (i:length(vlevels))){
      
      print (paste(ag[a], vlevels[i], vlevels[j]))
      
      #compare if correlation coefficient of each pair is 
      # significantly different from those of all other pairs
      if (i == j) {
        #rank = c(rank, NA)
        num_sig = NA
      }else{
        num_sig = get_num_sig_within(ag_cormat_repl, ag_cormat, ag, vlevels, n_repl, a, i, j)
        #rank = c(rank, num_sig)
      }
      
      ag_cormat[ag_cormat$Age_group == ag[a] &
                  (
                    (ag_cormat$Virus1 == vlevels[i] & ag_cormat$Virus2 == vlevels[j]) |
                    (ag_cormat$Virus1 == vlevels[j] & ag_cormat$Virus2 == vlevels[i])
                    
                  ), "rank"] <- num_sig
    }  
  }
}




#######################################################################

saveRDS(ag_cormat, "../result/ag_cormat_real_rmv_undetectable.rds")


