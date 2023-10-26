
#################################################################
# Load data
#################################################################

library(reshape2)


# Load sera

if(opt_na==0 ) {
  source("../../../Data/script/load_data/load_data.R")
  
  if (opt_ageGroup > 0){
    sera = transform(sera, Age_group = cut(Age, breaks = c(0, 4, 17, 44, 64, 90)))
    # transform function add X in front of column name if it starts with number, 
    # so remove this X 
    colnames(sera) = gsub("X", "", colnames(sera)) 
    ag = levels(sera$Age_group)
    sera = sera[sera$Age_group == ag[opt_ageGroup], ]
  }
  
}else if (opt_na==1) {
  
  source("../../../Data/script/load_data/load_na.R")
  sera = ella
  
  if (opt_ageGroup > 0){
    sera = transform(sera, Age_group = cut(Age, breaks = c(0, 4, 17, 44, 64, 90)))
    colnames(sera) = gsub("X", "", colnames(sera))
    ag = levels(sera$Age_group)
    sera = sera[sera$Age_group == ag[opt_ageGroup], ]
  }
  
}
  

if(opt_na == 0){
  age_bin = readRDS("age_bin_ha.rds")
  
}else if (opt_na == 1) {
  age_bin = readRDS("age_bin_na.rds")
}


sera = merge(sera, age_bin, by="Age")

# Pop : population size for each age in 2017


pop = read.csv("../data/np2014_d1.csv", sep=",", head=T)
pop = pop[pop$origin == 0 & pop$race == 0 & pop$sex == 0, ]
pop = pop[pop$year == 2017,]

i_minAge = which(colnames(pop) == "pop_1")
i_maxAge = which(colnames(pop) == "pop_90")
pop = pop[,c(i_minAge : i_maxAge)]

pop = data.frame(t(as.matrix(pop)))
colnames(pop) = "n"
pop$Age = as.numeric(gsub("pop_", "", rownames(pop)))

if(group == c("bin")){
  
  pop = merge(pop, age_bin, by="Age")
  
  pop = pop %>%
    group_by(bin) %>%
    summarise(n = sum(n))
  
  pop = data.frame(pop)
}

#################################################################
# Calculate susceptibility
#################################################################



source("utils_susceptibility.R")

# bootstrapping

bootstrapping = function(sera, group, viruses, threshold) {
  susc_boot = data.frame()
  
  if(threshold > 0){

    for (i in 1:n_repl){
      susc1 = susceptibility_frac_unprotected_boot(sera, group, viruses, threshold)
      susc1 = data.frame(susc1)
      susc_boot = rbind(susc_boot, susc1)
    }
  }else if (threshold == 0){
    
    for (i in 1:n_repl){
      susc1 = susceptibility_gmt_boot(sera, group, viruses, threshold)
      susc1 = data.frame(susc1)
      susc_boot = rbind(susc_boot, susc1)
    }
  }

  return (susc_boot)
}

# susceptibility difference and rank 
calc_difference_and_significance = function(susceptibility, susc_boot){
  
  diffs = data.frame()
  for (i in 1:length(susceptibility)) {
    for (j in 1:length(susceptibility)) {
      clade = names(susceptibility)[i]
      compare = names(susceptibility)[j]
      
      diff = susceptibility[i] - susceptibility[j]
  
      diff_boot = susc_boot[,i] - susc_boot[,j]
      diff_boot = sort(diff_boot)
      diff_0052 = diff_boot[ 0.95*length(diff_boot) ]
    
      sig = ifelse( diff > diff_0052 - diff, 1, 0)[[1]]
  
      diffs = rbind(diffs, data.frame(clade, compare, diff, diff_0052, sig) )
    }
  }
  return (diffs)
}

calc_difference_and_significance_pval = function(susceptibility, susc_boot){
  
  diffs = data.frame()
  for (i in 1:length(susceptibility)) {
    for (j in 1:length(susceptibility)) {
      clade = names(susceptibility)[i]
      compare = names(susceptibility)[j]
      
      diff = susceptibility[i] - susceptibility[j]
      
      diff_boot = susc_boot[,i] - susc_boot[,j]
      diff_boot = sort(diff_boot)
      
      diff_0052 = diff_boot[ 0.95*length(diff_boot) ]
      diff_0012 = diff_boot[ 0.99*length(diff_boot) ]
      diff_00012 = diff_boot[ 0.999*length(diff_boot) ]
      
      
      sig_005 = ifelse( diff > diff_0052 - diff, 1, 0)[[1]]
      sig_001 = ifelse( diff > diff_0012 - diff, 1, 0)[[1]]
      sig_0001 = ifelse( diff > diff_00012 - diff, 1, 0)[[1]]
      
      sig = sig_005
      
      diffs = rbind(diffs, data.frame(clade, compare, diff, 
                                      diff_0052, diff_0012, diff_00012, 
                                      sig_005, sig_001, sig_0001, sig) )
    }
  }
  return (diffs)
}

# add rank to df_susc

add_rank = function(df_susc, diffs){
  
  num_sig = c()
  for (v in 1:length(viruses)) {
    diffs_sub = diffs[diffs$clade == viruses[v], ]
    num_sig = c( num_sig, sum(diffs_sub$sig) )
  }
  df_susc$num_sig = num_sig
  df_susc$rank = as.integer( rank (-df_susc$num_sig) )
  
  
  df_susc$Clade = rownames(df_susc)
  df_susc$rank = as.factor(df_susc$rank)
  df_susc$Clade = factor(df_susc$Clade, levels = viruses)
  
  return (df_susc)
}


###################################################################
# relative inferred susceptibility to "susceptibility"

if (threshold > 0){

  susceptibility = susceptibility_frac_unprotected(sera, group, viruses, threshold)
  
}else if (threshold == 0){

  susceptibility = susceptibility_gmt(sera, group, viruses, threshold)
  
}

# bootsrap values to "susc_boot"
susc_boot = bootstrapping(sera, group, viruses, threshold)


#ci
ci1 = sapply(susc_boot, get_95CI1)
ci2 = sapply(susc_boot, get_95CI2)

# df_susc : colnames = susceptibility, ci, ci2, num_sig, rank, clade
df_susc = data.frame(susceptibility, ci1, ci2)

diffs = calc_difference_and_significance_pval(susceptibility, susc_boot)

df_susc = add_rank(df_susc, diffs)



if(opt_na == 0){
  fname = paste0(result_dir, "df_susc", threshold, "_", opt_ageGroup,  ".rds")
  
}else if(opt_na==1){
  fname = paste0(result_dir, "df_susc_na_", threshold, "_", opt_ageGroup,  ".rds")

}
saveRDS(df_susc, fname)

#write.csv(diffs, gsub("rds", "csv", fname), row.names=T)

