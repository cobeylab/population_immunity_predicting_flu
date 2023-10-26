library(dplyr)

# opt_na = 1: for NA
# opt_na = 0: for HA

opt_na = 1

age_bin = data.frame()

for(oAG in 1:5){
  opt_ageGroup = oAG

if(opt_na==0 ) {
    source("../../../Data/script/load_data/load_data.R")
  
    sera = sera[,c("Sample_ID", "Age", "3C3.A")]
    
    if (opt_ageGroup > 0){
      sera = transform(sera, Age_group = cut(Age, breaks = c(0, 4, 17, 44, 64, 90)))
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
  
  
  
  tb_ns = sera %>%
    group_by(Age) %>%
    summarise(num_sample = length(Age))
  
  
  ns = (tb_ns$num_sample)
  idx = 1
  idxs = c()
  bin = c()
  for (i in 1:length(ns)) {
    bin = c(bin, ns[i])
    
    if(sum(bin) >= 8){
      
      idxs = c(idxs, idx)
      
      #reset
      bin = c()
      idx = idx + 1
      
    }else{
      
      idxs = c(idxs, idx)
      
      
    }
    
  }
  
  
  last_idx = idxs[length(idxs)]
  last_bin = which(idxs == last_idx)
  if ( sum(ns[last_bin]) < 8 ) {
    idxs[last_bin] = last_idx - 1
  }
  
  age_bin = rbind(age_bin, data.frame(tb_ns$Age,  ns, opt_ageGroup*20+idxs))
  
}

colnames(age_bin) = c("Age", "num_samp", "bin")


if(opt_na == 0){
  saveRDS(age_bin, "age_bin_ha.rds")
}else if(opt_na == 1){
  saveRDS(age_bin, "age_bin_na.rds")
}