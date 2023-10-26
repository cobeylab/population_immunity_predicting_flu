library(ggplot2)
library(dplyr)
library(reshape2)
source("analysis_util.R")

ag_brk = c(0, 4, 17, 44, 64, 90)

##############################################################
# Data: ag_sera (age group + HA titers)
################################################################


source("../../../Data/script/load_data/load_data.R")

# levels for viruses used for correlation analysis

vlevels = c("X3C3.A", "X3C2.A", "N121K_S144K",
            "N171K",  "N121K_N171K", "N121K_T135K_N171K",
            "T131K_R142K", "T131K_R142K_R261Q")

#observed data

sera0 = sera

sera0 = transform(sera0, Age_group = cut(Age, breaks = c(0, 4, 17, 44, 64, 90)))

sera0$label = sera0$Age_group
sera0$label = gsub("\\(0,4\\]", "1-4", sera0$label)
sera0$label = gsub("\\(4,17\\]", "5-17", sera0$label)
sera0$label = gsub( "\\(17,44\\]", "18-44", sera0$label)
sera0$label = gsub( "\\(44,64\\]", "45-64", sera0$label)
sera0$label = gsub("\\(64,90\\]", "65-90", sera0$label)

sera0$label = factor(sera0$label, levels = c("1-4", "5-17", "18-44", "45-64", "65-90"))

# finally, make data suitable for analysis
# age group + titers

ag_sera = sera0[, c("Age_group", vlevels)]

  
################################################################

#
rs = rowSums(ag_sera[ , c(2:ncol(ag_sera))], na.rm = T)
ag_sera = ag_sera[rs != 0, ]

##################################################################
#1. no scale

source("correlation_functions.R")

ag_cormat = get_ag_cormat_rp(ag_sera, "spearman")
ag_cormat$label = ag_cormat$Age_group 
ag_cormat$label = gsub("\\(0,4\\]", "1-4", ag_cormat$label)
ag_cormat$label = gsub("\\(4,17\\]", "5-17", ag_cormat$label)
ag_cormat$label = gsub( "\\(17,44\\]", "18-44", ag_cormat$label)
ag_cormat$label = gsub( "\\(44,64\\]", "45-64", ag_cormat$label)
ag_cormat$label = gsub("\\(64,90\\]", "65-90", ag_cormat$label)
ag_cormat$label = factor(ag_cormat$label, levels = c("1-4", "5-17", "18-44", "45-64", "65-90"))


# 2. correlation difference
source("correlation_function_bootstrap.R")

ag = unique(ag_sera$Age_group)
ag = sort(ag)


ag_cormat_repl = list()

n_repl=1000
for(r in 1:n_repl){
  
  ag_sera_resample = get_resamples(ag_sera, ag)
  ag_cormat_1 = get_ag_cormat_rp(ag_sera_resample, "spearman")
  ag_cormat_repl[[r]] = data.frame(ag_cormat_1)  
}


# p value from testing if the correlation is weaker than the other

#age group 1-4

pdf_all_1 = data.frame(i=numeric(0), j=numeric(0), ks = numeric(0), ms = numeric(0), ps = numeric(0) )

for(a in 1:1){
  
  #For each age group
  
  for (i in 1:length(vlevels)){
    print (i)
    for (j in (i:length(vlevels))){
      print (j)
      #compare if correlation coefficient of each pair is 
      # significantly different from those of all other pairs
      if (i == j) {
        next
      }else{
        ps = get_num_sig_pval_within(ag_cormat_repl, ag_cormat, ag, vlevels, n_repl, a, i, j)
        if(nrow(ps) == 0) {
          next
        }
        pdf = data.frame(i, j, ps)
        pdf_all_1 = rbind(pdf_all_1, pdf)
      }
      
    }  
  }
}
#saveRDS(pdf_all_1, "../result/pdf_1.rds")
#write.csv(pdf_all_1, "../result/pdf_1.csv", row.names=T)


#age group 5-17

pdf_all_2 = data.frame(i=numeric(0), j=numeric(0), ks = numeric(0), ms = numeric(0), ps = numeric(0) )

for(a in 2:2){
  
  #For each age group
  
  for (i in 1:length(vlevels)){
    print (i)
    for (j in (i:length(vlevels))){
      print (j)
      #compare if correlation coefficient of each pair is 
      # significantly different from those of all other pairs
      if (i == j) {
        next
      }else{
        ps = get_num_sig_pval_within(ag_cormat_repl, ag_cormat, ag, vlevels, n_repl, a, i, j)
        if(nrow(ps) == 0) {
          next
        }
        pdf = data.frame(i, j, ps)
        pdf_all_2 = rbind(pdf_all_2, pdf)
      }
      
    }  
  }
}
#saveRDS(pdf_all_2, "../result/pdf_2.rds")
#write.csv(pdf_all_2, "../result/pdf_2.csv", row.names=T)


#age group 18-44

pdf_all_3 = data.frame(i=numeric(0), j=numeric(0), ks = numeric(0), ms = numeric(0), ps = numeric(0) )

for(a in 3:3){
  
  #For each age group
  
  for (i in 1:length(vlevels)){
    print (i)
    for (j in (i:length(vlevels))){
      print (j)
      #compare if correlation coefficient of each pair is 
      # significantly different from those of all other pairs
      if (i == j) {
        next
      }else{
        ps = get_num_sig_pval_within(ag_cormat_repl, ag_cormat, ag, vlevels, n_repl, a, i, j)
        if(nrow(ps) == 0) {
          next
        }
        pdf = data.frame(i, j, ps)
        pdf_all_3 = rbind(pdf_all_3, pdf)
      }
      
    }  
  }
}
#saveRDS(pdf_all_3, "../result/pdf_33.rds")
#write.csv(pdf_all_3, "../result/pdf_33.csv", row.names=T)


#age group 45-64

pdf_all_4 = data.frame(i=numeric(0), j=numeric(0), ks = numeric(0), ms = numeric(0), ps = numeric(0) )

for(a in 4:4){
  
  #For each age group
  
  for (i in 1:length(vlevels)){
    print (i)
    for (j in (i:length(vlevels))){
      print (j)
      #compare if correlation coefficient of each pair is 
      # significantly different from those of all other pairs
      if (i == j) {
        next
      }else{
        ps = get_num_sig_pval_within(ag_cormat_repl, ag_cormat, ag, vlevels, n_repl, a, i, j)
        if(nrow(ps) == 0) {
          next
        }
        pdf = data.frame(i, j, ps)
        pdf_all_4 = rbind(pdf_all_4, pdf)
      }
      
    }  
  }
}
#saveRDS(pdf_all_4, "../result/pdf_4.rds")
#write.csv(pdf_all_4, "../result/pdf_4.csv", row.names=T)


#age group 65-90

pdf_all_5 = data.frame(i=numeric(0), j=numeric(0), ks = numeric(0), ms = numeric(0), ps = numeric(0) )

for(a in 5:5){
  
  #For each age group
  
  for (i in 1:length(vlevels)){
    print (i)
    for (j in (i:length(vlevels))){
      print (j)
      #compare if correlation coefficient of each pair is 
      # significantly different from those of all other pairs
      if (i == j) {
        next
      }else{
        ps = get_num_sig_pval_within(ag_cormat_repl, ag_cormat, ag, vlevels, n_repl, a, i, j)
        if(nrow(ps) == 0) {
          next
        }
        pdf = data.frame(i, j, ps)
        pdf_all_5 = rbind(pdf_all_5, pdf)
      }
      
    }  
  }
}
#saveRDS(pdf_all_5, "../result/pdf_5.rds")
#write.csv(pdf_all_5, "../result/pdf_5.csv", row.names=T)



