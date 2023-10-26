

get_resamples = function(ag_sera, ag){
  
  ag_sera_resample = data.frame()
  for (i in 1:length(ag)){
    
    one_ag = ag_sera[ag_sera$Age_group == ag[i],]
    resampled =  sample.int(nrow(one_ag),  nrow(one_ag), replace=T)
    one_ag_resample = one_ag[resampled, ]
    
    ag_sera_resample = rbind(ag_sera_resample, one_ag_resample)
  }

  return (ag_sera_resample)  
}



get_rho_all = function(ag_sera_resample, ag) {
  
  diffs = list()
  
  diffs_greater = c()
  diffs_less = c()
  
  for(a in 2:length(ag)){
    
    ag1_sera_resample = ag_sera_resample[ag_sera_resample$Age_group == ag[1],]
    ag2_sera_resample = ag_sera_resample[ag_sera_resample$Age_group == ag[a],]
    
    for( i in c( 2: (ncol(ag1_sera_resample)-1)  ) ){
      for( j in c( i:ncol(ag1_sera_resample) ) ) {
        
        ct1 = cor.test(ag1_sera_resample[,i], ag1_sera_resample[,j],
                       method="spearman")
        ct2 = cor.test(ag2_sera_resample[,i], ag2_sera_resample[,j],
                       method="spearman")
        
        diffs_greater = c(diffs_greater, ct2$estimate[[1]] - ct1$estimate[[1]])
        diffs_less = c(diffs_less, ct1$estimate[[1]] - ct2$estimate[[1]])
        
      }
    }
  }
  
  diffs$greater = diffs_greater
  diffs$less = diffs_less
  return(diffs)
}

get_diffs = function(ag_sera_resample, ag) {

  diffs = list()
  
  diffs_greater = c()
  diffs_less = c()
  
  for(a in 2:length(ag)){

    ag1_sera_resample = ag_sera_resample[ag_sera_resample$Age_group == ag[1],]
    ag2_sera_resample = ag_sera_resample[ag_sera_resample$Age_group == ag[a],]
    
    for( i in c( 2: (ncol(ag1_sera_resample)-1)  ) ){
      for( j in c( i:ncol(ag1_sera_resample) ) ) {
        
        ct1 = cor.test(ag1_sera_resample[,i], ag1_sera_resample[,j],
                       method="spearman")
        ct2 = cor.test(ag2_sera_resample[,i], ag2_sera_resample[,j],
                       method="spearman")

        diffs_greater = c(diffs_greater, ct2$estimate[[1]] - ct1$estimate[[1]])
        diffs_less = c(diffs_less, ct1$estimate[[1]] - ct2$estimate[[1]])
        
      }
    }
  }
  
  diffs$greater = diffs_greater
  diffs$less = diffs_less
  return(diffs)
}



get_num_sig_within = function(ag_cormat_repl, ag_cormat, ag, vlevels, n_repl, a, i, j){
  
  num_sig = 0
  
  for (k in 1:length(vlevels)){
    for(m in k:length(vlevels)) {
      
      if(k==m) {
        next
      }
      
      null_dist_less = c()
      for(r in (1:n_repl)){
        
        age_group_1 = ag_cormat_repl[[r]][ag_cormat_repl[[r]]$Age_group == ag[a],]
        
        group_1 = age_group_1[age_group_1$Virus1 == vlevels[i] & age_group_1$Virus2 == vlevels[j],]
        group_2 = age_group_1[age_group_1$Virus1 == vlevels[k] & age_group_1$Virus2 == vlevels[m],]
        
        diff_less = group_2[,"r"] - group_1[,"r"]
        
        null_dist_less = c(null_dist_less, diff_less)
        
      }
      
      ag1 = data.frame( ag_cormat[ag_cormat$Age_group == ag[a], ] )
      
      r1 = ag1[ag1$Virus1 == vlevels[i] & ag1$Virus2 == vlevels[j], "r"]
      r2 = ag1[ag1$Virus1 == vlevels[k] & ag1$Virus2 == vlevels[m], "r"]
      
      null_dist_less = null_dist_less - (r2-r1)
      null_dist_less = sort(null_dist_less)
      
      if(null_dist_less[as.integer( n_repl*0.95 )] < (r2-r1) ){
        num_sig = num_sig + 1
      }
      
    }
  }

  return (num_sig)
  
}



get_num_sig_pval_within = function(ag_cormat_repl, ag_cormat, ag, vlevels, n_repl, a, i, j){
  
  num_sig = 0
  ks = c()
  ms = c()
  ps = c()
  
  for (k in 1:length(vlevels)){
    for(m in k:length(vlevels)) {
      
      if(k==m) {
        next
      }
      
      null_dist_less = c()
      for(r in (1:n_repl)){
        
        age_group_1 = ag_cormat_repl[[r]][ag_cormat_repl[[r]]$Age_group == ag[a],]
        
        group_1 = age_group_1[age_group_1$Virus1 == vlevels[i] & age_group_1$Virus2 == vlevels[j],]
        group_2 = age_group_1[age_group_1$Virus1 == vlevels[k] & age_group_1$Virus2 == vlevels[m],]
        
        diff_less = group_2[,"r"] - group_1[,"r"]
        
        null_dist_less = c(null_dist_less, diff_less)
        
      }
      
      ag1 = data.frame( ag_cormat[ag_cormat$Age_group == ag[a], ] )
      
      r1 = ag1[ag1$Virus1 == vlevels[i] & ag1$Virus2 == vlevels[j], "r"]
      r2 = ag1[ag1$Virus1 == vlevels[k] & ag1$Virus2 == vlevels[m], "r"]
      
      null_dist_less = null_dist_less - (r2-r1)
      null_dist_less = sort(null_dist_less)
      

      
      if(null_dist_less[as.integer( n_repl*0.95 )] < (r2-r1) ){
        
        for (p in (n_repl*0.95):n_repl){
          if (null_dist_less[n_repl] < (r2-r1)) {
            pval = 0
            break
          }
          
          if (null_dist_less[p] < (r2-r1)){
            next
          }else {
            pval = (n_repl - p + 1)/(n_repl)
            #print (paste0 (p, " ", null_dist_less[p], " ", (r2-r1), " ", pval))
            
            break
          }
          
        }
        num_sig = num_sig + 1
        ks = c(ks, k)
        ms = c(ms, m)
        ps = c(ps, pval)
      }
      
    }
  }
  
  df = data.frame(ks, ms, ps)
  return (df)
  
}


