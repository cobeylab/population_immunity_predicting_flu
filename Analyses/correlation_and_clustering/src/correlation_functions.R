#library(cocor)


get_cormat = function(sera, cols, method){
  #sera: matrix of titers, each virus at each column
  cormat = cor(sera[,cols], method = method, use="pairwise.complete.obs")
  cormat[upper.tri(cormat)] = NA
  m_cormat = melt(cormat, na.rm=T)
  m_cormat
}

get_ag_cormat = function(ag_sera, method){
  #sera: age group in the first column and then titers by virus for each column
  ag_cormat = ag_sera %>%
    group_by(Age_group) %>%
    do(get_cormat(., c(2:ncol(ag_sera)), method))
  
  ag_cormat
}

plot_cor_heatmap = function(ag_cormat, title){
  #ag_cormat : output of get_ag_cormat
  p = ggplot(ag_cormat, aes(Var1, Var2, fill=value)) +
    geom_tile(col="white") +
    scale_fill_gradient2(low="white", mid = "light blue", high="blue",  
                         midpoint = 0.8, name="Correlation coefficient (r)") +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    coord_fixed() + 
    facet_wrap(~Age_group) +
    ylab("") + xlab("") +
    ggtitle(title)
  
  print (p)
}

scale_each_row = function(sera) {
  #sera: matrix of titers, each virus at each column
  tz = c()
  for(i in 1:nrow(sera)) {
    row1 = sera[i,]
    if(sd(row1, na.rm=T) == 0) {
      row1 = rep(0, 8)
    }else{
      row1 = t( scale(t(row1), center=T, scale=T ) )
    }
    tz = rbind(tz, row1)
  }
  tz
}

cor_test_rp = function(ag_sera, cols, method){
  
  titers = ag_sera[,cols]
  cortest = c()
  
  for(i in c(1: (ncol(titers))) ) {
    for(j in c( (i) : (ncol(titers)) )){
      
      ct = cor.test(titers[[i]], as.vector(titers[[j]]), method=method)
      ct1 = c(colnames(titers)[i], colnames(titers)[j], ct$estimate, ct$p.value)
      cortest = rbind(cortest, ct1)
      
    } 
  }
  colnames(cortest) = c("Virus1", "Virus2", "r", "p")
  
  data.frame(cortest)
}

get_ag_cormat_rp = function(ag_sera, method){
  #sera: age group in the first column and then titers by virus for each column
  ag_cormat = ag_sera %>%
    group_by(Age_group) %>%
    do(cor_test_rp(., c(2:ncol(ag_sera)), method))
  
  ag_cormat$r = as.numeric(ag_cormat$r)
  ag_cormat$p = as.numeric(ag_cormat$p)
  
  ag_cormat$Virus1 = factor(ag_cormat$Virus1, levels = vlevels)
  ag_cormat$Virus2 = factor(ag_cormat$Virus2, levels = rev(vlevels))
  
  ag_cormat
}


cor_test_rp_rmv_pair_undetectable = function(ag_sera, cols, age_group, method){
  
  titers = ag_sera[ag_sera$Age_group == age_group, cols]

  is = c()
  js = c()
  coeffs = c()
  ps = c()
  Age_group = c()

  for(i in c(1: (ncol(titers))) ) {
    for(j in c( (i) : (ncol(titers)) )){
      
      v1 = c()
      v2 = c()
      for (k in 1:nrow(titers)){
        #print (paste("i:", i, " j:", j, " k:", k))
        if(is.na(titers[k,i]) | is.na(titers[k,j])){
          next
        }
        if(titers[k,i] == 0 & titers[k,j]==0){
          next
        }
        v1 = c(v1, titers[k,i])
        v2 = c(v2, titers[k,j])
      }
      
      ccc = cor.test(v1, v2, method = method, use="complete.obs")
      is = c(is, colnames(titers)[i])
      js = c(js, colnames(titers)[j])
      coeffs = c(coeffs, ccc$estimate)
      ps = c(ps, ccc$p.value)
      Age_group = c(Age_group, age_group)
      
    } 
  }
  pair_df = data.frame(Age_group, is, js, coeffs, ps)
  
  pair_df
}

colv = colnames(ag_sera)[2:ncol(ag_sera)]

ff11 = formula(paste("~", colv[1], "+", colv[1], "|", colv[1], "+", colv[1] ))
ff12 = formula(paste("~", colv[1], "+", colv[2], "|", colv[1], "+", colv[2] ))
ff13 = formula(paste("~", colv[1], "+", colv[3], "|", colv[1], "+", colv[3] ))
ff14 = formula(paste("~", colv[1], "+", colv[4], "|", colv[1], "+", colv[4] ))
ff15 = formula(paste("~", colv[1], "+", colv[5], "|", colv[1], "+", colv[5] ))
ff16 = formula(paste("~", colv[1], "+", colv[6], "|", colv[1], "+", colv[6] ))
ff17 = formula(paste("~", colv[1], "+", colv[7], "|", colv[1], "+", colv[7] ))
ff18 = formula(paste("~", colv[1], "+", colv[8], "|", colv[1], "+", colv[8] ))

ff22 = formula(paste("~", colv[2], "+", colv[2], "|", colv[2], "+", colv[2] ))
ff23 = formula(paste("~", colv[2], "+", colv[3], "|", colv[2], "+", colv[3] ))
ff24 = formula(paste("~", colv[2], "+", colv[4], "|", colv[2], "+", colv[4] ))
ff25 = formula(paste("~", colv[2], "+", colv[5], "|", colv[2], "+", colv[5] ))
ff26 = formula(paste("~", colv[2], "+", colv[6], "|", colv[2], "+", colv[6] ))
ff27 = formula(paste("~", colv[2], "+", colv[7], "|", colv[2], "+", colv[7] ))
ff28 = formula(paste("~", colv[2], "+", colv[8], "|", colv[2], "+", colv[8] ))

ff33 = formula(paste("~", colv[3], "+", colv[3], "|", colv[3], "+", colv[3] ))
ff34 = formula(paste("~", colv[3], "+", colv[4], "|", colv[3], "+", colv[4] ))
ff35 = formula(paste("~", colv[3], "+", colv[5], "|", colv[3], "+", colv[5] ))
ff36 = formula(paste("~", colv[3], "+", colv[6], "|", colv[3], "+", colv[6] ))
ff37 = formula(paste("~", colv[3], "+", colv[7], "|", colv[3], "+", colv[7] ))
ff38 = formula(paste("~", colv[3], "+", colv[8], "|", colv[3], "+", colv[8] ))

ff44 = formula(paste("~", colv[4], "+", colv[4], "|", colv[4], "+", colv[4] ))
ff45 = formula(paste("~", colv[4], "+", colv[5], "|", colv[4], "+", colv[5] ))
ff46 = formula(paste("~", colv[4], "+", colv[6], "|", colv[4], "+", colv[6] ))
ff47 = formula(paste("~", colv[4], "+", colv[7], "|", colv[4], "+", colv[7] ))
ff48 = formula(paste("~", colv[4], "+", colv[8], "|", colv[4], "+", colv[8] ))

ff55 = formula(paste("~", colv[5], "+", colv[5], "|", colv[5], "+", colv[5] ))
ff56 = formula(paste("~", colv[5], "+", colv[6], "|", colv[5], "+", colv[6] ))
ff57 = formula(paste("~", colv[5], "+", colv[7], "|", colv[5], "+", colv[7] ))
ff58 = formula(paste("~", colv[5], "+", colv[8], "|", colv[5], "+", colv[8] ))

ff66 = formula(paste("~", colv[6], "+", colv[6], "|", colv[6], "+", colv[6] ))
ff67 = formula(paste("~", colv[6], "+", colv[7], "|", colv[6], "+", colv[7] ))
ff68 = formula(paste("~", colv[6], "+", colv[8], "|", colv[6], "+", colv[8] ))

ff77 = formula(paste("~", colv[7], "+", colv[7], "|", colv[7], "+", colv[7] ))
ff78 = formula(paste("~", colv[7], "+", colv[8], "|", colv[7], "+", colv[8] ))

ff88 = formula(paste("~", colv[8], "+", colv[8], "|", colv[8], "+", colv[8] ))

compare_cor = function(lsera, alt){
  cc = c()

  cc = c(cc, 1)
  cc = c(cc, cocor( ff12, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff13, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff14, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff15, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff16, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff17, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff18, alternative = alt,lsera)@fisher1925$p.value)
  
  cc = c(cc, 1)
  cc = c(cc, cocor( ff23, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff24, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff25, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff26, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff27, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff28, alternative = alt,lsera)@fisher1925$p.value)
 
  cc = c(cc, 1)
  cc = c(cc, cocor( ff34, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff35, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff36, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff37, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff38, alternative = alt,lsera)@fisher1925$p.value)
  
  cc = c(cc, 1)
  cc = c(cc, cocor( ff45, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff46, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff47, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff48, alternative = alt,lsera)@fisher1925$p.value)
  
  cc = c(cc, 1)
  cc = c(cc, cocor( ff56, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff57, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff58, alternative = alt,lsera)@fisher1925$p.value)
  
  cc = c(cc, 1)
  cc = c(cc, cocor( ff67, alternative = alt,lsera)@fisher1925$p.value)
  cc = c(cc, cocor( ff68, alternative = alt,lsera)@fisher1925$p.value)
 
  cc = c(cc, 1)
  cc = c(cc, cocor( ff78, alternative = alt,lsera)@fisher1925$p.value)

  cc = c(cc, 1)
  
  cc
}


compare_cor_within = function(lsera_g, alt, i, j){
  
  cc = c()
  
  for (k in 1:length(colv)) {
    for (m in k:length(colv)){
      if (k == m){
        cc = c(cc, 1)
        
      }else if(i == k & j == m){
        cc = c(cc, 1)
        
      }else if ( (i != k) & (j != m) & (i != m) & (j != k) ) {#non-overlapping variables
        fml = formula(paste("~", colv[i], "+", colv[j], "|", colv[k], "+", colv[m] ))
        cc = c(cc, cocor( fml, alternative = alt, lsera_g)@silver2004$p.value)
        
      }else {#overlapping variables
        fml = formula(paste("~", colv[i], "+", colv[j], "|", colv[k], "+", colv[m] ))
        cc = c(cc, cocor( fml, alternative = alt, lsera_g)@hittner2003$p.value)
      }
    }
  }
  
  cc
}



