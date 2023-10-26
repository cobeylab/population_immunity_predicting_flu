
targeting_difference = function(sera) {
  
  ag_brk = c(0,4,17,44,64,90)
  
  sera = transform(sera, Age_group = cut(Age, breaks=ag_brk))
  sera = sera[ !is.na(sera[,2]) * !is.na(sera[,3]), ]
  
  sera$diff = sera[,2] - sera[,3]
  
  #within age group by indiv
  
  age_group = sort(unique(sera$Age_group))
  
  
  
  diff_indiv = c()
  
  for (i in 1:length(age_group)) {
    
    ag1 = sera[sera$Age_group == age_group[i],]
    
    for (j in c(1: (nrow(ag1))-1 )){
      for (k in c( (j+1) : nrow(ag1) )){
        diff_indiv = c(diff_indiv, abs(ag1$diff[j] - ag1$diff[k]))
      }
    }
  }
  
  sum(diff_indiv)/length(diff_indiv)
  
  ########################
  
  diff_ag = c()
  what1 = c()
  what2 = c()
  
  for (i in c( 1: (length(age_group)-1) ) ) {
    ag1 = sera[sera$Age_group == age_group[i],]
    
    for (h in c( (i+1) : length(age_group) )){
      ag2 = sera[sera$Age_group == age_group[h],]
      
      for (j in c(1: (nrow(ag1)) )){
        for (k in c( 1: nrow(ag2) )){
          what1 = c(what1, ag1$diff[j])
          what2 = c(what2, ag2$diff[k])
          diff_ag = c(diff_ag, abs(ag1$diff[j] - ag2$diff[k]))
          
        }
      }
    }
  
  }
  sum(diff_ag)/length(diff_ag)
  
  tt = t.test(diff_indiv, diff_ag)
  return (tt)
  
}
###############################################################

source("../../../Data/script/load_data/load_data.R")
ag_brk = c(0,4,17,44,64,90)

ag_sera = transform(sera, Age_group = cut(Age, breaks=ag_brk))

ag_sera$diff1 =  ag_sera$X3C3.A - ag_sera$X3C2.A
ag_sera$diff2 =  ag_sera$X3C3.A - ag_sera$N171K
ag_sera$diff3 =  ag_sera$X3C3.A - ag_sera$N121K_N171K
ag_sera$diff4 =  ag_sera$X3C3.A - ag_sera$N121K_T135K_N171K
ag_sera$diff5 =  ag_sera$X3C3.A - ag_sera$T131K_R142K
ag_sera$diff6 =  ag_sera$X3C3.A - ag_sera$T131K_R142K_R261Q
ag_sera$diff7 =  ag_sera$X3C3.A - ag_sera$N121K_S144K


#tt1 = oneway.test(diff1 ~ Age_group, ag_sera, var.equal = F)

test1 = aov(diff1 ~ Age_group, ag_sera)
test2 = aov(diff2 ~ Age_group, ag_sera)
test3 = aov(diff3 ~ Age_group, ag_sera)
test4 = aov(diff4 ~ Age_group, ag_sera)
test5 = aov(diff5 ~ Age_group, ag_sera)
test6 = aov(diff6 ~ Age_group, ag_sera)
test7 = aov(diff7 ~ Age_group, ag_sera)

tt1 = summary(test1)
tt2 = summary(test2)
tt3 = summary(test3)
tt4 = summary(test4)#different
tt5 = summary(test5)
tt6 = summary(test6)#different
tt7 = summary(test7)#different

tt = data.frame(tt1[[1]][1, 3], tt1[[1]][2, 3], tt1[[1]][1, 4], tt1[[1]][1, 5]) 
colnames(tt) = c("a", "b", "c", "d")

tt = rbind(tt, data.frame(a=tt2[[1]][1, 3], b=tt2[[1]][2, 3], c=tt2[[1]][1, 4], d=tt2[[1]][1, 5]))
tt = rbind(tt, data.frame(a=tt3[[1]][1, 3], b=tt3[[1]][2, 3], c=tt3[[1]][1, 4], d=tt3[[1]][1, 5]))
tt = rbind(tt, data.frame(a=tt4[[1]][1, 3], b=tt4[[1]][2, 3], c=tt4[[1]][1, 4], d=tt4[[1]][1, 5])) 
tt = rbind(tt, data.frame(a=tt5[[1]][1, 3], b=tt5[[1]][2, 3], c=tt5[[1]][1, 4], d=tt5[[1]][1, 5]))
tt = rbind(tt, data.frame(a=tt6[[1]][1, 3], b=tt6[[1]][2, 3], c=tt6[[1]][1, 4], d=tt6[[1]][1, 5]))
tt = rbind(tt, data.frame(a=tt7[[1]][1, 3], b=tt7[[1]][2, 3], c=tt7[[1]][1, 4], d=tt7[[1]][1, 5]))


colnames(tt) = c("Sum of squares (Age group)", "Sum of squares (residuals)", "F", "p")
rownames(tt) = c("3C.2A", "A1", "A1-2", "A1-3", "A2", "A2-2", "A3")
write.csv(tt, "tt.csv")

library(xtable)
xtable(tt, type="latex")

print ( xtable(tt, type="latex"), file = "../stat_test/ANOVA.tex" )

