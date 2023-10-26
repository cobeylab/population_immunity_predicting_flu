library(reshape2)
library(ggplot2)
library(vegan)
library(gdata)
#library(stargazer)
#library(lfe)


###################################################################
assign_cluster = function(sera8, cl_fit, fig_dir, analysis){
  
  #png(filename= paste0(fig_dir, "Clustering_Kmean_", analysis, ".png"))
  #plot(cl_fit, sortg = TRUE, grpmts.plot = TRUE)
  #dev.off()
  
  best_cl_membership = get_best_cluster(cl_fit)
  
  cl1 = names(sort(summary(best_cl_membership$cluster), decreasing=T))
  cl2 = paste("Group", seq(1, length(unique(best_cl_membership$cluster))) )
  for(i in 1:length(cl2)) {
    #best_cl_membership$cluster = gsub( paste0("^", cl1[i], "$"), cl2[i], best_cl_membership$cluster)
    best_cl_membership$cluster = gsub(as.character(i), paste0("Group ",as.character(i)),  best_cl_membership$cluster)
    
  }
 
  sera8 = merge(sera8, best_cl_membership, by.x="Sample_ID", by.y="cl_sample_id")
  sera8$cluster = as.factor(sera8$cluster)
  sera8$cluster =relevel(sera8$cluster, "Group 1")
  
  m_sera8 = melt(sera8, id.vars= c("Sample_ID", "Age", "cluster"))
  colnames(m_sera8)[c(4,5)] = c("Virus", "Titer")
  m_sera8$Virus = factor(m_sera8$Virus, 
                         levels = c("3C3.A", "3C2.A", "N171K", "N121K_N171K", "N121K_T135K_N171K",
                                    "T131K_R142K", "T131K_R142K_R261Q", "N121K_S144K"))
  
  #p = ggplot(m_sera8, aes(x=Age, y=Titer)) +
  #  geom_jitter(aes(col=Virus), width=0.2, height=0.2, alpha=0.2) +
  #  geom_smooth(aes(col=Virus), alpha=0, size=1.5) +
  #  facet_wrap(~cluster) +
  #  scale_color_discrete(labels = plot_x_label, name="") 
  
  #print (p)
  #ggsave(paste0(fig_dir, "Clustering_titer_", analysis,  ".png"), width=5, height=3)
  
  sera8
}


age_analysis = function(sera8, fig_dir, analysis){
  m_sera8 = melt(sera8, id.vars= c("Sample_ID", "Age", "cluster"))
  colnames(m_sera8)[c(4,5)] = c("Virus", "Titer")
  m_sera8$Virus = factor(m_sera8$Virus, 
                         levels = c("3C3.A", "3C2.A", "N171K", "N121K_N171K", "N121K_T135K_N171K",
                                    "T131K_R142K", "T131K_R142K_R261Q", "N121K_S144K"))
  
  
  
  sm = summary(lm(Age ~ cluster, data=sera8))
  coef = sm$coefficients[2, 1]
  p = sm$coefficients[2, 4]
  
  p = ggplot(m_sera8, aes(y=Age, x=cluster)) +
    geom_jitter(alpha=0.5, width=0.2, height=0.2) +
    geom_boxplot(alpha=0)
  print (p)
  ggsave(paste0(fig_dir, "Clustering_age_mean_", analysis, ".png"))
  
  write.csv(sm$coefficients, paste0("../table/Clustering_age_mean_", analysis, ".csv"))
  
  #compare distribution
  
  c1 = sera8[sera8$cluster == "Group 1", ]$Age
  c2 = sera8[sera8$cluster == "Group 2", ]$Age
  
  kt = ks.test(c1, c2)
  p = ggplot(m_sera8) +
    geom_histogram(aes(x=Age)) +
    facet_wrap(~cluster)
  print (p)
  ggsave(paste0(fig_dir, "Age_dist_", analysis, ".png"))
  
  write.csv(kt$p.value, paste0("../table/Clustering_age_dist_", analysis, ".csv"))
}

###################################################################
#linear regression for individual subject

fit_each_individual = function(vars_sorted, model, coef_names, printFit){
  subjects = as.character(unique(vars_sorted$Sample_ID))
  
  #fit model to each subject and save coefficient in z
  z = c()
  sid = c()
  for(i in 1:length(subjects)){
    subject = subset(vars_sorted, Sample_ID == subjects[i])
    sm = summary(glm(model, data=subject))
    params = sm$coefficients[,"Estimate"]
    params = as.vector(params)
    if(printFit == 1){
      print(sm)
    }
    if(all(params==0)){
      next
    }
    z = rbind(z, params)
    sid = c(sid, subjects[i])
  }
  
  #make z usable for clustering
  z = as.data.frame(z)
  rownames(z) = sid
  colnames(z) = coef_names
  
  return (z)
}
###################################################################
#K means clustering

Calinski_clustering = function(z, itc_for_clustering, centering_and_scaling){
  if(itc_for_clustering == 0){
    z_for_clustering = z[,2:ncol(z)]
  }else{
    z_for_clustering = z
  }
  
  if(centering_and_scaling == 1){
    fit = cascadeKM(scale(z_for_clustering,center=TRUE,scale=TRUE), 1, 7,ite=1500)
  }else{
    fit = cascadeKM(z_for_clustering, 1, 7,ite=1500)
  }
  
  return(fit)
}

get_best_cluster = function(fit){
  best_num_cl =as.numeric(which.max(fit$results[2,]))
  best_cl_name = paste(best_num_cl, "groups")
  
  cl_sample_id = rownames(fit$partition)
  
  cluster = fit$partition[,best_cl_name]
  best_cl_membership = as.data.frame(cbind(cl_sample_id, cluster)) 
  
  return(best_cl_membership)
}

get_4_cluster = function(fit){
  best_num_cl = 4
  best_cl_name = paste(best_num_cl, "groups")
  
  cl_sample_id = rownames(fit$partition)
  
  cluster = fit$partition[,best_cl_name]
  best_cl_membership = as.data.frame(cbind(cl_sample_id, cluster)) 
  
  return(best_cl_membership)
}

###################################################################
#fit by group


fit_each_group = function(vars_sorted, best_cl_membership, model_for_fitting_group){
  group_names = sort(unique(best_cl_membership$cluster))
  results = list()
  
  for (i in 1:length(group_names)){
    ids_this_group = best_cl_membership[best_cl_membership$cluster == group_names[i],]$cl_sample_id
    
    #fit all individuals in this group together again
    vars_this_group = vars_sorted[vars_sorted$Sample_ID %in% ids_this_group,]
    sm_this_group = summary(glm(model_for_fitting_group, data=vars_this_group))
    
    results[[i]] = sm_this_group
  }
  
  return (results)
}

fit_each_group_return_raw = function(vars_sorted, best_cl_membership, model_for_fitting_group){
  group_names = sort(unique(best_cl_membership$cluster))
  results = list()
  
  for (i in 1:length(group_names)){
    ids_this_group = best_cl_membership[best_cl_membership$cluster == group_names[i],]$cl_sample_id
    
    #fit all individuals in this group together again
    vars_this_group = vars_sorted[vars_sorted$Sample_ID %in% ids_this_group,]
    fit = glm(model_for_fitting_group, data=vars_this_group)
    
    results[[i]] = fit
  }
  
  return (results)
}

get_ages_each_group = function(best_cl_membership){
  group_names = sort(unique(best_cl_membership$cluster))
  
  ages_each_group = c()
  for (i in 1:length(group_names)){
    ids_this_group = best_cl_membership[best_cl_membership$cluster == group_names[i],]$cl_sample_id
    
    #for testing mean group age: age of each group
    ages = sera[sera$Sample_ID %in% ids_this_group,]$Age
    ages_each_group = rbind(ages_each_group, 
                            data.frame(rep(group_names[i], length(ages)), ages))
  }
  colnames(ages_each_group) = c("group", "ages")

  return(ages_each_group)
}


print_fit_by_group = function(fit_by_group){
  for (i in 1:length(fit_by_group)){
    cat("\n")
    print (paste("group", i))
    sm_this_group = fit_by_group[[i]]
    result = (sm_this_group$coefficients[1:7,])
    pvals = as.vector(result[,4])
    stars = get_stars(pvals)
    result = as.data.frame(result)
    result$siginficance = stars
    print (paste("aic:", sm_this_group$aic))
    print (result)
  }
}


###################################################################
#other utilities

save_coefficients = function(fit, mean_age){
  coefs = (fit$coefficients)
  pvals = (get_stars(fit$pval))
  muts = rownames(fit$coefficients)
  ages = rep(mean_age, nrow(coefs))
  result = data.frame(muts, ages, coefs, pvals)
  return(result)
}


save_coefficients_glm = function(fit, mean_age){
  sm = summary(fit)
  muts = linear_model_coef[2:length(linear_model_coef)]
  mut_coefs = sm$coefficients[muts,]
  coefs = mut_coefs[,1]
  pvals = mut_coefs[,4]
  stars = get_stars(pvals)
  ages = rep(mean_age, length(coefs))
  result = data.frame(muts, ages, coefs, pvals, stars)
  return(result)
}


get_stars = function(pval){
  stars = c()
  pvals = as.numeric(as.character((pval)))
  
  for (i in 1:length(pval)){
    if (pval[i] < 0.0001){
      star = "***"
    }else if (pval[i] < 0.001) {
      star = "**"
    }else if (pval[i] < 0.05) {
      star = "*"
    }else if (pval[i] < 0.1) {
      star = "."
    }else {
      star = ""
    }
    stars = c(stars, star)
  }
  return (stars)
}
