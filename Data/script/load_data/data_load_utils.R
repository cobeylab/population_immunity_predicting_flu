sample_or_load = function(filename, model, data, chains = 4, cores = 4) {
  #Ed's code
  if(file.exists(filename)) {
    readRDS(filename)
  }
  else {
    fit <- stan(model, data=data, chains = chains, cores = cores)
    saveRDS(fit, filename)
    fit
  }
}

get_idx = function(var_group){
  idx = c() #index of columns with variations
  for(i in 1:ncol(var_group)){
    if (var(var_group[,i], na.rm=T)!=0){
      idx = c(idx, i)
    }
  }
  
  id = c() #index of columns not duplicated
  if(!is.null(idx)){
    
    uniq = !duplicated(t(var_group))
    for (i in 1:length(uniq)){
      if(uniq[i] & (i %in% idx)) {
        id = c(id, i)
      }
    } 
    dim(id) = length(id)
    
  }else{
    id = numeric()
  }
  
  id
}

get_correlated_id = function(shared_mutImp_3a){
  cor3a = cor(shared_mutImp_3a)
  cor3a[upper.tri(cor3a)] = NA
  
  i_rm3a = c()
  for(i in 1:ncol(cor3a)) {
    for(j in 1:nrow(cor3a)) {
      if (is.na(cor3a[j,i]) == T ){
        next
      }
      if( (cor3a[j,i] > 0.9) & (i != j) ){
        i_rm3a = c(i_rm3a, j)
      }
    }
  }
  i_rm3a = i_rm3a[!duplicated(i_rm3a)]
  
  i_rm3a
}

mutation_times_impriting_prob = function(mut, imp_prob, prefix_imp_prob, prefix_new_var){
  p_imp_prob_mut = data.frame(matrix(nrow=nrow(mut), ncol=0))
  for (i in 1:ncol(mut)){
    site = strsplit(colnames(mut)[i], "_")[[1]][2]
    id = grep(paste0("\\b",prefix_imp_prob,site,"\\b"), colnames(imp_prob))
    p_imp_prob_mut = data.frame(p_imp_prob_mut, imp_prob[,id]*mut[,i])
  }
  colnames(p_imp_prob_mut) = paste0(prefix_new_var, "_", colnames(mut)) 
  
  p_imp_prob_mut
}

add_epitope_to_mutImp_colnames = function(mut, imp_prob, prefix_imp_prob, prefix_new_var){
# Mapping mutImp to mutation. For example, m142 will be mapped to A_142_K and A_142_G
# Epitope site where each mutation belong will be also marked to colnames(mutImp)
  
  p_imp_prob_mut = data.frame(matrix(nrow=nrow(mut), ncol=0))

  for (i in 1:ncol(mut)){
    site = strsplit(colnames(mut)[i], "_")[[1]][2]
    id = grep(paste0("\\b",prefix_imp_prob,site,"\\b"), colnames(imp_prob))
    p_imp_prob_mut = data.frame(p_imp_prob_mut, imp_prob[,id])
  }
  colnames(p_imp_prob_mut) = paste0(prefix_new_var, "_", colnames(mut)) 
  
  p_imp_prob_mut
}

