library(reshape2)

data_cleaning = function(sera, rmvdup){
  #Remove unnecessary test viruses and columns which are with
  #"ND_90(Some other test viruses)","Micronic" and "birth"
  ex_idx = grep("ND_90", colnames(sera))
  if(length(ex_idx) == 0){
    sera = sera
  }else{
    sera = sera[,-c(ex_idx)]
  }  
  
  ex_idx = grep("Micronic", colnames(sera))
  if(length(ex_idx) == 0){
    sera = sera
  }else{
    sera = sera[,-c(ex_idx)]
  }  
  
  ex_idx = grep("birth", colnames(sera))
  if(length(ex_idx) == 0){
    sera = sera
  }else{
    sera = sera[,-c(ex_idx)]
  }  
  
  ex_idx = grep("Order", colnames(sera))
  if(length(ex_idx) == 0){
    sera = sera
  }else{
    sera = sera[,-c(ex_idx)]
  }  
  
  #remove rows with NA
  sera = na.omit(sera)
  
  #remove duplicated subjects
  if(rmvdup == 1){
    sera = sera[!duplicated(sera$Sample_ID),]
  }
  
  
  return (sera)
}

log_titer = function(sera) {
  #divide by 10 and log2
  tv_idx = grep("ND90", colnames(sera))
  sera[,tv_idx] = log2(sera[,tv_idx]/10)
  
  #add mean titer of each individual
  #sera$mean_titer = rowMeans(sera[,tv_idx])
  
  return (sera)
}

log_titer_varname = function(sera, varname) {
  #divide by 10 and log2
  tv_idx = grep(varname, colnames(sera))
  sera[,tv_idx] = log2(sera[,tv_idx]/10)
  
  #add mean titer of each individual
  #sera$mean_titer = rowMeans(sera[,tv_idx])
  
  return (sera)
}

change_virus_names = function(sera){
  
  colnames(sera) = gsub("ND90_", "", colnames(sera))
  colnames(sera) = gsub("Colo14_", "", colnames(sera))
  colnames(sera) = gsub("FRNT_", "", colnames(sera))  
  
  sera
}

melt_data = function(sera) {
  #melt data
  tv_idx = grep("ND90", colnames(sera))
  tv_names = colnames(sera)[tv_idx]
  m_sera = melt(sera, id.var = c("Sample_ID", "Age"), measure.vars = tv_names)
  
  #change column names
  var_idx = grep("variable", colnames(m_sera))
  val_idx = grep("value", colnames(m_sera))
  colnames(m_sera)[var_idx] = "test_virus"
  colnames(m_sera)[val_idx] = "titer"
  
  return (m_sera)
}

short_tvnames = function(df, tv){
  df[,tv] = gsub("ND90_", "", df[,tv])
  df[,tv] = gsub("Colo14_", "", df[,tv])
  df[,tv] = gsub("WT", "3C.2a", df[,tv])
  df[,tv] = factor(df[,tv], levels=c('3C.3a', '3C.2a', 'N171K', 'N121K_N171K', 'N121K_T135K_N171K', 'N121K_S144K', 'T131K_R142K' ))
  
  return (df[,tv])
}


