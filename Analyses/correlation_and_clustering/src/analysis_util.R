
read_predicted_values = function(fname, df, value="prdd") {
  if (  grepl("prdd", fname) ){

    prdd = readRDS(fname)
    return (prdd)
    
  }else if (  grepl("pred", fname) ){
    
    pred = readRDS(fname)
    return (pred)
    
  }else if (grepl("fit", fname) ){
  
    fit = readRDS(fname)
    df$pred = summary(fit, pars = c("H_reg"))$summary[,"mean"]
    
    if(value == "pred") {
      return (df$pred)
    }
    
    sigma = summary(fit, pars = c("sigma"))$summary[,"mean"]
    sigma_na = summary(fit, pars = c("sigma_na"))$summary[,"mean"]
    
    df$prdd = ifelse(df$is_NA == 0, 
                     floor(df$pred + rnorm(nrow(df), 0, sigma)),
                     floor(df$pred + rnorm(nrow(df), 0, sigma_na)))
    
    df$prdd = ifelse(df$prdd < 0, 0, df$prdd)
    return (df$prdd)
  }
  
}
