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


impute_continuous_titers <- function(log2_titer){
  original_titer <- 2^log2_titer
  
  # For any titers coded as 1:10, recode as 1:1
  original_titer[abs(original_titer) - 10 < 1e-5] <- 1
  
  # Sample continuous value between original titer and next dilution
  # (Next dilution is 20 if original titer coded as 1:1, *2 otherwise)
  next_dilutions <- original_titer*ifelse(original_titer ==1, 20, 2)
  
  continuous_titers <- runif(length(original_titer),
                             min = original_titer, max = next_dilutions)
  
  return(log(continuous_titers, base = 2))
}