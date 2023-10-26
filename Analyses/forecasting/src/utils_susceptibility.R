

get_frac_under_threshold = function(x, threshold) {
  x = as.vector(x)
  x = x[!is.na(x)]
  return ( sum(x < threshold)/length(x) )
}

susceptibility_frac_unprotected = function(sera, group, viruses, threshold) {
  
  get_susceptibility = function(x) {
    susc = x*frac$popfrac
    susc = sum(susc, na.rm=T)
    return (susc)
  }
  
  frac = sera %>%
    group_by_(.dots = group) %>%
    summarize_at(viruses, get_frac_under_threshold, threshold ) 
  
  frac = data.frame(frac)
  frac = merge(frac, pop, by= group[1] )
  frac$popfrac = frac$n / sum(frac$n)
 
  colnames(frac) = gsub("X", "", colnames(frac))
  
  susceptibility = sapply(frac[,viruses], get_susceptibility)
  
  return (susceptibility)  
}

#################################################################
get_susc_gmt = function(x){
  x = as.vector(x)
  x = x[!is.na(x)]
  susc_gmt = (mean(x))
  return (susc_gmt)
}

susceptibility_gmt = function(sera, group, viruses, threshold) {
  
  get_susceptibility = function(x) {
    susc = gmt$popfrac*x
    susc = -1 * sum(susc, na.rm=T)
    return (susc)
  }
  
  gmt = sera %>%
    group_by_(.dots = group) %>%
    summarize_at(viruses, get_susc_gmt) 
  
  gmt = data.frame(gmt)
  gmt = merge(gmt, pop, by= group[1] )
  gmt$popfrac = gmt$n / sum(gmt$n)
  
  colnames(gmt) = gsub("X", "", colnames(gmt))

  susceptibility = sapply(gmt[,viruses], get_susceptibility)
  
  return (susceptibility)  
}


######################################################################################


get_95CI1 = function (x) {
  x = sort(x)
  
  ci1 = x[ as.integer( 0.025*length(x) ) ]
  
  print (ci1)
}
get_95CI2 = function (x) {
  x = sort(x)
  ci2 = x[ as.integer( 0.975*length(x) ) ]
  
  print (ci2)
}

susceptibility_frac_unprotected_boot = function(sera, group, viruses, threshold) {
  
  get_susceptibility = function(x) {
    susc = x*frac$popfrac
    susc = sum(susc, na.rm=T)
    return (susc)
  }
  
  frac = sera %>%
    group_by_(.dots = group) %>%
    do( sample_n(., size=length(.), replace=T ) ) %>%
    summarize_at(viruses, get_frac_under_threshold, threshold ) 
  
  frac = data.frame(frac)
  frac = merge(frac, pop, by= group[1] )
  frac$popfrac = frac$n / sum(frac$n)

  colnames(frac) = gsub("X", "", colnames(frac))
  
  susceptibility = lapply(frac[,viruses], get_susceptibility)
  
  return (susceptibility)  
}

susceptibility_gmt_boot = function(sera, group, viruses, threshold) {
  
  get_susceptibility = function(x) {
    susc = gmt$popfrac*x
    susc = -1 * sum(susc, na.rm=T)
    return (susc)
  }
  
  gmt = sera %>%
    group_by_(.dots = group) %>%
    do( sample_n(., size=length(.), replace=T ) ) %>%
    summarize_at(viruses, get_susc_gmt) 
  
  gmt = data.frame(gmt)
  gmt = merge(gmt, pop, by= group[1] )
  gmt$popfrac = gmt$n / sum(gmt$n)
  
  colnames(gmt) = gsub("X", "", colnames(gmt))
  
  susceptibility = lapply(gmt[,viruses], get_susceptibility)
  
  return (susceptibility)  
}
