library(ggplot2)
library(dplyr)
library(reshape2)
library(lsa)

#read titer data
source("../../../Data/script/load_data/load_data.R")

######################################################################


cosine_theme = theme_bw() + theme(
  axis.text.x = element_text(size=8),
  axis.text.y = element_text(size=8),
  axis.title = element_text(size=9),
  strip.text = element_text(size=8),
  legend.title=element_text(size=8),
  legend.text = element_text(size=8),
  legend.key.size = unit(0.35, "cm"),
  panel.grid = element_blank(),
  strip.background = element_rect(fill = "white", colour = "white"))

################################################################

# First, estimate expected Cosine similarity of 8-dimensional vectors
# with elements uniformly drawn from discrete set of existing dilutions
null_cosine_sim <- replicate(100000,
                             {
                               x <- 2^sample(0:9, size = 8, replace = T)*10
                               y <- 2^sample(0:9, size = 8, replace = T)*10
                               
                               cosine(log(x, base = 2),
                                      log(y, base = 2))
                             }
) %>% mean()
# ==== MV: Generalizing some previous code by KK
plot_similarity_vs_age_single_pairing <- function(ages, similarities){
  p_one = ggplot(df, aes(x=ages, y=similarities)) +
    geom_point() +
    geom_smooth(method = "lm", col="brown1", alpha=0) +
    cosine_theme +
    xlab("Age (years)")
}

plot_similarity_vs_age_replicate_pairings <- function(predicts){
  
  # mean of splines from randomized data
  lm.MEAN = data.frame(xseq, rowMeans(predicts[,c(2:n_replicate_pairings +1)]))
  colnames(lm.MEAN) = c("x", "y")
  
  predicts_sorted = t(apply(predicts,1,sort))
  
  # CI of splines from randomized data
  lm.CI = data.frame(xseq, predicts_sorted[,round(n_replicate_pairings * 0.025)+1],
                     predicts_sorted[,round(n_replicate_pairings * 0.975)+1])
  
  colnames(lm.CI) = c("x", "ci1", "ci2")
  
  ggplot() +
    geom_ribbon(data = lm.CI, aes(ymin = ci1, ymax = ci2, x=xseq), fill = "grey", alpha=0.5) +
    geom_smooth(aes_auto(lm.MEAN), data=lm.MEAN, stat="identity", col="black") +
    cosine_theme + 
    xlab("Age (years)")
}
# ====

# Remove individuals with all-undetectable titers. Remove individuals with NA

sera = sera[rowSums(is.na(sera))==0, ]
sera = sera[rowSums(sera[ ,c(3:10)])!=0, ]

# Do this n_replicate_pairing times:
# For a randomly drawn individual, randomly match another individual with +/- 3 age difference
# Do each drawing without replacement
# Do this until no more pair can be matched 

xseq = seq(from=1, to=90, length=90)
cosine_vs_age_slopes = c()
cosine_predicts = data.frame(xseq)
cor_predicts = data.frame(xseq)

# For this many replicate pairings
n_replicate_pairings <- 1000

# When people have values at LOD, perform this many imputations
n_replicate_imputations <- 1000


for (r in 1:n_replicate_pairings) {
  cosine_similarities = c()
  spearman_correlations <- c()
  ages = c()
  IDs = sera$Sample_ID
  
  # Pair individuals of similar age, compute cosine similarity, Spearman correlation
  for (i in 1:(length(sera$Sample_ID))/2 ){
    
    if (length(IDs) == 0) {
      break
    }
    indiv1 = sample(IDs, 1)
    
    IDs = IDs[ !IDs == indiv1]
    indiv1_age = sera[sera$Sample_ID == indiv1,]$Age  
    age_window = sera[sera$Age >= indiv1_age - 3 & sera$Age <= indiv1_age + 3, ]$Sample_ID
    age_window = intersect(age_window, IDs)
    
    try(
      {
        if (length(age_window) == 0) {
          continue
        }
        indiv2 = sample(age_window, 1)
        IDs = IDs[ !IDs == indiv2]
        indiv2_age = sera[sera$Sample_ID == indiv2,]$Age  
        
        titer1 = as.numeric(sera[sera$Sample_ID == indiv1 , c(3:10)])
      
        titer2 = as.numeric(sera[sera$Sample_ID == indiv2 , c(3:10)])
        
        # Re-calculate log2 titers without dividing by 10
        titer1 <- log(2^titer1*10, base = 2)
        titer2 <- log(2^titer2*10, base = 2)
        
        #If any values are at the limit of detection (i.e., originally coded as 10)
        if(any(titer1 == log(10, base = 2)) | any(titer2 == log(10, base = 2))){
          # For any titers coded as 10 (LOD), impute a random uniform value between 1-20
          # Do this n_replicate_imputations times for this pair of people, compute mean cosine similarity

          replicate(n = n_replicate_imputations,
                    {
                    # For any titers coded as 10 (LOD), impute a random uniform value between 1-20
                    titer1[titer1 == log(10, base = 2)] <- log(runif(n = length(titer1[titer1 == log(10, base = 2)]),
                                                                     1, 20),
                                                               base = 2)

                    titer2[titer2 == log(10, base = 2)] <- log(runif(n = length(titer2[titer2 == log(10, base = 2)]),
                                                                     1, 20),
                                                               base = 2)
                    
                     tibble(cosine_sim = cosine(titer1, titer2)[1] - null_cosine_sim,
                            # Note that Spearman correlation will be NA if an individual has the same titer to all viruses
                            spearman_cor = cor.test(titer1, titer2, method = 'spearman')$estimate)
                    }, simplify = F

          ) %>%
            bind_rows() %>% 
            summarise(across(everything(), mean)) -> mean_stats_across_imputations
          
          cosine_similarity <- mean_stats_across_imputations$cosine_sim
          spearman_cor <- mean_stats_across_imputations$spearman_cor
            
        }else{
          cosine_similarity <- cosine(titer1, titer2) - null_cosine_sim
          spearman_cor <- cor.test(titer1, titer2, method = 'spearman')$estimate
        }
        
        age = (indiv1_age + indiv2_age)/2
        cosine_similarities = c(cosine_similarities, cosine_similarity)
        spearman_correlations <- c(spearman_correlations, spearman_cor)
        ages = c(ages, age)
      },silent=T
    )
    
  }

  df = data.frame(cosine_similarities, spearman_correlations, ages)
  
  cosine_vs_age_model = lm(cosine_similarities ~ ages, data=df)
  cor_vs_age_model = lm(spearman_correlations ~ ages, data=df)
  
  cosine_predicts <- cbind(cosine_predicts,
                           predict(cosine_vs_age_model, newdata = data.frame(ages = xseq),
                                   se=TRUE)$fit)
  
  cosine_vs_age_slopes = c(cosine_vs_age_slopes, summary(cosine_vs_age_model)$coefficients[2,1])
  
  cor_predicts <- cbind(cor_predicts,
                        predict(cor_vs_age_model, newdata = data.frame(ages = xseq),
                                se=TRUE)$fit)
  
  print (r)
}

cosine_vs_age_slopes = sort(cosine_vs_age_slopes)
write.csv(cosine_vs_age_slopes, "../result/cosine_vs_age_slopes_lm_no_replacement_1000_age_window_3.csv")
quantile(cosine_vs_age_slopes, c(0.025, 0.975))

saveRDS(cosine_predicts, "../result/cosine_predicts_lm_no_replacement_1000_age_window_3.rds")


#One example of cosine similarity vs age (single pairing)
saveRDS(df, "../result/cosine_similarity_lm_one_example_age_window_3.rds")

plot_similarity_vs_age_single_pairing(ages = ages, similarities = cosine_similarities) +
  ylab("Normalized cosine similarity")
ggsave("../figure/cosine_similarity_lm_one_example_age_window_3.png", height=2.7, width=4.5)
ggsave("../figure/cosine_similarity_lm_one_example_age_window_3.pdf", height=2.7, width=4.5)

write.csv(data.frame(unlist(summary(cosine_vs_age_model)$coefficients)),
          "../stat_test/cosine_similarity_lm_one_example_age_window_3.csv", row.names=T)

#One example of cosine similarity vs age (single pairing)
plot_similarity_vs_age_single_pairing(ages = ages, similarities = spearman_correlations) +
  ylab("Spearman correlation")
ggsave("../figure/spearman_cor_lm_one_example_age_window_3.png", height=2.7, width=4.5)
ggsave("../figure/spearman_cor_lm_one_example_age_window_3.pdf", height=2.7, width=4.5)

write.csv(data.frame(unlist(summary(cor_vs_age_model)$coefficients)),
          "../stat_test/spearman_cor_lm_one_example_age_window_3.csv", row.names=T)

# Cosine similarity vs age: mean relationship across replicate pairings

plot_similarity_vs_age_replicate_pairings(predicts = cosine_predicts) +
  ylab('Normalized cosine similarity')
ggsave("../figure/cosine_similarity_lm_1000_age_window_3.png", height=2.7, width=4.5)
ggsave("../figure/cosine_similarity_lm_1000_age_window_3.pdf", height=2.7, width=4.5)

plot_similarity_vs_age_replicate_pairings(predicts = cor_predicts) +
  ylab('Spearman correlation')
ggsave("../figure/spearman_cor_lm_1000_age_window_3.png", height=2.7, width=4.5)
ggsave("../figure/spearman_cor_lm_1000_age_window_3.pdf", height=2.7, width=4.5)



