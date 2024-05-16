library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(reshape2)
library(lsa)
library(stringr)

#read titer data
source("../../../Data/script/load_data/load_data.R")
source("analysis_util.R")

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
                               # Randomly sample titer from the observed dilutions
                               x <- 2^sample(0:9, size = 8, replace = T)*10
                               y <- 2^sample(0:9, size = 8, replace = T)*10
                               
                               # Randomly impute continue titer values
                               imputed_log2_x <- impute_continuous_titers(log(x, base= 2))
                               imputed_log2_y <- impute_continuous_titers(log(y, base =2))
                               
                               cosine(imputed_log2_x,
                                      imputed_log2_y)
                             }
) %>% mean()
# ==== MV: Generalizing some previous code by KK
plot_similarity_vs_age_single_pairing <- function(ages, similarities){
  p_one = ggplot(single_pairing_similarity_df, aes(x=ages, y=similarities)) +
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

plot_pairs <- function(single_pairing_similarity_df, sera, pairs){
  
  long_form_serum_data <- as_tibble(sera) %>%
    pivot_longer(cols = !any_of(c('Sample_ID', 'Age')),
                 names_to = 'clade', values_to = 'log2_of_titer_by_10') %>%
    # Compute log2 titer without the division by 10
    mutate(log2_titer = log(2^log2_of_titer_by_10*10, base = 2))
  
  unique_values <- unique(long_form_serum_data$log2_titer)

  wide_format_pairs <-
    lapply(as.list(pairs),
         FUN = function(pair){
           
           pair_label <- single_pairing_similarity_df %>%
             filter(pairs == pair) %>%
             mutate(label = paste0(pairs, "\n",
                                   'Mean age = ', ages, "\n",
                                   'Cosine sim. = ', round(cosine_similarities,3), "\n",
                                   'Spearman cor = ', round(spearman_correlations, 3))) %>%
             pull(label)
           
           pair <- str_split(pair, ';')[[1]]
           
           long_form_serum_data %>%
             filter(Sample_ID %in% pair) %>%
             select(clade, Sample_ID, log2_titer) %>%
             pivot_wider(names_from = Sample_ID, values_from = log2_titer) %>%
             rename_with(.cols = all_of(pair), .fn = function(x){paste0('person_',which(pair == x))}) %>%
             mutate(pair = paste0(pair, collapse = ';'),
                    pair_label = pair_label) %>%
             select(pair, everything()) %>%
             # For each point in scatterplot, make size a functioin of how many values are at that point
             mutate(position = paste(person_1, person_2, sep = ',')) %>%
             group_by(position) %>%
             mutate(point_color = n()) %>%
             ungroup() %>%
             select(-position)
         }
         
         ) %>%
    bind_rows()
  
  scatterplot <- wide_format_pairs %>%
    ggplot(aes(x = person_1, y = person_2)) +
    geom_point(shape = 21, size = 3, aes(fill = factor(point_color))) +
    geom_text(aes(label = point_color, color = point_color < 5), size = 2) +
    scale_x_continuous(breaks = unique_values, limits = c(range(unique_values)),
                       labels = ~ 2^(.x) ) +
    scale_y_continuous(breaks = unique_values, limits = c(range(unique_values)),
                       labels = ~ 2^(.x) ) +
    facet_wrap('pair_label') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = -90, vjust = 0, size = 9),
          axis.text.y = element_text(size = 9),
          strip.text = element_text(size = 9),
          axis.title = element_text(size = 9),
          legend.position = 'none') +
    xlab('Person 1') +
    ylab('Person 2') +
    scale_fill_brewer(name = 'n values') +
    scale_color_manual(values = c('white','black'))
  
  long_format_pairs <- wide_format_pairs %>%
    clade_relabeller() %>%
    select(pair, clade, pair_label, matches('person')) %>%
    pivot_longer(cols = matches('person'), names_to = 'person', values_to = "log2_titer")
    
  
  
  landscape_plot <- long_format_pairs %>%
    ggplot(aes(x = clade, y = log2_titer, color = person, group = person)) + 
    geom_point() +
    geom_line() +
    theme_cowplot() +
    facet_wrap('pair_label') +
    scale_y_continuous(breaks = unique_values, limits = c(range(unique_values)),
                       labels = ~ 2^(.x) )  +
    scale_x_discrete(labels = function(x){str_remove(x,'3C\\.')}) +
    theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5, size = 10),
          legend.position = 'none') +
    xlab('Clade') +
    ylab('Titer')
  
  return(list(scatterplot = scatterplot, landscape_plot = landscape_plot))
  
  
}

clade_relabeller <- function(data){
  data %>%
    mutate(clade = case_when(
      clade == '3C3.A' ~ '3C.3A',
      clade == '3C2.A' ~ '3C.2A',
      clade == 'N171K' ~ '3C.2A1-1',
      clade == 'N121K_N171K'~ '3C.2A1-2',
      clade == 'N121K_T135K_N171K' ~ '3C.2A1-3',
      clade == 'T131K_R142K'~ '3C.2A2-1',
      clade == 'T131K_R142K_R261Q' ~ '3C.2A2-2',
      clade == 'N121K_S144K' ~ '3C.2A3'
    )) %>%
    mutate(clade = factor(
      clade,
      levels = c('3C.3A', '3C.2A','3C.2A1-1', '3C.2A1-2','3C.2A1-3','3C.2A2-1',
                 '3C.2A2-2','3C.2A3')))
}

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
n_replicate_imputations <- 100


for (r in 1:n_replicate_pairings) {
  cosine_similarities = c()
  spearman_correlations <- c()
  ages = c()
  pairs <- c()
  n_values_below_LOD <- c() # Number of values below LOD across both people in the pair 
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
        
        n_values_below_LOD <- c(n_values_below_LOD,
                          sum(round(2^titer1) == 10) +
                            sum(round(2^titer2) == 10))
        
        # Randomly impute continuous titer values
        replicate(n = n_replicate_imputations,
                  {
                    continuous_titers1 <- impute_continuous_titers(titer1)
                    continuous_titers2 <- impute_continuous_titers(titer2)
                    
                    tibble(cosine_sim = (cosine(continuous_titers1,  continuous_titers2)[1] - null_cosine_sim) / (1 - null_cosine_sim),
                           # Note that Spearman correlation will be NA if an individual has the same titer to all viruses
                           spearman_cor = cor.test(continuous_titers1,  continuous_titers2, method = 'spearman')$estimate)
                  }, simplify = F
                  
        ) %>%
          bind_rows() %>% 
          summarise(across(everything(), mean)) -> mean_stats_across_imputations
        
      
        cosine_similarity <- mean_stats_across_imputations$cosine_sim
        spearman_cor <- mean_stats_across_imputations$spearman_cor
            
        
        age = (indiv1_age + indiv2_age)/2
        cosine_similarities = c(cosine_similarities, cosine_similarity)
        spearman_correlations <- c(spearman_correlations, spearman_cor)
        ages = c(ages, age)
        pairs <- c(pairs, paste0(indiv1, ';', indiv2))
      },silent=T
    )
    
  }

  single_pairing_similarity_df = tibble(cosine_similarities, spearman_correlations, ages, pairs, n_values_below_LOD)
  
  cosine_vs_age_model = lm(cosine_similarities ~ ages, data=single_pairing_similarity_df)
  cor_vs_age_model = lm(spearman_correlations ~ ages, data=single_pairing_similarity_df)
  
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
saveRDS(single_pairing_similarity_df, "../result/cosine_similarity_lm_one_example_age_window_3.rds")

plot_similarity_vs_age_single_pairing(ages = ages, similarities = cosine_similarities) +
  ylab("Normalized cosine similarity")
ggsave("../figure/cosine_similarity_lm_one_example_age_window_3.png", height=2.7, width=4.5)
ggsave("../figure/cosine_similarity_lm_one_example_age_window_3.pdf", height=2.7, width=4.5)

write.csv(data.frame(unlist(summary(cosine_vs_age_model)$coefficients)),
          "../stat_test/cosine_similarity_lm_one_example_age_window_3.csv", row.names=T)

#One example of cosine similarity vs age (single pairing)
plot_similarity_vs_age_single_pairing(ages = ages, similarities = spearman_correlations) +
  ylab("Spearman correlation") + ylim(-1,1)
ggsave("../figure/spearman_cor_lm_one_example_age_window_3.png", height=2.7, width=4.5)
ggsave("../figure/spearman_cor_lm_one_example_age_window_3.pdf", height=2.7, width=4.5)

write.csv(data.frame(unlist(summary(cor_vs_age_model)$coefficients)),
          "../stat_test/spearman_cor_lm_one_example_age_window_3.csv", row.names=T)

# Cosine similarity vs age: mean relationship across replicate pairings

plot_similarity_vs_age_replicate_pairings(predicts = cosine_predicts) +
  ylab('Normalized cosine similarity') +
  ylim(0,1)
ggsave("../figure/cosine_similarity_lm_1000_age_window_3.png", height=2.7, width=4.5)
ggsave("../figure/cosine_similarity_lm_1000_age_window_3.pdf", height=2.7, width=4.5)

plot_similarity_vs_age_replicate_pairings(predicts = cor_predicts) +
  ylab('Spearman correlation') +
  ylim(-1,1)
ggsave("../figure/spearman_cor_lm_1000_age_window_3.png", height=2.7, width=4.5)
ggsave("../figure/spearman_cor_lm_1000_age_window_3.pdf", height=2.7, width=4.5)

# For a single random pairing (the last in the loop), make scatterplots for some pairs

# 12 pairs with the highest cosine similarity
top_12_cosine_sim <- plot_pairs(single_pairing_similarity_df, sera = sera,
           pairs = single_pairing_similarity_df %>%
             arrange(desc(cosine_similarities)) %>%
             slice(1:12) %>% pull(pairs)
)

top_12_cosine_sim$scatterplot + ggtitle("12 pairs with the highest cosine similarity")
top_12_cosine_sim$landscape_plot + ggtitle("12 pairs with the highest cosine similarity")

# 12 pairs with the lowest cosine similarity
bottom_12_cosine_sim <- plot_pairs(single_pairing_similarity_df, sera = sera,
           pairs = single_pairing_similarity_df %>%
             arrange(cosine_similarities) %>%
             slice(1:12) %>% pull(pairs)
)

bottom_12_cosine_sim$scatterplot + ggtitle("12 pairs with the lowest cosine-similarity")
bottom_12_cosine_sim$landscape_plot + ggtitle("12 pairs with the lowest cosine-similarity")

# 12 pairs with the highest Spearman cor
top_12_spearman_cor <- plot_pairs(single_pairing_similarity_df, sera = sera,
                                pairs = single_pairing_similarity_df %>%
                                  arrange(desc(spearman_correlations)) %>%
                                  slice(1:12) %>% pull(pairs)
)

top_12_spearman_cor$scatterplot +  ggtitle("12 pairs with the highest Spearman correlation") 
top_12_spearman_cor$landscape_plot +  ggtitle("12 pairs with the highest Spearman correlation") 

# 12 pairs with the lowest Spearman cor

bottom_12_spearman_cor <- plot_pairs(single_pairing_similarity_df, sera = sera,
                                  pairs = single_pairing_similarity_df %>%
                                    arrange((spearman_correlations)) %>%
                                    slice(1:12) %>% pull(pairs)
)

bottom_12_spearman_cor$scatterplot +  ggtitle("12 pairs with the lowest Spearman correlation") 
bottom_12_spearman_cor$landscape_plot +  ggtitle("12 pairs with the lowest Spearman correlation") 


# Plotting cosine similarity as a function of n values below LOD
single_pairing_similarity_df %>%
  ggplot(aes(x = n_values_below_LOD, y = cosine_similarities)) + 
  geom_point() +
  theme_cowplot() +
  xlab('Number of values below the LOD across both people') +
  ylab('Cosine similarity\n(normalized, with random imputations)')


# Linear model linking cosine similarity to spearman correlation
cossim_vs_spearman <- lm(cosine_similarities ~ spearman_correlations,
                         data = single_pairing_similarity_df)

single_pairing_similarity_df <- single_pairing_similarity_df %>% 
  mutate(residual = cossim_vs_spearman$residuals,
         rank_abs_residual = rank(abs(residual))) %>%
  mutate(selected = case_when(
    rank_abs_residual <= 12 | (rank_abs_residual >= (max(rank_abs_residual) - 12 + 1)) ~ T,
    T ~ F
  ))


# Plotting cosine similarity vs spearman correlation
single_pairing_similarity_df  %>%
  ggplot(aes(x = spearman_correlations, y = cosine_similarities)) + 
  theme_cowplot() +
  geom_point(aes(color = selected)) +
  xlab('Spearman correlation') +
  ylab('Cosine similarity') +
  geom_smooth(method = 'lm', color = 'black') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('black','red'))


# Plot the pairs with the greatest absolute residuals
plot_pairs(single_pairing_similarity_df, sera = sera,
           pairs = single_pairing_similarity_df %>%
             filter(selected) %>%
             arrange(desc(residual)) %>%
             slice(1:12) %>%
             pull(pairs))

# Plots with the smallest absolute residuals
plot_pairs(single_pairing_similarity_df, sera = sera,
           pairs = single_pairing_similarity_df %>%
             filter(selected) %>%
             arrange(residual) %>%
             slice(1:12) %>%
             pull(pairs))
