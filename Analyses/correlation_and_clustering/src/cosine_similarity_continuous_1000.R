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



# Remove individuals with all-undetectable titers. Remove individuals with NA

sera = sera[rowSums(is.na(sera))==0, ]
sera = sera[rowSums(sera[ ,c(3:10)])!=0, ]

# Do this 1000 times:
# For a randomly drawn individual, randomly match another individual with +/- 3 age difference
# Do each drawing without replacement
# Do this until no more pair can be matched 

xseq = seq(from=1, to=90, length=90)
slops = c()
pvals = c()
predicts = data.frame(xseq)

for (r in 1:1000) {
  cosine_similarities = c()
  ages = c()
  IDs = sera$Sample_ID
  
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
        
        cosine_similarity = cosine(titer1, titer2)
        age = (indiv1_age + indiv2_age)/2
        cosine_similarities = c(cosine_similarities, cosine_similarity)
        ages = c(ages, age)
      },silent=T
    )
    
  }

  df = data.frame(cosine_similarities, ages)
  
  model = lm(cosine_similarities ~ ages, data=df)
  pred = predict(model, newdata = data.frame(ages = xseq), se=TRUE)
  y = pred$fit
  predicts = cbind(predicts, y)
  
  sm_fit = summary(model)
  slop = sm_fit$coefficients[2,1]
  slops = c(slops, slop)
  p = sm_fit$coefficients[2,4]
  pvals = c(pvals, p)
  print (r)
}

slops = sort(slops)
write.csv(slops, "../result/slops_lm_no_replacement_1000_age_window_3.csv")
slops[25]
slops[975]

saveRDS(predicts, "../result/predicts_lm_no_replacement_1000_age_window_3.rds")


#mean of splines from randomized data
lm.MEAN = data.frame (xseq, rowMeans(predicts[,c(2:1001)]))
colnames(lm.MEAN) = c("x", "y")

#CI of splines from randomized data
predicts_sorted = t(apply(predicts,1,sort))
lm.CI = data.frame(xseq, predicts_sorted[,25+1], predicts_sorted[,975+1])
colnames(lm.CI) = c("x", "ci1", "ci2")

#One example of cosine similarity 
saveRDS(df, "../result/cosine_similarity_lm_one_example_age_window_3.rds")

p_one = ggplot(df, aes(x=ages, y=cosine_similarities)) +
  geom_point() +
  geom_smooth(method = "lm", col="brown1", alpha=0) +
  cosine_theme +
  xlab("Age (years)") + ylab("Cosine similarity")
p_one
ggsave("../figure/cosine_similarity_lm_one_example_age_window_3.png", height=2.7, width=4.5)
ggsave("../figure/cosine_similarity_lm_one_example_age_window_3.pdf", height=2.7, width=4.5)

sm_df = data.frame(unlist(sm_fit$coefficients))
write.csv(sm_df, "../stat_test/cosine_similarity_lm_one_example_age_window_3.csv", row.names=T)

#####################################################################33

p_lm_1000 = ggplot() +
  geom_ribbon(data = lm.CI, aes(ymin = ci1, ymax = ci2, x=xseq), fill = "grey", alpha=0.5) +
  geom_smooth(aes_auto(lm.MEAN), data=lm.MEAN, stat="identity", col="black") +
  cosine_theme + 
  ylim(c(0,1)) + xlim(c(1,90)) +
  xlab("Age (years)") + ylab("Cosine similarity") 
p_lm_1000
ggsave("../figure/cosine_similarity_lm_1000_age_window_3.png", height=2.7, width=4.5)
ggsave("../figure/cosine_similarity_lm_1000_age_window_3.pdf", height=2.7, width=4.5)

##########################################################################3
# Now randomly redistribute ages and do the same thing


slops_shuffled = c()
pvals_shuffled = c()
predicts_shuffled = data.frame(xseq)


for (r in 1:1000){
  
  sera_shuffled = transform(sera, Age=sample(sera$Age))
  cosine_similarities_shuffled = c()
  ages_shuffled = c()
  IDs = sera_shuffled$Sample_ID
  
  
  for (i in 1:(length(sera_shuffled$Sample_ID))/2 ){
    
    if (length(IDs) == 0) {
      break
    }
    indiv1 = sample(IDs, 1)
    
    IDs = IDs[ !IDs == indiv1]
    indiv1_age = sera_shuffled[sera_shuffled$Sample_ID == indiv1,]$Age  
    age_window = sera_shuffled[sera_shuffled$Age >= indiv1_age - 3 & 
                                 sera_shuffled$Age <= indiv1_age + 3, ]$Sample_ID
    age_window = intersect(age_window, IDs)
    
    try(
      {
        if (length(age_window) == 0) {
          continue
        }
        indiv2 = sample(age_window, 1)
        IDs = IDs[ !IDs == indiv2]
        indiv2_age = sera_shuffled[sera_shuffled$Sample_ID == indiv2,]$Age  
        
        titer1 = as.numeric(sera_shuffled[sera_shuffled$Sample_ID == indiv1 , c(3:10)])
        titer2 = as.numeric(sera_shuffled[sera_shuffled$Sample_ID == indiv2 , c(3:10)])
        
        cosine_similarity = cosine(titer1, titer2)
        age = (indiv1_age + indiv2_age)/2
        cosine_similarities_shuffled = c(cosine_similarities_shuffled, cosine_similarity)
        ages_shuffled = c(ages_shuffled, age)
      },silent=T
    )
    
  }
  
  df_shuffled = data.frame(cosine_similarities_shuffled, ages_shuffled)
  model_shuffled = lm(cosine_similarities_shuffled ~ ages_shuffled, data=df_shuffled)
  pred_shuffled = predict(model_shuffled, newdata = data.frame(ages_shuffled = xseq), se=TRUE)
  y_shuffled = pred_shuffled$fit
  predicts_shuffled = cbind(predicts_shuffled, y_shuffled)
  
  sm_fit_shuffled = summary(model_shuffled)
  slop_shuffled = sm_fit_shuffled$coefficients[2,1]
  slops_shuffled = c(slops_shuffled, slop_shuffled)
  p = sm_fit_shuffled$coefficients[2,4]
  pvals_shuffled = c(pvals_shuffled, p)
  print (r)
}


slops_shuffled = sort(slops_shuffled)
write.csv(slops_shuffled, "../result/slops_shuffled_lm_no_replacement_1000_age_window_3.csv")
slops_shuffled[25]
slops_shuffled[975]

saveRDS(predicts_shuffled, "../result/predicts_shuffled_lm_no_replacement_1000_age_window_3.rds")

#mean of splines from randomized data
lm.MEAN_shuffled = data.frame (xseq, rowMeans(predicts_shuffled[,c(2:1001)]))
colnames(lm.MEAN_shuffled) = c("x", "y")

#CI of splines from randomized data
predicts_sorted_shuffled = t(apply(predicts_shuffled,1,sort))
lm.CI_shuffled = data.frame(xseq, predicts_sorted_shuffled[,25+1], predicts_sorted_shuffled[,975+1])
colnames(lm.CI_shuffled) = c("x", "ci1", "ci2")


p_shuffled_lm_1000 = ggplot() +
  geom_ribbon(data = lm.CI_shuffled, aes(ymin = ci1, ymax = ci2, x=xseq), fill = "grey", alpha=0.5) +
  geom_smooth(aes_auto(lm.MEAN_shuffled), data=lm.MEAN_shuffled, stat="identity", col="black") +
  cosine_theme + 
  ylim(c(0,1)) + xlim(c(1,90)) +
  xlab("Age (years)") + ylab("Cosine similarity") 
p_shuffled_lm_1000

ggsave("../figure/cosine_similarity_shuffled_lm_1000_age_window_3.png", height=2.7, width=4.5)



################################################
library(ggpubr)

#ggarrange(p_one, p_lm_1000, p_shuffled_lm_1000,
#          nrow=3, ncol=1,
#          labels=c("A", "B", "C"))
#ggsave("../figure/cosine_similarity_lm_age_window_3.png", height=8, width=4.5)




