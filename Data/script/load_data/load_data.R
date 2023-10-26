
source("../../../Data/script/load_data/convert_sera.R")


load_frnt = function(file, prefix, rmvdup = 1){
  sera = read.table(file, sep=",", head=T, stringsAsFactors = F)
  sera = data_cleaning(sera, rmvdup)
  sera = log_titer_varname(sera, prefix)
  sera = change_virus_names(sera)
  
  sera
}

melt_frnt = function(sera){
  m_sera = melt(sera, id.vars = c("Sample_ID", "Age"))
  tv_id = grep("variable", colnames(m_sera))
  t_id = grep("value", colnames(m_sera))
  colnames(m_sera)[c(tv_id, t_id)] = c("Test_virus", "Titer")
  
  m_sera
}

load_and_merge_mutations = function(file, m_sera){
  #Import x variables: mutations in test virus
  xvar_muts = read.table(file, sep=",", head=T)
  vars = merge(m_sera, xvar_muts, by="Test_virus")
  
  vars
} 


load_and_merge_xvars = function(file, sera_mutations){
  imp = read.table(file, head=T, sep=",")
  imp$Age = 2017 - imp$y
  imp$Test_virus = imp$Virus
  
  #Merge input data
  vars_sorted_imp = merge(sera_mutations, imp, by = c("Test_virus", "Age"))
  
  vars_sorted_imp
}

#####################################################################
#filenames

seraf1 = "../../../Data/sera/data/BFRNT_for1718.csv"
seraf2 = "../../../Data/sera/data/FRNT_WT_clade3.csv"

#prefix in the column names to recognize test strain
prefix1 = "ND90"
prefix2 = "FRNT"


sera1 = load_frnt(seraf1, prefix1, rmvdup = 0)
sera2 = load_frnt(seraf2, prefix2, rmvdup = 0)
sera = merge(sera1, sera2, by=c("Sample_ID", "Age"), all=T)

#####################################################################
# WT and WT2 are titers that are obtained from testing the same serum samples to the same reference strain(3C.2A) at different time
# WT2 were obtained when 3C.2A2-2 were tested, to confirm that there is no difference in titers by testing time (WT vs. WT2)
# Move WT2 to WT if matched (WT == WT2). If not matched, remove that replicate, as Sigrid did.

sera$match =  ifelse(sera$WT == sera$WT2, 1, 0)

#test if titer is different by time using serum samples tested to WT at different time points (WT v. WT2)
#timetest = na.omit(sera)
#t.test(timetest$WT, timetest$WT2)
#sum(sera$match == 1, na.rm=T)
#sum(sera$match == 0, na.rm=T)
#unmatched = sera[sera$match == 0, c("WT", "WT2")]
#unmatched = na.omit(unmatched)

sera[is.na(sera$match), "match"] = 1 #match is na if there is only WT or WT2: those individuals should be kept
sera$WT = ifelse((sera$match == 1 & !is.na(sera$WT2)), sera$WT2, sera$WT)
sera  = sera[, -grep("WT2", colnames(sera))]

sera = sera[sera$match == 1,]
sera = sera[,-grep("match", colnames(sera))]

##############################################################################
# There are some replicated IDs that have two sets of full vector of titers
# Use first replicates when there are duplicated individuals
# For testing the effect of test time, duplicated individuals are used
sera = sera[!duplicated(sera$Sample_ID), ]


m_sera = melt(sera, 
              id.vars = c("Sample_ID", "Age"),
              variable.name = "Test_virus",
              value.name = "Titer")
m_sera = na.omit(m_sera)

df = m_sera

df = transform(df, Age_group = cut(Age, breaks = c(0, 4, 17, 44, 64, 90)))

################################################################################
#merge sera, mutation and covariates

df = na.omit(df)

#sort df by Sample_ID
#df = df[order(df[,"Sample_ID"]), ]

# this data is for HA, not for NA

df$is_NA = 0

# Change virus names for data frame and sera
#unique(df$Test_virus)
df$Test_virus = gsub( "3C.3a", "3C3.A", df$Test_virus )
df$Test_virus = gsub( "WT", "3C2.A", df$Test_virus )
df$Test_virus = factor( df$Test_virus, levels = c("3C3.A", "3C2.A", "N171K", "N121K_N171K", "N121K_T135K_N171K",
                                                  "T131K_R142K", "T131K_R142K_R261Q", "N121K_S144K"))

colnames(sera) = gsub("3C.3a", "3C3.A", colnames(sera))
colnames(sera) = gsub("WT", "3C2.A", colnames(sera))


