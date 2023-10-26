library(reshape2)

#NA ELLA
sera_na = read.table("../../../Data/sera/data/BFRNT dataset for analysis_v2.csv", sep=",", 
                     head=T, stringsAsFactors = F)
ella = data.frame(sera_na[, c("Sample_ID", "Age")], sera_na[, grep("ELLA", colnames(sera_na))])
ella = na.omit(ella)
colnames(ella) = gsub("ELLA_", "", colnames(ella))
ella$`3c2.A` = log2(ella$`3c2.A`/10)
ella$A2 = log2(ella$A2/10)
m_ella = melt(ella, id.vars = c("Sample_ID", "Age"))
colnames(m_ella)[grep("variable", colnames(m_ella))] = c("Test")
colnames(m_ella)[grep("value", colnames(m_ella))] = c("Titer")

colnames(m_ella) = gsub("Test", "Test_virus", colnames(m_ella))

df_na = m_ella
df_na = na.omit(df_na)
df_na$is_NA = 1
df_na = transform(df_na, Age_group = cut(Age, breaks = c(0, 4, 17, 44, 64, 90)))

df_na$Test_virus = factor( df_na$Test_virus, levels = c("3c2.A", "A2"))









