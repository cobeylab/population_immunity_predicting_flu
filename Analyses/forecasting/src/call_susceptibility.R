######################################################################
#The user needs to set options as below.
######################################################################
# opt_na
# opt_na = 0: HA titers
# opt_na = 1: NA titers
#
# threshold
# threshold = 0: using GMT to calculate relative susceptibility
# threshold = 1: using cutoff of 1:20 to calculate relative susceptibility
# threshold = 2: using cutoff of 1:40 
# threshold =3: using cutoff of 1:80 
# threshold = 4: using cutoff of 1:160 
# 
# opt_ageGroup
# opt_ageGroup = 0: using all population
# opt_ageGroup = 1: 1-4 years old 
# opt_ageGroup = 2: 5-17 years old 
# opt_ageGroup = 3: 18-44 years old
# opt_ageGroup = 4: 45-64 years old
# opt_ageGroup = 5: 65-90 years old


#################################################################
# Parameters


if("dplyr" %in% (.packages())){
  detach("package:dplyr", unload=TRUE) 
  detach("package:plyr", unload=TRUE) 
} 

library(dplyr)
fig_dir = "../fig/"
result_dir = "../result/" 

group = c("bin")

n_repl = 1000 #number of repeats for bootstrapping
opt_na = 0 #opt_na = 0: ha; opt_na = 1: na

##########################################################


if(opt_na==0){
  # set the order of reference viruses for analysis and result
  viruses = c("3C3.A", "3C2.A", "N171K", "N121K_N171K", "N121K_T135K_N171K",
              "T131K_R142K", "T131K_R142K_R261Q", "N121K_S144K")
  
}else if (opt_na==1){

  # set the order of reference viruses for analysis and result
    viruses = c("3c2.A", "A2")
  
}

###########################################################
#threshold 0: using GMT to calculate relative susceptibility


threshold = 0
opt_ageGroup = 0
source("susceptibility.R")

threshold = 0
opt_ageGroup = 1
source("susceptibility.R")

threshold = 0
opt_ageGroup = 2
source("susceptibility.R")

threshold = 0
opt_ageGroup = 3
source("susceptibility.R")

threshold = 0
opt_ageGroup = 4
source("susceptibility.R")

threshold = 0
opt_ageGroup = 5
source("susceptibility.R")



###########################################################
# threshold = 1: using cut off 1:20

threshold = 1
opt_ageGroup = 0
#source("susceptibility.R")

threshold = 1
opt_ageGroup = 1
#source("susceptibility.R")

threshold = 1
opt_ageGroup = 2
#source("susceptibility.R")

threshold = 1
opt_ageGroup = 3
#source("susceptibility.R")

threshold = 1
opt_ageGroup = 4
#source("susceptibility.R")

threshold = 1
opt_ageGroup = 5
#source("susceptibility.R")

###########################################################
# threshold = 2: using cut off 1:40

threshold = 2
opt_ageGroup = 0
source("susceptibility.R")

threshold = 2
opt_ageGroup = 1
source("susceptibility.R")

threshold = 2
opt_ageGroup = 2
source("susceptibility.R")

threshold = 2
opt_ageGroup = 3
source("susceptibility.R")

threshold = 2
opt_ageGroup = 4
source("susceptibility.R")

threshold = 2
opt_ageGroup = 5
source("susceptibility.R")

###########################################################
# threshold = 3: using cut off 1:80

threshold = 3
opt_ageGroup = 0
#source("susceptibility.R")

threshold = 3
opt_ageGroup = 1
#source("susceptibility.R")

threshold = 3
opt_ageGroup = 2
#source("susceptibility.R")

threshold = 3
opt_ageGroup = 3
#source("susceptibility.R")

threshold = 3
opt_ageGroup = 4
#source("susceptibility.R")

threshold = 3
opt_ageGroup = 5
#source("susceptibility.R")


###############################################
# threshold = 4: using cutoff 1:160

threshold = 4
opt_ageGroup = 0
source("susceptibility.R")

threshold = 4
opt_ageGroup = 1
source("susceptibility.R")

threshold = 4
opt_ageGroup = 2
source("susceptibility.R")

threshold = 4
opt_ageGroup = 3
source("susceptibility.R")

threshold = 4
opt_ageGroup = 4
source("susceptibility.R")

threshold = 4
opt_ageGroup = 5
source("susceptibility.R")


#######################################

