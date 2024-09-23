import re

############################################################################################################
# define clades : clades are exclusive from each other
############################################################################################################

#North East

Region1 = ["Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont"]
Region2 = ["New Jersey", "New York"]
Region3 = ["Delaware", "District of Columbia", "Maryland", "Pennsylvania", "Virginia", "West Virginia"]

NE =  Region1 + Region2 + Region3
#NA clades

NAsites = [468,161,  140,339,  329,  220,303,176,386,  245,380]

N329S_K220N_V303I = re.compile("....SNIIPN.") #A1b/135N, A1b/135K(A1_3)

N329S             = re.compile("....SKVIPN.")
N329S_I176M       = re.compile("....SKVMPN.") # some of basal nodes of A2/re 
N329S_I176M_P386S = re.compile("....SKVMSN.") #A2/re (it is also our test virus)
	
L140I             = re.compile("..INNKVIPN.") # Some of A1, Some of 3C3.A
D339N             = re.compile("..LNNKVIPN.")

others            = re.compile(".NLDN....N.")

P468L_N161S       = re.compile("LS.NN....N.") #A3

gly245            = re.compile(".........N.")
nonGly245         = re.compile(".........S.")
anyNA             = re.compile("...........")
not245N380V       = re.compile(".........S.")

#HA clades 

HAsites = [128, 131, 135, 138, 142, 144, 159, 160, 121, 261, 171] 

c2a = 	re.compile("TTTA.SYTNRN")

A1 = 	re.compile("TTTA.SYTNRK")
A1_2 = 	re.compile("TTTA.SYTKRK")
A1_3 = 	re.compile("TTKA.SYTKRK")

A2 = 	re.compile("TKTAKSYTNRN")
A2_2 = 	re.compile("TKTAKSYTNQN")

A3 = 	re.compile("TTTA.KYTKRN")
c3a = 	re.compile("ATTS..SKNRN")

















