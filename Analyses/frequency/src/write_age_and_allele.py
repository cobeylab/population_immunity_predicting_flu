import csv
import numpy
import os
import re
from copy import deepcopy as deepcopy
from age_distribution_by_clade_functions import sequences_of_the_season 
from age_distribution_by_clade_functions import metadata_age 
from age_distribution_by_clade_functions import metadata_age_state 
from age_distribution_by_clade_functions import write_age 
from age_distribution_by_clade_functions import write_clade
from age_distribution_by_clade_functions import write_age_loc
from age_distribution_by_clade_functions import metadata_age_allele
from age_distribution_by_clade_functions import write_age_allele
from age_distribution_by_clade_functions import write_age_labs
from age_distribution_by_clade_functions import metadata_origin
#from new_variant import allele_frequency_increase

from define_virus import *

HAclades = [c2a, A1, A1_2, A1_3, A2, A2_2, A3, c3a]
HAnames = ["c2a", "A1", "A1_2", "A1_3", "A2", "A2_2", "A3", "c3a"]
	
NAclades = [anyNA, gly245, nonGly245, not245N380V, N329S_K220N_V303I, N329S, N329S_I176M, N329S_I176M_P386S,  L140I, D339N, others, P468L_N161S]
NAnames = ["anyNA", "gly245", "nonGly245", "not245N380V", "NA_N329S_K220N_V303I", "NA_N329S", "NA_N329S_I176M",
 "NA_N329S_I176M_P386S",  "NA_L140I", "NA_D339N", "NA_others", "NA_P468L_N161S"]
 

NA_cids = [[] for i in range(len(NAclades))]


############################################################################################################
#write age for each clade
############################################################################################################


full_season = 1

outf_dir = "../data_clade_assigned/"
outf_dir_NA = "../data_clade_assigned/NA_"
outf_dir_US = "../data_clade_assigned/US_"
outf_dir_NE = "../data_clade_assigned/NE_"
outf_dif_allele = "../data_allele/"

#HA
segment = "HA"
sites = HAsites
mClades = HAclades
cNames = HAnames


# season = 2017-18
season = 2017
HAf_name = "../data/gisaid_1718_submitted_upto_20220111_align_cut_std_AA.fas"metaf_name = "../meta/meta_1718_submitted_upto_20220111.csv"

clades_17, seqs_17, numSeqs_17 = sequences_of_the_season(season, sites, HAf_name, full_season, mClades)

ages_17, ids_17, dates_17 = metadata_age_state(season, clades_17, metaf_name, segment, NE)
write_age(season, ids_17, ages_17, dates_17, cNames + ['other'], outf_dir_NE)

ages_17, ids_17, dates_17, locs_17 = metadata_age(season, clades_17, metaf_name, segment, USCanada=1)
write_age_loc(season, ids_17, ages_17, dates_17, cNames + ['other'], locs_17, outf_dir_NA)

ages_17, ids_17, dates_17, locs_17 = metadata_age(season, clades_17, metaf_name, segment, USCanada=0)
write_age_loc(season, ids_17, ages_17, dates_17, cNames + ['other'], locs_17, outf_dir_US)


# season = 2016-17
season = 2016
HAf_name = "../data/gisaid_1617_submitted_upto_20220111_align_cut_std_AA.fas"
metadir = "../meta/meta_1617_submitted_upto_20220111.csv"

clades_16, seqs_16, numSeqs_16 = sequences_of_the_season(season, sites, HAf_name, full_season, mClades)	

ages_16, ids_16, dates_16  = metadata_age_state(season, clades_16, metadir, segment, NE)
write_age(season, ids_16, ages_16, dates_16, cNames+ ['other'], outf_dir_NE)

ages_16, ids_16, dates_16, locs_16  = metadata_age(season, clades_16, metadir, segment, USCanada=1)
write_age_loc(season, ids_16, ages_16, dates_16, cNames+ ['other'], locs_16, outf_dir_NA)

ages_16, ids_16, dates_16, locs_16  = metadata_age(season, clades_16, metadir, segment, USCanada=0)
write_age_loc(season, ids_16, ages_16, dates_16, cNames+ ['other'], locs_16, outf_dir_US)


