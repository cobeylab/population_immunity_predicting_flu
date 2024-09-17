
from random import sample


#meta files
meta_files = ["../frequency/meta/2012_meta_HA.csv", "../frequency/meta/2013_meta_HA.csv", "../frequency/meta/2014_meta_HA.csv",
				"../frequency/meta/2015_meta_HA.csv", "../frequency/meta/2016_meta_HA.csv", 
				"../frequency/meta/meta_1617_submitted_upto_20200109.csv", "../frequency/meta/meta_1718_submitted_upto_20200109.csv"]

#HA fasta files 
files_ha = ["../frequency/data/gisaid_1217_submitted_upto_20170605_align_cut_std.fas",
			"../frequency/data/gisaid_1617_submitted_upto_20200109_align_cut_std.fas", 
			"../frequency/data/gisaid_1718_submitted_upto_20200109_align_cut_std.fas"]
			

#NA fasta files  
files_na = ["../frequency/data/gisaid_na_1216_align_cut_std.fas",
			"../frequency/data/gisaid_na_1617_sub_upto_20200306_align_cut_std.fas", 
			"../frequency/data/gisaid_na_1718_sub_upto_20200306_align_cut_std.fas"]


####################################################################################################

#save all meta data of viruses with HA and NA sample
def get_good_strain_names (meta_files):
	
	good_names = []
	for file in meta_files:
		f = open(file, "r")
		
		for line in f:
			if line.find("Isolate_Id") >= 0:
				continue
			each = line.split(",")
			if len(each) < 10:
				continue
			if (("|" in each[4]) & ("|" in each[6])):
				name = each[4].split(" | ")[1]
				good_names.append(name)
				
	return good_names
	

#put sequences into list by season
def seq_by_season(files, s_start, s_end, good_names, seasons):
	for file in files:

		f = open(file, "r")

		for line in f:
						
			if line.find(">") >= 0:
				found = 0 
				ymd = line.split("|")[2]
				ymd = ymd.split("/")
				
				name = line.split("|")[1]
				if name not in good_names:
					continue
					
				try:
					y = ymd[0]
					m = ymd[1]
					d = ymd[2]
					
					if int(m)>=10:
						season = int(y) + 1
					elif int(m)<10:
						season = int(y)

					season = season - 2013
					found = 1
					head = line
				except IndexError:
					continue
				
			else:
				if found == 1:
					
					line = line.split("\n")[0]
					ambiguous = 0
					for s in line:
						if s != 'a' and s!= 'c' and s!= 'g' and s!='t':
							ambiguous = 1
					
					if ambiguous == 1:
						continue
					if line.find("-") >= 0 or line.find("?") >= 0:
						continue
					else:
						if season >= s_start and season <= s_end:
							seasons[season].append(head+line+"\n")
							
	return seasons
						
						
def sample_ha_na (num_samp, season_ha, season_na):
	
	ha_names = []
	for seq in season_ha:
		name = seq.split("|")[1]
		ha_names.append(name)
		
	na_names = []
	for seq in season_na:
		name = seq.split("|")[1]
		na_names.append(name)
		
	samples_ha = []
	samples_na = []
	
	sample_names = []
	for n in range(0, num_samp):
		ha_sample_name = sample(ha_names, 1)[0]
		while (ha_sample_name not in na_names): #for selected ha, na should also be selected
			ha_sample_name = sample(ha_names, 1)[0] #if there is no na sequence, sample ha sequence again
			
		while (ha_sample_name in sample_names): #if selected ha sequenc is duplicated, sample ha sequence again
			ha_sample_name = sample(ha_names, 1)[0]
			while (ha_sample_name not in na_names):
				ha_sample_name = sample(ha_names, 1)[0]
			
		sample_names.append(ha_sample_name) 
		
		ha_idx = ha_names.index(ha_sample_name)
		na_idx = na_names.index(ha_sample_name)
		
		samples_ha.append(season_ha[ha_idx])
		samples_na.append(season_na[na_idx])
		
		
	return samples_ha, samples_na
	


def write_fasta_and_date_file (outf_name, datef_name, samples):
	
	outf = open(outf_name, "w")
	datef = open(datef_name, "w")		
	
	for s in range(len(samples)):
		
		for seq in samples[s]:
		
			head = seq.split("\n")[0]
			head = head.replace(";", "")
			head = head.replace(" ", "_")

			line = seq.split("\n")[1]
			
			seq = head + "\n" + line + "\n"
			outf.write(seq)

			name = head.split(">")[1]
			ymd = head.split("|")[2]
			datef.write(name + "\t" + ymd + "\n")

	outf.close()
	datef.close()


	
###########################################################################################


seasons_ha = [[], [], [], [], [], []]
seasons_na = [[], [], [], [], [], []]

print ("reading meta file")
good_names = get_good_strain_names (meta_files)

print ("saving HA sequences")
seasons_ha = seq_by_season(files_ha, 0, 5, good_names, seasons_ha)

print ("saving NA sequences")
seasons_na = seq_by_season(files_na, 0, 5, good_names, seasons_na)
						
print ("now sampling")
samples_ha_all = []
samples_na_all = []
for s in range(len(seasons_ha)):
	if s < 4:
		samples_ha, samples_na = sample_ha_na (20, seasons_ha[s], seasons_na[s])
		samples_ha_all.append(samples_ha)
		samples_na_all.append(samples_na)

	elif s >= 4:
		samples_ha, samples_na = sample_ha_na (100, seasons_ha[s], seasons_na[s])
		samples_ha_all.append(samples_ha)
		samples_na_all.append(samples_na)


########################################################################################
# Now write on file

print ("now writing")


ha16_outf_name = "../genealogy/sample_for_genealogy_1216.fasta"
ha16_datef_name = "../genealogy/sample_for_genealogy_1216_dates.txt"

na16_outf_name = "../genealogy/na_sample_for_genealogy_1216.fasta"
na16_datef_name = "../genealogy/na_sample_for_genealogy_1216_dates.txt"

ha17_outf_name = "../genealogy/sample_for_genealogy_1217.fasta"
ha17_datef_name = "../genealogy/sample_for_genealogy_1217_dates.txt"

na17_outf_name = "../genealogy/na_sample_for_genealogy_1217.fasta"
na17_datef_name = "../genealogy/na_sample_for_genealogy_1217_dates.txt"


#write_fasta_and_date_file(ha16_outf_name, ha16_datef_name, samples_ha_all[:-1])
#write_fasta_and_date_file(ha16_outf_name, ha16_datef_name, samples_ha_all[:-1])

#write_fasta_and_date_file(na16_outf_name, na16_datef_name, samples_na_all[:-1])
#write_fasta_and_date_file(na16_outf_name, na16_datef_name, samples_na_all[:-1])

write_fasta_and_date_file(ha17_outf_name, ha17_datef_name, samples_ha_all)
write_fasta_and_date_file(ha17_outf_name, ha17_datef_name, samples_ha_all)

write_fasta_and_date_file(na17_outf_name, na17_datef_name, samples_na_all)
write_fasta_and_date_file(na17_outf_name, na17_datef_name, samples_na_all)

