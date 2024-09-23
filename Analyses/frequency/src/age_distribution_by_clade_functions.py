import csv
import numpy
import os
import re
from copy import deepcopy as deepcopy

def determine_season(y, m, full_season):
	if full_season == 1:
		if m < 10:
			season = y-1
		else:
			season = y
		return season
	elif full_season == 0:
		if m <= 3:
			season = y-1
		elif m > 3 and m < 10:
			return None
		else:
			season = y
		return season
		
def determine_clade(sites, seq, id, seqs, mClades):
	segs = ''
	all_c = []
	
	for s in sites:
		segs += seq[s-1]
	for c in range(len(mClades)):
		if re.match(mClades[c], segs) != None:

			return c

	return -1

	
def determine_clade_c(sites, seq, id, clades, seqs, mClades):
	segs = ''
	all_c = []
	for s in sites:
		segs += seq[s-1]
	for c in range(len(mClades)):
		if re.match(mClades[c], segs) != None:
			#clades[c].append(id)
			#seqs[c].append(seq)
			all_c.append(c)

	return all_c
######################################################
#sequences
	
def sequences_of_the_season(season, sites, inf_name, full_season, mClades):
	clades = [[] for i in range(len(mClades))] + [[]]
	seqs = [[] for i in range(len(mClades))] + [[]]


	i_id = 0
	i_date = 2

	numSeq = 0
	inf = open(inf_name, "rU")
	for line in inf:
		if line.find(">") >= 0:
			each = line.split("\n")[0].split(">")[1].split("|")
			id = each[i_id]
			date = each[i_date]
			try:
				y = int(date.split("/")[0])
			except ValueError:
				y = ''
				continue
			try:
				m = int(date.split("/")[1])
			except ValueError:
				try:
					m = int(date.split("/")[1].split("_")[0])
				except ValueError:
					m = ''
					continue
			s = determine_season(y, m, full_season)
		else:
			if y=='' or m=='':
				continue

			seq = line
			if s == season:
				numSeq += 1

				c = determine_clade(sites, seq, id, seqs, mClades)
				clades[c].append(id)
				seqs[c].append(seq)
				
	return clades, seqs, numSeq

def sequences_into_clade_with_date(season, sites, inf_name, full_season, mClades):
	clades = [[] for i in range(len(mClades))]
	seqs = [[] for i in range(len(mClades))]
	dates = [[] for i in range(len(mClades))]
	
	i_id = 0
	i_date = 2

	numSeq = 0
	inf = open(inf_name, "rU")
	for line in inf:
		if line.find(">") >= 0:
			each = line.split("\n")[0].split(">")[1].split("|")
			id = each[i_id]
			date = each[i_date]
			try:
				y = int(date.split("/")[0])
			except ValueError:
				continue
				
			try:
				m = int(date.split("/")[1])
			except ValueError:
				try:
					m = int(date.split("/")[1].split("_")[0])
				except ValueError:
					m = ''
					continue
					
			try:
				d = int(date.split("/")[2])
			except ValueError:
				d = ''
				continue
			except IndexError:
				d = ''
				continue
				
			s = determine_season(y, m, full_season)
		else:
			if m=='' or d=='':
				continue

			seq = line
			if s == season:
				numSeq += 1
				all_c = determine_clade_c(sites, seq, id, clades, seqs, mClades)
				for c in all_c:
					clades[c].append(id)
					seqs[c].append(seq)
					dates[c].append(date)
				
	return clades, seqs, dates
	
######################################################
# Metadata

#use most updated metadata for season 2017
	
def metadata_age(season, clades, metaf_name, segment, USCanada):

	c_age = "Host_Age"
	c_au = "Host_Age_Unit"

	c_loc = "Location"
	c_date = "Collection_Date"
	
	if segment == "HA" :

		c_id = "HA Segment_Id"
	elif segment == "NA" :
		c_id = "NA Segment_Id"


	ages = [[] for i in range(len(clades)+1)]
	ids = [[] for i in range(len(clades)+1)]
	dates = [[] for i in range(len(clades)+1)]

	locs = [[] for i in range(len(clades)+1)]

	count = 0
	
	#for file in os.listdir(metadir):
	#	if file.endswith(".csv"):
	#		metaf_name = os.path.join(metadir, file)
	metaf = open(metaf_name, "rU")
		
	reader = csv.DictReader(metaf)
	for row in reader:
		try:

			id = row[c_id].split(" | ")[0].split("EPI")[1]

			loc = row[c_loc]
			coldat = row[c_date]
			
			if USCanada == 1:
				if loc.find("United States") < 0 and loc.find("Canada") < 0:
					continue
			elif USCanada == 0:
				if loc.find("United States") < 0:
					continue

			try:
				age = int(float(row[c_age]))
				au = row[c_au]
				if au == 'M':
					age = age//12
			except ValueError:

				#print (row[c_age])
				age = "NA"
				#continue
			
		except IndexError:
			
			continue
		

		for c in range(len(clades)):
			if id in clades[c]:
				ages[c].append(str(age))
				ids[c].append(id)
				dates[c].append(coldat)

				locs[c].append(loc)


	return ages, ids, dates, locs
	
def metadata_age_prevSeasons(season, clades, metadir, segment, USCanada):

	c_age = "Host_Age"
	c_au = "Host_Age_Unit"

	c_loc = "Location"
	c_date = "Collection_Date"
	
	if segment == "HA" :

		c_id = "HA Segment_Id"
	elif segment == "NA" :
		c_id = "NA Segment_Id"


	ages = [[] for i in range(len(clades)+1)]
	ids = [[] for i in range(len(clades)+1)]
	dates = [[] for i in range(len(clades)+1)]
	count = 0

	
	for file in os.listdir(metadir):
		if file.endswith(".csv"):
			metaf_name = os.path.join(metadir, file)
			metaf = open(metaf_name, "rt", encoding = 'utf8')
				
			reader = csv.DictReader(metaf)
			for row in reader:
				try:

					id = id = row[c_id].split(" | ")[0].split("EPI")[1]
					loc = row[c_loc]
					coldat = row[c_date]
					
					if USCanada == 1:
						if loc.find("United States") < 0 and loc.find("Canada") < 0:
							continue

					try:
						age = int(float(row[c_age]))
						au = row[c_au]
						if au == 'M':
							age = age//12
					except ValueError:
						continue
					
				except IndexError:
					continue
				

				for c in range(len(clades)):
					if id in clades[c]:
						ages[c].append(str(age))
						ids[c].append(id)
						dates[c].append(coldat)


	return ages, ids, dates
	
	
def metadata_age_state(season, clades, metaf_name, segment, states):

	c_age = "Host_Age"
	c_au = "Host_Age_Unit"

	c_loc = "Location"
	c_date = "Collection_Date"
	
	if segment == "HA" :
		c_id = "HA Segment_Id"
	elif segment == "NA" :
		c_id = "NA Segment_Id"


	ages = [[] for i in range(len(clades))]
	ids = [[] for i in range(len(clades))]
	dates = [[] for i in range(len(clades))]
	
	locs = []
	#for file in os.listdir(metadir):	#	if file.endswith(".csv"):
			#metaf_name = os.path.join(metadir, file)
	metaf = open(metaf_name, "rU")
		
	reader = csv.DictReader(metaf)
	for row in reader:
		try:
			id = id = row[c_id].split(" | ")[0].split("EPI")[1]		
			loc = row[c_loc]
			state = row[c_loc].split(" / ")[2]
			coldat = row[c_date]
			
			if loc.find("United States") < 0 :
				continue
			
			
			if state not in states:
				continue
			locs.append(state)

			try:
				age = int(float(row[c_age]))
				au = row[c_au]
				if au == 'M':
					age = age//12
			except ValueError:

				age = "NA"
				#continue

			
		except IndexError:
			continue
		
		for c in range(len(clades)):
			if id in clades[c]:
				ages[c].append(str(age))
				ids[c].append(id)
				dates[c].append(coldat)

	print (set(locs))
	return ages, ids, dates
	
	
def metadata_origin(season, clades, metadir, segment):

	c_age = "Host_Age"
	c_au = "Host_Age_Unit"

	c_loc = "Location"
	c_date = "Collection_Date"
	
	c_origin = "Originating_Lab"
	
	if segment == "HA" :
		c_id = "HA Segment_Id"
	elif segment == "NA" :
		c_id = "NA Segment_Id"

	ages = [[] for i in range(len(clades))]
	ids = [[] for i in range(len(clades))]
	labs = [[] for i in range(len(clades))]
	
	for file in os.listdir(metadir):
		if file.endswith(".csv"):
			metaf_name = os.path.join(metadir, file)
			metaf = open(metaf_name, "rU")
				
			reader = csv.DictReader(metaf)
			for row in reader:
				try:
					id = row[c_id].split(" | ")[0].split("EPI")[1]
					loc = row[c_loc]
					lab = row[c_origin]
					
					if loc.find("United States") < 0 and loc.find("Canada") < 0:
						continue

					try:
						age = int(float(row[c_age]))
						au = row[c_au]
						if au == 'M':
							age = age//12
					except ValueError:
						continue
					
				except IndexError:
					continue
				
				for c in range(len(clades)):
					if id in clades[c]:
						ages[c].append(str(age))
						ids[c].append(id)
						labs[c].append(lab)

	return ages, ids, labs

def metadata_age_allele(season, clades, seqs, metadir, increased, segment):
	c_age = "Host_Age"
	c_au = "Host_Age_Unit"

	c_loc = "Location"
	c_date = "Collection_Date"
	
	if segment == "HA" :
		c_id = "HA Segment_Id"
	elif segment == "NA" :
		c_id = "NA Segment_Id"
		
	ages = [[] for i in range(len(clades))]
	ids = [[] for i in range(len(clades))]
	alleles = [[] for i in range(len(clades))]
	dates = [[] for i in range(len(clades))]
		
	for file in os.listdir(metadir):
		if file.endswith(".csv"):
			metaf_name = os.path.join(metadir, file)
			metaf = open(metaf_name, "rU")
				
			reader = csv.DictReader(metaf)
			for row in reader:
				try:
					id = row[c_id].split(" | ")[0].split("EPI")[1]
					loc = row[c_loc]
					date = row[c_date]
					
					if loc.find("United States") < 0 and loc.find("Canada") < 0:
						continue

					try:
						age = int(float(row[c_age]))
						au = row[c_au]
						if au == 'M':
							age = age//12
					except ValueError:
						continue
					
				except IndexError:
					continue
				
				
				for c in range(len(clades)):
					if id in clades[c]:
						ages[c].append(str(age))
						ids[c].append(id)
						dates[c].append(date)
						
						idx = clades[c].index(id)
						
						#alleles in segregating sites
						al = []
						for a in range(len(increased[c])):
							s = increased[c][a][0] - 1
							this_allele = seqs[c][idx][s]
							
							if this_allele == increased[c][a][1]:
								state = "ances_"
							elif this_allele == increased[c][a][2]:
								state = "mut_"
							else:
								state = "other_"
							al.append(state + seqs[c][idx][s])
						alleles[c].append(al)

	return ages, ids, alleles, dates
	
######################################################		


#def write_age(season, ids, ages, dates, cNames, outf_dir):
#	for c in range(len(cNames)):
#		outf = open(outf_dir + cNames[c]+"_season"+str(season) +".csv", "w")
#		outf.write("id,age,collection,clade,season\n")
#		pSeason = "Season" + str(season) + str(season+1)
#		for a in range(len(ages[c])):
#			outf.write(ids[c][a] + "," + ages[c][a] + "," + dates[c][a] + "," + cNames[c] + "," + pSeason + "\n")
#			
#		outf.close()
		
		
def write_age(season, ids, ages, dates, cNames, outf_dir):

	outf = open(outf_dir + "clade_assigned_season_"+str(season) +".csv", "w")
	outf.write("id,age,collection,clade,season\n")
	
	for c in range(len(cNames)):
		
		pSeason = "Season" + str(season) + str(season+1)
		for a in range(len(ages[c])):
			outf.write(ids[c][a] + "," + ages[c][a] + "," + dates[c][a] + "," + cNames[c] + "," + pSeason + "\n")
			

	outf.close()
	
def write_age_loc(season, ids, ages, dates, cNames, locs, outf_dir):

	outf = open(outf_dir + "clade_assigned_season_"+str(season) +".csv", "w")
	outf.write("id,age,collection,clade,location,season\n")
	
	for c in range(len(cNames)):
		
		pSeason = "Season" + str(season) + str(season+1)
		for a in range(len(ages[c])):
			outf.write(ids[c][a] + "," + ages[c][a] + "," + dates[c][a] + "," + cNames[c] + "," + locs[c][a].split("/")[0] +","+ pSeason + "\n")
			
	outf.close()

	
def write_clade(season, ids, cNames,  outf_dir):
	outf = open(outf_dir + "clade_assigned_season_"+str(season) +".csv", "w")
	outf.write("id,age,collection,clade,location,season\n")
	
	for c in range(len(cNames)):
		
		pSeason = "Season" + str(season) + str(season+1)
		for a in range(len(ids[c])):
			outf.write(ids[c][a] + "," + "NA" + "," + "NA" + "," + cNames[c] + "," + "NA" +","+ pSeason + "\n")
			
	outf.close()



def write_date(season, clades_17, dates_17, cNames, outf_dir_date, full_season):
	outf = open(outf_dir_date + "clades_with_dates_" + str(season) +"_" + str(full_season) + ".csv", "w")
	outf.write("id,clade,date,season\n")
	for c in range(len(cNames)):
		for i in range(len(clades_17[c])):
			oneline = clades_17[c][i] + "," + cNames[c] + "," + dates_17[c][i] +"," + str(season)
			outf.write(oneline+"\n")
	outf.close()

def write_age_labs(season, ids, ages, labs, cNames, outf_dir):
	for c in range(len(cNames)):
		outf = open(outf_dir + cNames[c]+"_season"+str(season) +".csv", "w")
		outf.write("id,age,clade,labs,season\n")
		pSeason = "Season" + str(season) + str(season+1)
		for a in range(len(ages[c])):
			outf.write(ids[c][a] + "," + ages[c][a] + "," + cNames[c] + "," + labs[c][a] + "," + pSeason + "\n")
			
		outf.close()
		
def write_age_allele(season, ids, ages, increased, alleles, dates, cNames, outf_dir):
	for c in range(len(cNames)):
		outf = open(outf_dir + "alleles_" + cNames[c]+"_season"+str(season) +".csv", "w")
		head_alleles = ''
		for i in range(len(increased[c])):
			head_alleles += str(increased[c][i][0]) + ","
		outf.write("id,age,clade,"+head_alleles+"season,date\n")
		pSeason = "Season" + str(season) + str(season+1)
		for a in range(len(ages[c])):
			outf.write(ids[c][a] + "," + ages[c][a] + "," + cNames[c] + ",")
			val_alleles = ''
			for i in range(len(alleles[c][a])):
				val_alleles += alleles[c][a][i] + ","
			outf.write(val_alleles + pSeason + "," + dates[c][a] + "\n")
			
		outf.close()		
######################################################	
	