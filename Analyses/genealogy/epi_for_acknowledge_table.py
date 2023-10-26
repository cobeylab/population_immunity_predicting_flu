inf = open("sample_for_genealogy_1217_clade.csv", "rU")
outf = open("epi_for_acknowledge_table_genealogy.txt", "w")

epi_ids = ''
for line in inf:
	epi_id = "EPI" + line.split(",")[0].split("|")[0] + " "
	epi_ids = epi_ids + epi_id
	
inf.close()

inf = open("na_sample_for_genealogy_1217_clade.csv", "rU")
for line in inf:
	epi_id = "EPI" + line.split(",")[0].split("|")[0] + " "
	epi_ids = epi_ids + epi_id
	
inf.close()

outf.write(epi_ids)
outf.close()
