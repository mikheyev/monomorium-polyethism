import pdb
"""
collect go terms from file, and convert them into sql-friendly csv
"""
infile = open("blast2go_export.txt")
goTable = open("go_table.csv","w")
goTable.write("isoform,GO\n")
blastTable = open("blast_table.csv","w")
blastTable.write("isoform,evalue,hit,species\n")

infile.readline() 
for line in infile:
	line=line.split("\t")
	blastTable.write(",".join(line[0:2]+['"'+line[2].split("[")[0].rstrip()+'"']+[line[3]])+"\n")
	for go in line[4].split(";"):
		goTable.write(",".join([line[0],go.strip()])+"\n")


