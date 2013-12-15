from Bio.Blast import NCBIXML
from glob import glob
import sets
import pdb

# read blast results from both species

s2m = {} # solenopsis matching monomorium, with solenopsis genes as keys
for file in glob("data/sinv_blast/sinv_groups/*xml"):
	infile = open(file)
	for rec in NCBIXML.parse(infile):
		if rec.alignments:
			s2m[str(rec.query)] = str(rec.alignments[0].hit_def).split()[0]
	infile.close()

m2s = {} # solenopsis matching monomorium, with solenopsis genes as keys
for file in glob("data/sinv_blast/mp_groups/*xml"):
	infile = open(file)
	for rec in NCBIXML.parse(infile):
		if rec.alignments:
#			pdb.set_trace()
			if str(rec.alignments[0].hit_def) in s2m and \
			s2m[str(rec.alignments[0].hit_def)] == str(rec.query).split()[0]:
				m2s[str(rec.alignments[0].hit_def)] = str(rec.query).split()[0]
	infile.close()

# combine results to find reciprocal best blast matches
common = set(s2m.keys()) & set(m2s.keys())


orthologs = open("data/sinv_blast/homologs.txt","w") #orthologous gene pairs
for gene in common:
	if s2m[gene] == m2s[gene]:
		orthologs.write("%s\t%s\n" % (gene, m2s[gene]))
orthologs.close()

