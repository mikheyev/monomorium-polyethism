from Bio import SeqIO
from glob import glob
import sets
import pdb

# read table of homologs
genes = {}
for line in open("data/sinv_blast/homologs.txt").readlines():
	line = line.rstrip().split()
	genes[line[0]] = line[2]

mp_genes = set(genes.values()) # set of mp orthologs

mp_prot = {}
mp_strand = {}
for file in glob("data/sinv_blast/mp_groups/*_prot.fa"):
	infile = open(file)
	for rec in SeqIO.parse(infile,"fasta"):
		if rec.id in mp_genes:
			mp_prot[rec.id] = rec
			mp_strand[rec.id] = rec.description.split()[1][0]
	infile.close()

mp_nucl = {}
for rec in SeqIO.parse(open("data/sinv_blast/mp_loci.fa"),"fasta"):
	if rec.id in mp_genes:
		# reverse complement genes on opposite strand
		if rec.id in mp_strand and mp_strand[rec.id] == '-':
			rec.seq = rec.seq.reverse_complement()
		mp_nucl[rec.id] = rec

si_nucl = {}
for rec in SeqIO.parse(open("data/sinv_blast/sinv_nucl.fa"),"fasta"):
	if str(rec.seq.translate())[:-1].find("*") == -1: #ignore pseudogenes
		if rec.id in genes:
			si_nucl[rec.id] = rec

# write nucleotide and protein sequence pares for homologous sequences
for si in genes:
	mp = genes[si]
	if mp in mp_nucl and si in si_nucl and mp in mp_prot:
		outfile_nucl =open("data/sinv_blast/nucl/"+mp+".fa","w")
		outfile_prots =open("data/sinv_blast/prots/"+mp+".fa","w")
		outfile_nucl.write(mp_nucl[mp].format("fasta"))
		outfile_nucl.write(si_nucl[si].format("fasta"))
		si_nucl[si].seq = si_nucl[si].seq.translate()
		outfile_prots.write(mp_prot[mp].format("fasta"))
		outfile_prots.write(si_nucl[si].format("fasta"))
		outfile_nucl.close()
		outfile_prots.close()
