from Bio import SeqIO

# Choose longest isoform form cufflinks assembly for comparative genomics

transcripts = {}

for rec in SeqIO.parse(open("ref/transcripts.fa"),"fasta"):
	name = rec.description.split("=")[1]
	if name in transcripts and transcripts[name][0] >= len(rec.seq):
		continue
	transcripts[name] = [len(rec.seq),rec.id]

for rec in SeqIO.parse(open("ref/transcripts.fa"),"fasta"):
	name = rec.description.split("=")[1]
	if transcripts[name][1] == rec.id:
		rec.description = rec.id
		rec.id = name
		print rec.format("fasta"),