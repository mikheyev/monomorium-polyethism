#takes a nucleotide file
#makes a codon alignment, assuming that the nucleotide sequence is already correctly aligned,
#then runs paml in pairwise mode to estimate dn/ds
from Bio.Phylo.PAML import codeml
from Bio import SeqIO
import os,sys,pdb

nucl_file = sys.argv[1]
fname = sys.argv[1].split("/")[-1] 
prot_file = "data/sinv_blast/prots/" + fname
outdir = "data/sinv_blast/alignments/"
os.system("muscle -quiet -in "+prot_file + " -out " + outdir + "prot/" + fname)
os.system("pal2nal.pl " + outdir + "prot/" + fname + " " + nucl_file + " -output paml > " + outdir + "nucl/" + fname.replace(".fa",".phy"))
#PAML analysis
cml = codeml.Codeml()
cml.alignment =  outdir + "nucl/" + fname.replace(".fa",".phy")
cml.working_dir = "./data/sinv_blast/paml/scratch_" + fname.split(".")[0]
seqnames = SeqIO.to_dict(SeqIO.parse(open(nucl_file),"fasta")).keys()
outfile = open("./data/sinv_blast/paml/" + fname.replace(".fa",".tre"),"w")
outfile.write("(" + seqnames[0][:10] + "," + seqnames[1] + ");")
outfile.close()
cml.tree = "./data/sinv_blast/paml/" + fname.replace(".fa",".tre")
cml.out_file = "./data/sinv_blast/paml/" + fname.replace(".fa",".txt")
cml.working_dir = "./data/sinv_blast/paml/scratch_" + fname.split(".")[0]
cml.set_options(seqtype=1,
        verbose=1,
        noisy=0,
        model=1,
        runmode=-2,
        Mgene=0,
        NSsites=[0],
        CodonFreq=2,
        cleandata=0)
cml.run(verbose=True)

# parse cml.out_file
for line in open(cml.out_file):
    if line.find("dN/dS=") > -1:
        line = line.split()
        print "%s\t%s" % (fname.split(".")[0], line[line.index("dN/dS=")+1])
