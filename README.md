# Behavioral and transcriptomic analysis of age-based polyethism in the pharaoh's ant (*Monomorium pharaonis*)

## Reference based transcriptome assembly and gene expression counts
After filtering adaptor (*rename.sh*) and low quality sequence (*trim.sh*), concatenate results (*concat.sh*) and assemble genome reference using ABYSS (*abyss.sh*), using a range of k-mers from 51 to 91. k=69 produced the longest N50, and was chosen for future work.

	n		n:500	n:N50	min	N80		N50		N20		max		sum
	216388	42794	4628	500	5605	16262	36932	246554	283e6	contigs
	207026	36001	3980	500	6743	18960	43332	246554	283.7e6	scaffolds

Map transcriptional data to the reference genome using tophat, and infer loci and isoforms using cufflinks (*tophat.sh* and *cuff.sh*). After merging all the libraries using *cuffmerge*, extract reference transcriptome using cufflink's *gffread* utility.

Now use *rsem.sh* to map reads to extracted loci, providing a genes to isoforms table based on cufflinks output. This produces files with counts, which are then collected (*collect_counts.py* and *collect_fpkm.py*) entered into a SQL database.

Mapping stats for rsem

	id  	mapped		percent	library
	MP11	4802838 	33.73%	age6
	MP3 	4670833 	34.87%	age3
	MP5 	5785570 	38.60%	age12
	MP9 	7269589 	52.91%	nurse
	MP17	7767423 	59.50%	nurse
	MP22	11480354	67.10%	groom
	MP21	9769906 	69.22%	troph
	MP23	10854923	73.04%	age18
	MP1 	10664516	77.03%	age0
	MP16	11809890	82.97%	groom
	MP20	14648253	83.64%	age9
	MP19	14522217	84.37%	prot
	MP7 	11782737	85.74%	troph
	MP14	10697715	86.00%	age18
	MP8 	12053556	86.29%	age0
	MP2 	11246559	86.58%	age3
	MP10	10379494	86.64%	age6
	MP12	13243001	86.72%	age15
	MP4 	12984328	86.74%	carb
	MP13	11114399	86.91%	age15
	MP24	16818668	86.98%	prot
	MP18	11041174	87.14%	carb
	MP6 	12885559	87.29%	age12
	MP15	12227812	87.50%	age9

**Average: 74.9% +/- 17.5%**

## Transcriptome annotation and evolutionary rates analysis

Blastx assembled transcriptome vs NR database with evalue 0.00001. Run blast2go on the blast results to obtain GO term annotations. These are collected (*collect_go.py*) and uploaded to a SQL database


First, we need to get just one locus for every isoform. For this, we just choose the longest transcript.

 	python choose_isoform.py > ref/longest_isoforms.fa

This produces 22385 loci.

### Evolutionary rates vs fire ants (*Solenopsis invicta*)
Using *longest_isoforms.fa* conduct reciprocal best blast vs the *S. invicta* OGSv2.2.3 using blastx and and tblastn with evalue evalue 1e-10. These blast results (which were split into groups of 100 sequences for speed), can also be used to predict protein codeing genes using [OrfPredictor](http://www.ncbi.nlm.nih.gov/pubmed/15980561):

	#predict proteins using blast results
	for i in *.fa ; do perl /apps/MikheyevU/sasha/ORFPredictor/OrfPredictor_web3.pl $i `basename $i .fa`.xml 1 both test@test.com 1e-10 `basename $i .fa`_prot.fa test; done

Parse blast data using *combine_homologs.py* to extract reciprocal best blast hits, and use *combine_sequences.py* to create fasta files containing the *S. invicta* and *M. pharaonis* homologs.

Run *align.py* (parallelized by *align.sh*) to create a codon alignment and compute dN/dS ratios using PAML.

A similar approach was used to find honeybee (*Apis mellifera* homologs), except OGSv1.0 was used and the evalue was relaxed to 1e-5.

## Statistics and plotting

*monomorium.R* contains all of the statistical analysis.

## genome annotation

### Round 1

Used training set from Wasmannia and Vollenhovia for first pass. After the run was complete, extracted every prediction that had good transcriptional support

	maker2zff  -c 1 -e 1 -o 1

### Round 2

#### SNAPP training

	fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
	fathom genome.ann genome.dna -validate > validate.log 2>&1
	fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
	fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
	# 3328 genes
	mkdir params; cd params
	forge ../export.ann ../export.dna > ../forge.log 2>&1
	cd ..; hmm-assembler.pl vol params/ > mp.hmm

#### Augustus training

	perl zff2augustus_gbk.pl > genes.gb
	etraining --species=generic --stopCodonExcludedFromCDS=false genes.gb 2> train.err
	cat train.err | perl -pe 's/.*in sequence (\S+): .*/$1/' > badgenes.lst
	# no bad genes
#	filterGenes.pl badgenes.lst round1.gb > genes.gb
	randomSplit.pl genes.gb 200
	autoAug.pl --species=ph --genome ../../../ref/Mp.fa --trainingset genes.gb   # 2673987
	optimize_augustus.pl --cpus=12  --UTR=on --species=mph genes.gb.train
	augustus --species=vem genes.gb.test >train_out.txt
	grep -A 22 Evaluation  train_out.txt
	

	

