#!/bin/bash
#$ -q short
#$ -j y
#$ -cwd
#$ -l h_vmem=10G
#$ -l virtual_free=10G
#$ -N kg
. $HOME/.bashrc
export TEMPDIR=/genefs/MikheyevU/temp
export TMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

# load tombo modules

module load R/3.1.0 

#kmergenie --diploid -l 45 -k 95 -o muuler ../data/attamyces/reads/strain1/mueller_merged.tmp.fq
#set +o posix
kmergenie -o data/kmer data/assembly/reads/trim1.fastq.gz

#kmergenie --diploid -o strain1 <(zcat ../data/attamyces/reads/strain1/*fq.gz)
