#!/bin/bash
#$ -q shortP
#$ -j y
#$ -cwd
#$ -N bt
#$ -pe smp 12
#$ -l h_vmem=4G
#$ -l virtual_free=4G

. $HOME/.bashrc

bowtie2 --sam-rg ID:mp --sam-rg LB:Truseq --sam-rg SM:mp --sam-rg PL:ILLUMINA  -p 8 -x ref/Mp -1 data/assembly/reads/trim1.fastq.gz -2 data/assembly/reads/trim2.fastq.gz |samtools view -Su - | novosort --ram 24G -t 4 -i -o data/remapped.bam -t /genefs/MikheyevU/temp -

