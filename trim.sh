#!/bin/bash
#$ -q long
#$ -j y
#$ -cwd
#$ -N trim
. $HOME/.bashrc
a=(split/1??)
b=(split/2??)
f=${a["SGE_TASK_ID"-1]}
r=${b["SGE_TASK_ID"-1]} 

condetri_v2.2.pl -sc=33 -fastq1=$f -fastq2=$r
#condetri_v2.2.pl -sc=33 -fastq1=R1.fastq -fastq2=R2.fastq

