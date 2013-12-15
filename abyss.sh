#!/bin/bash
#$ -N mp_abyss
#$ -q long
#$ -j y
##$ -m ea
#$ -cwd
##$ -l mf=40G
##$ -l h_vmem=40G
. $HOME/.bashrc
rm -r k$SGE_TASK_ID
mkdir k$SGE_TASK_ID
dir=$(pwd)
cd k$SGE_TASK_ID
abyss-pe name=mc j=10 k=$SGE_TASK_ID in='/genefs/MikheyevU/sasha/monomorium/trim1.fastq /genefs/MikheyevU/sasha/monomorium/trim2.fastq' se='/genefs/MikheyevU/sasha/monomorium/trimu.fastq'
