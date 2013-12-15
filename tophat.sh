#!/bin/bash
#$ -q long
#$ -j y
#$ -cwd
#$ -N tophat
. $HOME/.bashrc
a=(*1.fastq.gz)
#SGE_TASK_ID=1
b=$(basename ${a["SGE_TASK_ID"-1]} "_R1.fastq.gz")

tophat2 -z pigz -p 12 -o $b ../Mp "$b"_R1.fastq.gz "$b"_R2.fastq.gz