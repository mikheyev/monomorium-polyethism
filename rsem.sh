#!/bin/bash
#$ -q long
#$ -j y
#$ -cwd
#$ -N rsem
#$ -l h_vmem=4G
#$ -l virtual_free=4G
. $HOME/.bashrc
export TEMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp
a=(../*1.fastq)
b=(../*2.fastq)
base=$(basename ${a["SGE_TASK_ID"-1]} "_R1.fastq")
f=${a["SGE_TASK_ID"-1]}
r=${b["SGE_TASK_ID"-1]} 

rsem-calculate-expression -p 4 --paired-end $f $r mp  $base   