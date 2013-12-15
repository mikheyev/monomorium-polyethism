#!/bin/bash
#$ -q long
#$ -j y
#$ -cwd
#$ -N trim
. $HOME/.bashrc

#cat split/*_trim1.fastq | awk 'BEGIN {RS="@VADER"; ORS="@VADER"; FS="\n";OFS="\n"} $1~"1:N:0:AGTCAA" {sub(/ 1:N:0:AGTCAA/,"/1"); print}'  > trim1.fastq 

#cat split/*_trim2.fastq | awk 'BEGIN {RS="@VADER"; ORS="@VADER"; FS="\n";OFS="\n"} $1~"2:N:0:AGTCAA" {sub(/ 2:N:0:AGTCAA/,"/2"); print}'  > trim2.fastq 

cat split/*_unpaired.fastq | awk 'BEGIN {RS="@VADER"; ORS="@VADER"; FS="\n";OFS="\n"} $1~"AGTCAA" '  > trimu.fastq 