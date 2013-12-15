#!/bin/bash
#$ -q long
#$ -j y
#$ -cwd
#$ -N trim
. $HOME/.bashrc

#cat reads/*R2.fastq |cutadapt -a AGATCGGAAGAGC - | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | awk '{OFS=""} NF==2 {$2="/2";print;next} {print}' > R2.fastq 

cat reads/*R1.fastq |cutadapt -a AGATCGGAAGAGC - | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" | awk '{OFS=""} NF==2 {$2="/1";print;next} {print}' > R1.fastq 
