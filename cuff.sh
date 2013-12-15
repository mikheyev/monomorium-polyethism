#!/bin/bash
#$ -q long
#$ -j y
#$ -cwd
#$ -N cufflinks
. $HOME/.bashrc
a=(*/a*bam)
#SGE_TASK_ID=1
b=$(dirname ${a["SGE_TASK_ID"-1]})

cufflinks -p 12 -o $b ${a["SGE_TASK_ID"-1]}

