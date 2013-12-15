#!/bin/bash
#$ -N paml
#$ -q short
#$ -l h_vmem=4G
#$ -l virtual_free=4G
#$ -j y
##$ -m ea
#$ -cwd
. $HOME/.bashrc
a=(data/sinv_blast/nucl/XLO*fa) # 8282
python align.py ${a["SGE_TASK_ID"-1]}
