#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N R_phyloseq
#$ -o out_err_files/R_phyloseq_$JOB_ID.log
#$ -q mThC.q
#$ -pe mthread 3
#$ -l mres=3G,h_data=3G,h_vmem=3G

module load tools/R
Rscript BCS_phyloseq.R $NSLOTS
