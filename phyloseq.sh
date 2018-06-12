#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N R_phyloseq
#$ -o out_err_files/R_phyloseq_$JOB_ID.log
#$ -q sThM.q
#$ -pe mthread 6
#$ -l mres=9G,h_data=9G,h_vmem=9G,himem

module load tools/R
Rscript BCS_phyloseq.R $NSLOTS
