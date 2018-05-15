#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N R_phyloseq
#$ -o out_err_files/R_phyloseq_$JOB_ID.log
#$ -q mThM.q
#$ -pe mthread 1
#$ -l mres=16G,h_data=16G,h_vmem=16G,himem

module load tools/R
Rscript BCS_phyloseq.R $NSLOTS
