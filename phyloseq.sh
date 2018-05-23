#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N R_phyloseq
#$ -o out_err_files/R_phyloseq_$JOB_ID.log
#$ -q mThC.q
#$ -pe mthread 6
#$ -l mres=4G,h_data=4G,h_vmem=4G

module load tools/R
Rscript BCS_phyloseq.R $NSLOTS
