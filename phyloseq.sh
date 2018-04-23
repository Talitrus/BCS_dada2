#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N R_phyloseq
#$ -o out_err_files/R_phyloseq.log
#$ -q mThM.q
#$ -pe mthread 8
#$ -l mres=6G,h_data=6G,h_vmem=6G,himem

module load tools/R
Rscript BCS_phyloseq.R $NSLOTS
