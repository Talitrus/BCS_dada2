#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N LULU
#$ -o out_err_files/LULU_$JOB_ID.log
#$ -q mThM.q
#$ -l mres=16G,h_data=16G,h_vmem=16G,himem

module use /data/genomics/nguyenbn/modulefiles
module load tools/R

# Takes ~1.5 days
# generate curated OTU table
Rscript LULU.R
