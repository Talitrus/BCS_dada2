#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N LULU
#$ -o out_err_files/LULU_$JOB_ID.log
#$ -q mThC.q
#$ -l mres=5G,h_data=5G,h_vmem=5G

module use /data/genomics/nguyenbn/modulefiles
module load tools/R

# Takes ~1.5 days
# generate curated OTU table
Rscript LULU.R
