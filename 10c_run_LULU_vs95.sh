#!/bin/bash
#SBATCH -J LULU
#SBATCH -o out_err_files/LULU_%A.out
#SBATCH -e out_err_files/LULU_%A.err
#SBATCH --nodes=1
#SBATCH -t 5-00:00:00
#SBATCH -p defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

module load R/3.5.3

# Takes ~1.5 days
# generate curated OTU table
Rscript LULU_vs95.R
