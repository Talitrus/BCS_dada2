#!/bin/bash
#SBATCH -J R_rename
#SBATCH -o out_err_files/R_%A.out
#SBATCH -e out_err_files/R_%A.err
#SBATCH --nodes=1
#SBATCH -t 01:00:00
#SBATCH -p tiny,short,defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

module load R/3.5.3
Rscript rename_postQC.R
