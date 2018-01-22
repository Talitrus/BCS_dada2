#!/bin/bash
#SBATCH -J R
#SBATCH -o out_err_files/R_%A.out
#SBATCH -e out_err_files/R_%A.err
#SBATCH --nodes=1
#SBATCH -t 02:00:00
#SBATCH -p defq,short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

module load R/3.4.2
Rscript rename_postQC.R
