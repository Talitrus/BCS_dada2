#!/bin/bash
#SBATCH -J R
#SBATCH -o out_err_files/R_%A.out
#SBATCH -e out_err_files/R_%A.err
#SBATCH --nodes=1
#SBATCH -t 14-00:00:00
#SBATCH -p defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

module load R/3.4.2
Rscript phy_copy.R
