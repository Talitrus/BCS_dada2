#!/bin/bash
#SBATCH -J DADA2_filt
#SBATCH -o out_err_files/R_DADA2_filt_%A_%a.out
#SBATCH -e out_err_files/R_DADA2_filt_%A_%a.err
#SBATCH --nodes=1
#SBATCH -t 2-00:00:00
#SBATCH -p defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

module load R/3.5.3
Rscript dada2_analysis.R
