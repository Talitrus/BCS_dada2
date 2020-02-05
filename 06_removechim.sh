#!/bin/bash
#SBATCH -J DADA2_chimera
#SBATCH -o out_err_files/R_DADA2_chimera_%A.out
#SBATCH -e out_err_files/R_DADA2_chimera_%A.err
#SBATCH --nodes=1
#SBATCH -t 08:00:00
#SBATCH -p tiny,short,defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

module use /groups/cbi/modulefiles
module load R/3.5.3
module load vsearch
Rscript dada2_removechim.R
vsearch --fastx_filter uniques.fasta --relabel_sha1 --fastaout uniques_sha1.fasta
