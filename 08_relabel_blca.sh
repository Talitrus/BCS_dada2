#!/bin/bash
#SBATCH -p tiny,short,defq
#SBATCH -t 00:15:00
#SBATCH -J relabel
#SBATCH --nodes=1
#SBATCH -o out_err_files/relabel_%A.out
#SBATCH -e out_err_files/relabel_%A.err

#Must be run AFTER vsearch clustering
grep ">" uniques_sha1.fasta | sed -E "s/>//" > unique_sha1_labels.txt
module load R/3.5.3
Rscript relabel_blca.R
