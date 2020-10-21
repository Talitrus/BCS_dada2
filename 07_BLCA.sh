#!/bin/bash
#SBATCH -J BLCA
#SBATCH -o out_err_files/BLCA_%A.out
#SBATCH -e out_err_files/BLCA_%A.err
#SBATCH --nodes=1
#SBATCH -t 1-00:00:00
#SBATCH -p defq,short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

module load BLCA

sed 's/;.*;//g' uniques.fasta > uniques_trim.fasta
python /GWSPH/groups/cbi/Apps/BLCA/git/2.blca_main.py -i uniques_trim.fasta -r /lustre/groups/cbi/Users/bnguyen/Midori_COI/GB239/midori_blca/midori_kingdom_blca_taxonomy.txt -q /lustre/groups/cbi/Users/bnguyen/Midori_COI/GB239/midori_blca/midori_blca_blast -b 80
