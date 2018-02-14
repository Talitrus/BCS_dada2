#!/bin/bash
#SBATCH -J BLCA
#SBATCH -o out_err_files/BLCA_%A.out
#SBATCH -e out_err_files/BLCA_%A.err
#SBATCH --nodes=1
#SBATCH -t 2-00:00:00
#SBATCH -p defq,short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

module load BLCA
python /groups/cbi/shared/apps/BLCA/BLCA.git/2.blca_main.py -i uniques_trim_copy.fasta -r /groups/cbi/bryan/Midori_COI/midori_blca_taxonomy_tab.txt -q /groups/cbi/bryan/Midori_COI/midori_blca_blast -b 70 -c 0.75
