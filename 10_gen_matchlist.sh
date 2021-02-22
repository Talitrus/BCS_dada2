#!/bin/bash
#SBATCH -J gen_matchlist
#SBATCH -o out_err_files/matchlist_%A.out
#SBATCH -e out_err_files/matchlist_%A.err
#SBATCH --nodes=1
#SBATCH -t 8:00:00
#SBATCH -p tiny,short,defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

module load vsearch/2.14.1
module load blast+

CENTROIDS_FILE="vsearch_size_centroid97.fasta"
NSLOTS=$(nproc)
# Make BLAST database with centroids of clusters
makeblastdb -in $CENTROIDS_FILE -parse_seqids -dbtype nucl

# BLAST the centroids against themselves

blastn -db $CENTROIDS_FILE -outfmt '6 qseqid sseqid pident' -out matchlist.txt -qcov_hsp_perc 80 -perc_identity 84 -query $CENTROIDS_FILE -num_threads $NSLOTS


# Use VSEARCH to generate matchlist.txt
#vsearch --usearch_global $CENTROIDS_FILE --db $CENTROIDS_FILE --self --id .84 --iddef 1 --userout matchlist.txt --userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10 --threads $NSLOTS
