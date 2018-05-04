#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N matchlist
#$ -o out_err_files/matchlist_$JOB_ID.log
#$ -q sThC.q
#$ -pe mthread 4
#$ -l mres=1G,h_data=1G,h_vmem=1G


module use /data/genomics/nguyenbn/modulefiles
module load vsearch/2.8.0
#module load bioinformatics/blast/2.6.0 

CENTROIDS_FILE="vsearch_size_centroid97.fasta"

# Make BLAST database with centroids of clusters
#echo "makeblastdb -in $CENTROIDS_FILE -parse_seqids -dbtype nucl"

# BLAST the centroids against themselves

#echo "blastn -db $CENTROIDS_FILE -outfmt '6 qseqid sseqid pident' -out matchlist.txt -qcov_hsp_perc 80 -perc_identity 84 -query $CENTROIDS_FILE -num_threads $NSLOTS"


# Use VSEARCH to generate matchlist.txt
vsearch --usearch_global $CENTROIDS_FILE --db $CENTROIDS_FILE --self --id .84 --iddef 1 --userout matchlist.txt --userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10 --threads $NSLOTS

