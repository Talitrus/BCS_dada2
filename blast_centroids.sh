#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N BLAST
#$ -o out_err_files/BLAST_$JOB_ID.log
#$ -q sThC.q
#$ -pe mthread 4
#$ -l mres=2G,h_data=2G,h_vmem=2G

module load bioinformatics/blast/2.6.0 

CENTROIDS_FILE="vsearch_size_centroid97.fasta"

# Make BLAST database with centroids of clusters
makeblastdb -in $CENTROIDS_FILE -parse_seqids -dbtype nucl

# BLAST the centroids against themselves
blastn -db $CENTROIDS_FILE -outfmt '6 qseqid sseqid pident' -out matchlist.txt -qcov_hsp_perc 80 -perc_identity 84 -query $CENTROIDS_FILE -num_threads 4

