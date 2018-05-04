#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N vsearch_clust
#$ -o out_err_files/vs_clust_$JOB_ID.log
#$ -q sThC.q
#$ -pe mthread 4
#$ -l mres=2G,h_data=2G,h_vmem=2G

module use /data/genomics/nguyenbn/modulefiles
module load vsearch/2.8.0

UC_OUTFILE="vs_us_out_size97_OTUTABLE"
OUT_TAB_FILE="cluster_hits.tsv"



vsearch --threads $NSLOTS --cluster_size uniques_seqheaders.fasta --id 0.97 --sizein --centroids vsearch_size_centroid97.fasta --qmask none --dbmask none --uc $UC_OUTFILE

# Make a tab-separated file containing two columns: hit sequence, cluster centroid sequence.
grep -E "^H" $UC_OUTFILE | awk -F'\t' '{ print $9 "\t" $10 }' | sed -E 's/;size=[0-9]+;//g' >$OUT_TAB_FILE
