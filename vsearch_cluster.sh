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

UC_OUTFILE="vs_us_out_size97_UC"
OUT_TAB_FILE="cluster_hits.tsv"



vsearch --threads $NSLOTS --cluster_size DADA2_uniq_sha1.fasta --id 0.97 --sizein --centroids vsearch_size_centroid97.fasta --qmask none --dbmask none --uc $UC_OUTFILE

# Make a tab-separated file containing two columns: hit sequence, cluster centroid sequence.
grep -E "^H" $UC_OUTFILE | awk -F'\t' '{ print $9 "\t" $10 }' | sed -E 's/;size=[0-9]+;//g' >$OUT_TAB_FILE

# Map raw reads against OTU representatives concat.fasta is from LULU_prep.sh, which should be run before this.
vsearch --usearch_global DADA2_extracted_samples/concat.fasta --db vsearch_size_centroid97.fasta --id 0.97 --maxaccepts 0 --dbmask none --qmask none --uc FINAL_FOR_LULU.uc --threads $NSLOTS


# Then convert UC format to OTU table format with uc2otutab.py
python d5_py/uc2otutab.py FINAL_FOR_LULU.uc > FINAL_FOR_LULU_otutab.txt
