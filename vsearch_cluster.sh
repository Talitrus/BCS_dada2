#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N vsearch_clust
#$ -o out_err_files/vs_clust_$JOB_ID.log
#$ -q sThC.q
#$ -pe mthread 4
#$ -l mres=3G,h_data=3G,h_vmem=3G

module use /data/genomics/nguyenbn/modulefiles
module load vsearch/2.8.0

#UC_OUTFILE="vs_us_out_size97_UC"
OUT_TAB_FILE="cluster_hits.tsv"

# Variables
TMP_concat_derep=$(mktemp)
TMP_FASTA1=$(mktemp)
TMP_FASTA2=$(mktemp)
#UNSORTED_CENTROIDS=$(mktemp)
FINAL_CENTROIDS="vsearch_size_centroid97.fasta"
#LOW_ABUND_THRES=1


# Dereplicate (vsearch)
vsearch --threads $NSLOTS --derep_fulllength DADA2_extracted_samples/concat.fasta --sizein --sizeout --fasta_width 0 --minuniquesize 1 --output $TMP_concat_derep

# Remove extra semicolon added by VSEARCH
#sed -i -E "s/;barcodelabel=.+;;/;/g" $TMP_concat_derep
sed -i -E "s/;.+;/;/g" $TMP_concat_derep

# Sort by size
vsearch --threads $NSLOTS --fasta_width 0 --sortbysize $TMP_concat_derep --sizein --sizeout --output $TMP_FASTA1

# VSEARCH chimera checking
vsearch --threads $NSLOTS --uchime_denovo $TMP_FASTA1 --sizein --sizeout --dbmask none --qmask none --nonchimeras $TMP_FASTA2

# VSEARCH clustering
vsearch --threads $NSLOTS --cluster_size $TMP_FASTA2 --id 0.97 --sizein --centroids $FINAL_CENTROIDS --qmask none --dbmask none

# VSEARCH filtering and sorting
#vsearch --threads $NSLOTS --fasta_width 0 --sortbysize $UNSORTED_CENTROIDS --minsize $LOW_ABUND_THRES --output $FINAL_CENTROIDS --sizein

# Shouldn't be necessary to remove size annotation in sequence headers, but check for them anyway.
sed -i -E "s/;size=.*//g" $FINAL_CENTROIDS
# Make a tab-separated file containing two columns: hit sequence, cluster centroid sequence. NOT DOING THIS RIGHT NOW. NOT USED
#grep -E "^H" $UC_OUTFILE | awk -F'\t' '{ print $9 "\t" $10 }' | sed -E 's/;size=[0-9]+;//g' >$OUT_TAB_FILE

# Map raw reads against OTU representatives concat.fasta is from LULU_prep.sh, which should be run before this.
vsearch --usearch_global DADA2_extracted_samples/concat.fasta --db $FINAL_CENTROIDS --id 0.97 --maxaccepts 0 --dbmask none --qmask none --uc FINAL_FOR_LULU.uc --threads $NSLOTS


# Then convert UC format to OTU table format with uc2otutab.py
python d5_py/uc2otutab.py FINAL_FOR_LULU.uc > FINAL_FOR_LULU_otutab.txt

# remove tmp files used
rm -f $TMP_concat_derep $TMP_FASTA1 $TMP_FASTA2 $UNSORTED_CENTROIDS
