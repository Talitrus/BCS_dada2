#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N vsearch_clust
#$ -o out_err_files/vs_clust_$JOB_ID.log
#$ -q sThC.q

module use /data/genomics/nguyenbn/modulefiles
module load vsearch/2.8.0
module load tools/R

vsearch --fastx_filter uniques.fasta --relabel_sha1 --fastaout DADA2_uniq_sha1.fasta
echo "OTUId" > DADA2_sha1_labels
grep ">" DADA2_uniq_sha1.fasta | sed 's/>//' >> DADA2_sha1_labels.txt

# generate OTU table from seqtab_final.RDS
RScript LULU_prep.R $NSLOTS

cat DADA2_extracted_samples/*.fas > DADA2_extracted_samples/concat.fasta

# then run vsearc_cluster.sh and gen_matchlist.sh