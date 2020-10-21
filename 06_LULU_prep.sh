#!/bin/bash
#SBATCH -J LULU_prep
#SBATCH -o out_err_files/LULU_prep_%A.out
#SBATCH -e out_err_files/LULU_prep_%A.err
#SBATCH --nodes=1
#SBATCH -t 4-00:00:00
#SBATCH -p defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

module use /GWSPH/groups/cbi/Apps/_modulefiles
module load vsearch/2.14.1
module load R/3.5.3

#vsearch --fastx_filter uniques.fasta --relabel_sha1 --fastaout DADA2_uniq_sha1.fasta
#echo "OTUId" > DADA2_sha1_labels.txt
#grep ">" DADA2_uniq_sha1.fasta | sed 's/>//' >> DADA2_sha1_labels.txt

# generate OTU table from seqtab_final.RDS
Rscript LULU_prep.R
exit_status=$?
cat DADA2_extracted_samples/*.fas > DADA2_extracted_samples/concat.fasta

exit $exit_status
# then run vsearc_cluster.sh and gen_matchlist.sh
