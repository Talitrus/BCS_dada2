#!/bin/bash
#SBATCH -J flexbar
#SBATCH -o out_err_files/flexbar_%A_%a.out
#SBATCH -e out_err_files/flexbar_%A_%a.err
# assign array, then below = how many nodes you want.
#SBATCH --array=1-15
#SBATCH --nodes=1
# time stamp for the how long you expect the longest job in the array to take 
# each will run with that same time stamp specified)
#SBATCH -t 03:00:00
#SBATCH -p defq,short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

module load flexbar

cd ..
libname=$(sed -n "$SLURM_ARRAY_TASK_ID"p libraries.txt)


# move to the directory where the data files are located
cd $libname
rm *final*
ls *barcode*.fastq | grep -v 'unassigned' | sed -E 's/[12].fastq//' | sort -n | uniq > file_roots2.txt
list_name="file_roots2.txt"

while IFS='' read -r line || [[ -n "$line" ]]; do
	cmd_str="flexbar -r ${line}2.fastq -p ${line}1.fastq -t ${line}final -n $(nproc) -b ../ML_barcodes.fasta --barcode-trim-end LTAIL"
	echo "$cmd_str"
	$cmd_str
done < "$list_name"
FLEXBAR_EXIT_STATUS=$?


find . -type f -name '*flex_barcode_[0-9]_[0-9].fastq' -delete #remove first stage of flexbar run
find . -type f -name '*flex_barcode_[0-9][0-9]_[0-9].fastq' -delete
find . -name "*.fastq" -size -450k -delete #delete very small files from second stage of flexbar
mkdir logs
mv *.log logs/

mkdir -p unassigned
mv *unassigned*.fastq unassigned/
mv Undetermined*.fastq unassigned/

#Override exit status. Not the most elegant, but it should be fine.
exit $FLEXBAR_EXIT_STATUS
