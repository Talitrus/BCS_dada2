#!/bin/bash
#SBATCH -J flexbar
#SBATCH -o out_err_files/flexbar_%A_%a.out
#SBATCH -e out_err_files/flexbar_%A_%a.err
# assign array, then below = how many nodes you want.
#SBATCH --array=1-4
#SBATCH --nodes=1
# time stamp for the how long you expect the longest job in the array to take 
# each will run with that same time stamp specified)
#SBATCH -t 02:00:00
#SBATCH -p tiny,short,defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

module load flexbar

cd ..
libname=$(sed -n "$SLURM_ARRAY_TASK_ID"p libraries.txt)
ROOT_DIR="/groups/cbi/Users/bnguyen/bocas/BCS_18S"

# move to the directory where the data files are located
cd $libname
list_name="file_roots.txt"
rm *flex* #remove any old runs

while IFS='' read -r line || [[ -n "$line" ]]; do
	cmd_str="flexbar -r ${line}1_001.fastq* -p ${line}2_001.fastq* -t ${line}_flex -n $(nproc) -b $ROOT_DIR/scripts/helper_files/ML_barcodes.fasta --barcode-trim-end LTAIL" #no gzip used in output here because the files will get deleted shortly anyway.
	echo "$cmd_str"
	$cmd_str
done < "$list_name"

# Delete any files under 100k in size to prevent errors in flexbar2.sh
find . -name "*.fastq" -size -100k -delete
