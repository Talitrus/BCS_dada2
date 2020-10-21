#!/bin/bash
#SBATCH -J flexbar
#SBATCH -o out_err_files/flexbar_%A_%a.out
#SBATCH -e out_err_files/flexbar_%A_%a.err
# assign array, then below = how many nodes you want.
#SBATCH --array=18-19
#SBATCH --nodes=1
# time stamp for the how long you expect the longest job in the array to take 
# each will run with that same time stamp specified)
#SBATCH -t 03:00:00
#SBATCH -p tiny,short,defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

#module load flexbar
module load miniconda/miniconda3
source activate flexbar

cd ..
libname=$(sed -n "$SLURM_ARRAY_TASK_ID"p libraries.txt)


# move to the directory where the data files are located
cd $libname
list_name="file_roots.txt"
rm *flex* #remove any old runs
gunzip *.fastq.gz

while IFS='' read -r line || [[ -n "$line" ]]; do
	cmd_str="flexbar -r ${line}1_001.fastq -p ${line}2_001.fastq -t ${line}_flex -n $(nproc) --adapter-trim-end ANY -a ../truseq_adapters.fasta -ao 7 -b ../ML_barcodes.fasta --barcode-trim-end LTAIL"
	echo "$cmd_str"
	$cmd_str
done < "$list_name"

# Delete any files under 100k in size to prevent errors in flexbar2.sh
mkdir -p flexbar_filtered
find . -name "*.fastq" -size -100k -exec mv '{}' flexbar_filtered/ \; #move files instead of deleting them
#find . -name "*.fastq" -size -100k -delete
