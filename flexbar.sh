#!/bin/bash
#SBATCH -J flexbar
#SBATCH -o out_err_files/flexbar_%A_%a.out
#SBATCH -e out_err_files/flexbar_%A_%a.err
# assign array, then below = how many nodes you want.
#SBATCH --array=1-8
#SBATCH --nodes=1
# time stamp for the how long you expect the longest job in the array to take 
# each will run with that same time stamp specified)
#SBATCH -t 03:00:00
#SBATCH -p defq,short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

#We're going to just trim for adapters here. I haven't figured out an easy way to demultiplex Matt's pseudo dual-index primers yet. We will trim the primers (read 1) and index 2 + primers (read 2) off inside DADA2.
module load flexbar

cd ..
libname=$(sed -n "$SLURM_ARRAY_TASK_ID"p libraries.txt)


# move to the directory where the data files are located
cd $libname
list_name="file_roots.txt"
rm *flex* #remove any old runs

while IFS='' read -r line || [[ -n "$line" ]]; do
	#cmd_str="flexbar -r ${line}1_001.fastq -p ${line}2_001.fastq -b ../ML_barcodes.fasta --barcode-unassigned -t ${line}_flex -n $(nproc) --barcode-trim-end LTAIL --adapter-trim-end ANY -a ../illumina_adapters.fasta"
	cmd_str="flexbar -r ${line}1_001.fastq -p ${line}2_001.fastq -t ${line}_flex -n $(nproc) --adapter-trim-end ANY -a ../truseq_adapters.fasta -ao 7 -b ../ML_barcodes.fasta --barcode-trim-end LTAIL"
	echo "$cmd_str"
	$cmd_str
done < "$list_name"
