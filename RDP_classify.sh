#!/bin/bash
#SBATCH -J J_RDP
#SBATCH -o out_err_files/RDP_%A.out
#SBATCH -e out_err_files/RDP_%A.err
#SBATCH --nodes=1
#SBATCH -t 14-00:00:00
#SBATCH -p defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu

module load RDP
java -Xmx64g -jar /groups/cbi/shared/apps/RDPTools/classifier.jar classify -t ../../Midori_COI/trained/rRNAClassifier.properties -o BCS_RDP_output.txt -h BCS_RDP_hier.txt uniques.fasta
