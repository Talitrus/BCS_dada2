#!/bin/bash
#SBATCH -J swarm
#SBATCH -o out_err_files/swarm_%A.out
#SBATCH -e out_err_files/swarm_%A.err
#SBATCH --nodes=1
#SBATCH -t 01:00:00
#SBATCH -p tiny,short,defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnguyen@gwu.edu


module use /GWSPH/groups/cbi/Apps/_modulefiles
module load swarm
module load vsearch

#TEMP FILES
TMP_concat_derep=$(mktemp)
UCHIME="swarm_13.uc"
SWARM_REPS=$(mktemp)
SORTED_REPS="swarm_reps_sorted.fasta"
SWARM_STATS="swarm_13.stats"
SWARM_STRUCT="swarm_13.struct"
EXTRACT_DIR="DADA2_extracted_samples"
SWARM_EXTRACT_DIR="for_swarm/DADA2_extracted_samples"
SWARMS="uniques_13.swarms"
OTU_TABLE="swarm_13.otutable"

#change names to SHA1 with size labels
vsearch --fastx_filter DADA2_extracted_samples/concat.fasta --relabel_sha1 --sizein --sizeout --fastaout concat_sha1_size.fasta

# Dereplicate (vsearch)
vsearch --derep_fulllength concat_sha1_size.fasta --threads $(nproc) --sizein --sizeout --fasta_width 0 --output $TMP_concat_derep
#sed 's/;//g' uniques.fasta | sed 's/size=/_/g' > uniques_swarm.fasta
swarm -t $(nproc) -d 13 -z -w $SWARM_REPS -i $SWARM_STRUCT -s $SWARM_STATS -o $SWARMS $TMP_concat_derep

#Sort representatives
vsearch --sortbysize $SWARM_REPS --fasta_width 0 --output $SORTED_REPS -threads $(nproc)

# Chimera checking

vsearch --uchime_denovo $SORTED_REPS --uchimeout $UCHIME

# swarm requires *.fas files to not have barcode labels on them, so we will need to remove them using something like:
if [! -d "$SWARM_EXTRACT_DIR" ]
then
    mkdir -p $SWARM_EXTRACT_DIR
    #sed 's/;barcodelabel=.*//' $EXTRACT_DIR/BCS*.fas >$SWARM_EXTRACT_DIR/ uncomment this line to use this script
fi

python OTU_contingency_table_simple.py $SORTED_REPS $SWARM_STATS $SWARMS $UCHIME for_swarm/DADA2_extracted_samples/BCS*.fas > ${OTU_TABLE}

rm $TMP_concat_derep $SWARM_REPS
