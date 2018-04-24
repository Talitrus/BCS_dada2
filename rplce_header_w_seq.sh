#!/bin/bash
perl -0pe 's/>.+(;size=[0-9]+;)\n(.+)/>$2$1\n$2/g' uniques.fasta > uniques_seqheaders.fasta
