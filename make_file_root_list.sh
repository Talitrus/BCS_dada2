#!/bin/bash
while read f; do
	cd ../$f
	ls *.fastq | sed -E 's/[12]_001\.fastq//' | uniq > file_roots.txt
done <../libraries.txt
