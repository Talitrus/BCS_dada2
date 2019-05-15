#!/bin/bash

while read f; do
	cd ../$f
	mkdir -p unassigned
	mv Undeter* unassigned/
	ls *.fastq* | grep -v barcode | sed -E 's/[12]_001\..*//' | uniq > file_roots.txt
done <../libraries.txt
