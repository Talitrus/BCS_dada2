#!/bin/bash
while read f; do
	cd ../$f
	ls *.fastq* | grep -v barcode | sed -E 's/[12]_001\..*//' | uniq > file_roots.txt
done <../libraries.txt
