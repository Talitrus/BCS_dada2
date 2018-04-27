# BCS Pipeline Documentation
## Introduction
This README serves as a virtual lab notebook. We will document the scripts and the order that they're run in on here. This pipeline is built for use on a server that utilizes SLURM.

## Required software
- flexbar
- R 3.4.2 or later
  - DADA2
  - vegan
  - Plotly R package (properly configured for online use)
  - phyloseq
  - stringr
  - ggplot2
  - DESeq2
  - parallel
  - breakaway
- RDP Classifier
  - Java

## Setup
### Folder setup
Make a main project folder. Inside make subdirectories for each run (e.g. BCS1, BCS2, BCS3, BCS4, etc). The raw reads (FASTQ or FASTQ.GZ) can go into these subdirectories. Unzip any FASTQ.GZ files to maintain consistency across different folders. In the future, this step doesn't need to be done (.GZ files preferable), but minor changes need to be made to the flexbar scripts to accomodate this.
The main directory should also contain a tab-separated value file called `key.txt` that contains information corresponding library number, adapter index number (sample number from the Illumina run), and primer tag number to ML ID numbers and other sample metadata. Also in the main project folder, make another plaintext file call `libraries.txt` that contains the names of each of the run subdirectories (e.g. BCS1, BCS2, BCS3, etc), which one name on each line.
### Primer barcodes
Make a FASTA file (or multiple FASTA files if multiple barcoding schemes are used) containing the primer index barcodes in the main project folder. This is called `ML_barcodes.fasta` in these scripts.
### Adapter FASTA file
You will need to make a similar FASTA file containing all of the adapters that you wish flexbar to trim for. Here, it is called `truseq_adapters.fasta`.

### Script setup
This repository can be cloned into its own folder inside of the main project folder (name it "scripts" or something). Inside of your scripts folder, make a subdirectory called `out_err_files` to contain log files and error logs.

## Run Order
1. make\_file\_root\_list.sh
2. flexbar.sh
3. flexbar2.sh

## Cleanup
Move all unassigned reads into a new subdirectory called `unassigned`, if any are remaining for whatever reason after running `flexbar2.sh`. This should be performed automatically by the script now.

## R Analysis with DADA2 and RDP
### Description

`rplce\_header\_w\_seq.sh` replaces the sequence names in the unique sequence FASTA file (generated in `removechim.sh`).
`vsearch\_cluster.sh` then generates the clusters (UCLUST-formated output) and processes them in preparation for merging taxa in R with phyloseq.
`gen\_matchlist.sh` will search the centroids from clustering against themselves to generate the `matchlist.txt` file for LULU.
The `phyloseq.sh` and `BCS_phyloseq.R` scripts run community analyses and visualization (based on the Plotly tool). You may need to make major changes to this section for your own analysis or configure Plotly if you want to use the existing visualizations. For some reason, the RDP implementation in DADA2 doesn't perform well for us, possibly due to memory or scaling issues. We run RDP outside of R/DADA2 to get around this.

### Setup
R scripts can be placed in the same folder as the other scripts (in a folder inside of the folders with all the directories containing the sequenced libraries). Code for renaming files generated by flexbar2.sh will be contained here as well as R code for running DADA2 analysis.

### Run Order
0. rename.sh
1. filtertrim.sh
2. removechim.sh
3. RDP\_classify.sh (or BLCAb70.sh)
4. rplce\_header\_w\_seq.sh
5. vsearch\_cluster.sh
6. gen\_matchlist.sh
7. phyloseq.sh
