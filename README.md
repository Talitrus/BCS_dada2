# BCS Pipeline Documentation

This README serves as a virtual lab notebook. We will document the scripts and the order that they're run in on here.

## Setup
### Folder setup
Make a main project folder. Inside make subdirectories for each run (e.g. BCS1, BCS2, BCS3, BCS4, etc). The raw reads (FASTQ or FASTQ.GZ) can go into these subdirectories. Unzip any FASTQ.GZ files to maintain consistency across different folders. In the future, this step doesn't need to be done (.GZ files preferable), but minor changes need to be made to the flexbar scripts to accomodate this.
The main directory should also contain a tab-separated value file called `key.txt` that contains information corresponding library number, adapter index number (sample number from the Illumina run), and primer tag number to ML ID numbers and other sample metadata. Also in the main project folder, make another plaintext file call `libraries.txt` that contains the names of each of the run subdirectories (e.g. BCS1, BCS2, BCS3, etc), which one name on each line.

### Script setup
This repository can be cloned into its own folder inside of the main project folder (name it "scripts" or something). Inside of your scripts folder, make a subdirectory called `out_err_files` to contain log files and error logs.

## Run Order
1. make\_file\_root\_list.sh
2. flexbar.sh
3. flexbar2.sh

## Cleanup
Move all unassigned reads into a new subdirectory called `unassigned`, if any are remaining for whatever reason after running `flexbar2.sh`. This should be performed automatically by the script now.

## R Analysis with DADA2
### Description
The `phyloseq.sh` and `BCS_phyloseq.R` scripts run community analyses and visualization (based on the Plotly tool). You may need to make major changes to this section for your own analysis or configure Plotly if you want to use the existing visualizations.

### Setup
R scripts can be placed in the same folder as the other scripts (in a folder inside of the folders with all the directories containing the sequenced libraries). Code for renaming files generated by flexbar2.sh will be contained here as well as R code for running DADA2 analysis.

### Run Order
0. rename.sh
1. filtertrim.sh
2. removechim.sh
3. phyloseq.sh
