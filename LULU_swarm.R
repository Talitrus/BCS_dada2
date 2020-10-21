#!/usr/bin/env Rscript

library(lulu)

# Set working directory
setwd("/lustre/groups/cbi/Users/bnguyen/bocas/BCS_all/dada2_R/")

# Read in DADA2 (+VSEARCH) OTU table
# To use just raw DADA2 load in seqtab_final.rds and transpose it.
# Otherwise use "FINAL_FOR_LULU_otutab.txt"

dada2_vsearch_otutab <- read.delim(file = "swarm_13_trim.otutable", quote = "", row.names = 1) # Load OTU table from DADA2 + VSEARCH
names(dada2_vsearch_otutab) <- gsub("\\.(?=.*_)","-",names(dada2_vsearch_otutab), perl = TRUE) # At some point in the LULU pipeline, all of the hyphens got replaced with periods, so we're just going to go back in and replace the periods with hyphens again to maintain compatibility.

matchlist_LULU <- read.delim("swarm_matchlist.txt", stringsAsFactors = FALSE) # generate matchlist.txt with BLAST (preferably?) or VSEARCH
curated_LULU <- lulu(as.data.frame(dada2_vsearch_otutab), matchlist_LULU, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95)
saveRDS(curated_LULU, file = "LULU_curation_swarm.rds")
