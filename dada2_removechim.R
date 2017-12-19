library(dada2); packageVersion("dada2")
# Merge multiple runs (if necessary)
st1 <- readRDS("BCS1_seqtab.rds")
st4 <- readRDS("BCS4_seqtab.rds")
st8 <- readRDS("BCS8_seqtab.rds")
st9 <- readRDS("BCS9_seqtab.rds")
st.all <- mergeSequenceTables(st1, st4, st8, st9)
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy
# Write to disk
saveRDS(seqtab, "seqtab_final.rds") # CHANGE ME to where you want sequence table saved

tax <- assignTaxonomy(seqtab, "/groups/cbi/bryan/BOLD_inverts/BOLD_COI_training.fasta", multithread=TRUE)
saveRDS(tax, "tax_final.rds") # CHANGE ME ...
