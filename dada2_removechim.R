library(dada2); packageVersion("dada2")
# Merge multiple runs (if necessary)
st1 <- readRDS("/path/to/run1/output/seqtab.rds")
st2 <- readRDS("/path/to/run2/output/seqtab.rds")
st3 <- readRDS("/path/to/run3/output/seqtab.rds")
st.all <- mergeSequenceTables(st1, st2, st3)
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy
tax <- assignTaxonomy(seqtab, "/path/to/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
# Write to disk
saveRDS(seqtab, "/path/to/study/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
saveRDS(tax, "/pathto/study/tax_final.rds") # CHANGE ME ...