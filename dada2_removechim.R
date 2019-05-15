library(dada2); packageVersion("dada2")
library(openssl)
# Merge multiple runs (if necessary)
st2 <- readRDS("BCS2_seqtab.rds")
#st.all <- mergeSequenceTables(st2)

# Remove chimeras
seqtab <- removeBimeraDenovo(st2, method="consensus", multithread=TRUE)
# Assign taxonomy
# Write to disk
saveRDS(seqtab, "seqtab_final_original.rds") # CHANGE ME to where you want sequence table saved
#seqtab <- readRDS("seqtab_final.rds")
seqtab_sha1 <- seqtab
colnames(seqtab_sha1) <- sha1(colnames(seqtab_sha1))
saveRDS(seqtab_sha1, "seqtab_final.rds")
uniquesToFasta(seqtab, fout = "uniques.fasta")
tax <- assignTaxonomy(seqtab, "../../silva/dada2/silva_nr_v132_train_set.fa.gz", minBoot=70, multithread=TRUE, tryRC=TRUE)
saveRDS(tax, "silva_tax.rds") # CHANGE ME ...
