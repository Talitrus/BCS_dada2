library(dada2); packageVersion("dada2")
# Merge multiple runs (if necessary)
st1 <- readRDS("BCS1_seqtab.rds")
st4 <- readRDS("BCS4_seqtab.rds")
st8 <- readRDS("BCS8_seqtab.rds")
st9 <- readRDS("BCS9_seqtab.rds")
st10 <- readRDS("BCS10_seqtab.rds")
st6 <- readRDS("BCS6_seqtab.rds")
st11 <- readRDS("BCS11_seqtab.rds")
st3 <- readRDS("BCS3_seqtab.rds")
st12 <- readRDS("BCS12_seqtab.rds")
st.all <- mergeSequenceTables(st1, st4, st8, st9, st10, st11, st6, st3, st12)
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy
# Write to disk
saveRDS(seqtab, "seqtab_final.rds") # CHANGE ME to where you want sequence table saved
seqtab <- readRDS("seqtab_final.rds")
uniquesToFasta(seqtab, fout = "uniques.fasta")
#tax <- assignTaxonomy(seqtab, "../../Midori_COI/midori_COI_reformat.fasta", minBoot=70, multithread=TRUE, tryRC=TRUE)
#saveRDS(tax, "midori70RC_tax_final.rds") # CHANGE ME ...
