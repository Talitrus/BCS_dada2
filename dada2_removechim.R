library(dada2); packageVersion("dada2")
library(openssl)
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
st13 <- readRDS("BCS13_seqtab.rds")
st14 <- readRDS("BCS14_seqtab.rds")
st15 <- readRDS("BCS15_seqtab.rds")
st16 <- readRDS("BCS16_seqtab.rds")
st17 <- readRDS("BCS17_seqtab.rds")
st18 <- readRDS("BCS18_seqtab.rds")

st.all <- mergeSequenceTables(st1, st4, st8, st9, st10, st11, st6, st3, st12, st13, st14, st15, st16, st17, st18)
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy
# Write to disk
saveRDS(seqtab, "seqtab_final_original.rds") # CHANGE ME to where you want sequence table saved
#seqtab <- readRDS("seqtab_final.rds")
seqtab_sha1 <- seqtab
colnames(seqtab_sha1) <- sha1(colnames(seqtab_sha1))
saveRDS(seqtab_sha1, "seqtab_final.rds")
uniquesToFasta(seqtab, fout = "uniques.fasta")
#tax <- assignTaxonomy(seqtab, "../../Midori_COI/midori_COI_reformat.fasta", minBoot=70, multithread=TRUE, tryRC=TRUE)
#saveRDS(tax, "midori70RC_tax_final.rds") # CHANGE ME ...
