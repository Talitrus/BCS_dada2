library(dada2); packageVersion("dada2")
library(openssl)
# Merge multiple runs (if necessary)
st2 <- readRDS("BCS2_seqtab.rds")
st5 <- readRDS("BCS5_seqtab.rds")
st7 <- readRDS("BCS7_seqtab.rds")
st21 <- readRDS("BCS21_seqtab.rds")

st.all <- mergeSequenceTables(st2,st5,st7,st21)

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
tax <- assignTaxonomy(seqtab_sha1, "/groups/cbi/Users/bnguyen/db/dada2/silva_18s/silva_132.18s.99_rep_set.dada2.fa.gz", minBoot=70, multithread=TRUE, tryRC=TRUE) #test out assigning with SHA1 names? Will it work?
saveRDS(tax, "silva_tax.rds") # CHANGE ME ...
