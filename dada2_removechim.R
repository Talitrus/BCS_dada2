library(ggplot2)
library(magrittr)
library(plotly)
library(dada2); packageVersion("dada2")
library(openssl)

# Functions --------------------------------------

getAmpliconSizes <- function(st) {
    return(nchar(getSequences(st)))
}

plotAmpSizes <- function(st) {
    size_vec <- getAmpliconSizes(st)
    size_df <- data.frame(bp = size_vec)
    fig <- ggplot(size_df, aes(x = bp)) +
        geom_histogram(aes( y = ..density..)) +
        geom_density(alpha = 0.35, fill = "#6666FF")
    return(fig)
}

stSizeSelect <- function(st, l_cutoff, h_cutoff) {
    if (l_cutoff < h_cutoff) {
        st2 <- st[,nchar(colnames(seqtab)) %in% l_cutoff:h_cutoff]
        return(st2)
    } else {
        stop("Error: Lower cutoff must be an integer that is lower than the higher cutoff.")
    }
}



# Merge multiple runs (if necessary)
st2 <- readRDS("BCS2_seqtab.rds")
st5 <- readRDS("BCS5_seqtab.rds")
st7 <- readRDS("BCS7_seqtab.rds")
st21 <- readRDS("BCS21_seqtab.rds")

# Plot size distributions
st2.p <- plotAmpSizes(st2) %>%
    ggplotly()

st5.p <- plotAmpSizes(st5) %>%
    ggplotly()

st7.p <- plotAmpSizes(st7) %>%
    ggplotly()

st21.p <- plotAmpSizes(st21) %>%
    ggplotly()


api_create(st2.p, filename = "bocas/18S/dada2/BCS2_frag_sizes", sharing = "secret")
api_create(st5.p, filename = "bocas/18S/dada2/BCS5_frag_sizes", sharing = "secret")
api_create(st7.p, filename = "bocas/18S/dada2/BCS7_frag_sizes", sharing = "secret")
api_create(st21.p, filename = "bocas/18S/dada2/BCS21_frag_sizes", sharing = "secret")

st.all <- mergeSequenceTables(st2,st5,st7,st21)

# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=parallel::detectCores())
# Assign taxonomy
# Write to disk
saveRDS(seqtab, "seqtab_final_original.rds") # CHANGE ME to where you want sequence table saved
#seqtab <- readRDS("seqtab_final.rds")
seqtab_sha1 <- seqtab
colnames(seqtab_sha1) <- sha1(colnames(seqtab_sha1))
saveRDS(seqtab_sha1, "seqtab_final.rds")
uniquesToFasta(seqtab, fout = "uniques.fasta")
tax <- assignTaxonomy(seqtab, "/groups/cbi/Users/bnguyen/db/dada2/silva_18s/silva_132.18s.99_rep_set.dada2.fa.gz", minBoot=70, multithread=parallel::detectCores(), tryRC=TRUE)
rownames(tax) <- sha1(rownames(tax))
saveRDS(tax, "silva_tax.rds") # CHANGE ME ...
