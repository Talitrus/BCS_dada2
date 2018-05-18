# load libraries --------------

library(openssl)


# Functions ------------------------------------
extrSamDADA2 <- function(my_table) {
  out_path <- file.path(getwd(), "DADA2_extracted_samples")
  if(!file_test("-d", out_path)) dir.create(out_path)
  for (sampleX in seq(1:dim(my_table)[1])){
    sinkname <- file.path(out_path, paste0(rownames(my_table)[sampleX],".fas"))
    {
      sink(sinkname)
      for (seqX in seq(1:dim(my_table)[2])) {
        if (my_table[sampleX,seqX] > 0) {
          header <- paste0(">",openssl::sha1(colnames(my_table)[seqX]),";size=",
                           my_table[sampleX,seqX],";barcodelabel=",rownames(my_table)[sampleX],";","\n")
          cat(header)
          seqq <- paste0(colnames(my_table)[seqX],"\n")
          cat(seqq)
        }
      }
      sink()
    }
  }
}

# generate tab-delimited OTU table with taxa as rows
#DADA2_otutab <- readRDS("seqtab_final.rds")
#LULU_DADA2_otutab <- t(DADA2_otutab)
#rownames should already be SHA1 sums
#print(rownames(LULU_DADA2_otutab)[1:5])
#rownames(LULU_DADA2_otutab) <- sha1(rownames(LULU_DADA2_otutab))
#write.table(LULU_DADA2_otutab, file = "DADA2_sha1_otutab.tsv", sep = "\t", quote = FALSE)

extrSamDADA2(DADA2_otutab)
