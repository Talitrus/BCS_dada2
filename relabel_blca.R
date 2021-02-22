library(tidyverse)
make_blca_tax_table <- function(in_file, min_confidence = 0) {
  class.table <- read.table(textConnection(gsub("[\t:]", ";", readLines(in_file))), sep = ";", fill = TRUE, na.strings = c("", "NA"), col.names = c("SeqID", "level1", "Superkingdom", "l1conf", "level2", "Kingdom", "l2conf", "level3", "Phylum", "l3conf", "level4", "Class", "l4conf", "level5", "Order", "l5conf", "level6", "Family", "l6conf", "level7", "Genus", "l7conf", "level8", "species", "l8conf", "blank"))
  class.tib <- as.tibble(class.table[,1:25])
  rank.key <- c(3,6,9,12,15,18,21,24)
  names(rank.key) <- c("Superkingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  # Use confidence threshold to set low levels to NA
  if ((min_confidence > 0) & any(class.tib[,rank.key+1] < (min_confidence*100), na.rm = TRUE)) {
    ind_below_conf <- which(class.tib[,rank.key+1] < (min_confidence*100), arr.ind = TRUE)
    for (i in 1:nrow(ind_below_conf)) {
      ind_row <- ind_below_conf[i,1]
      ind_col <- ind_below_conf[i,2]
      class.tib[ind_row, ind_col*3] <- NA
    }
  }
  class.tib.mod <- mutate(class.tib, SeqInt = as.integer(substr(as.character(SeqID), 3, nchar(as.character(SeqID))))) #Assumes all sequence IDs start with two letters before the informative number.
  return(select(arrange(class.tib.mod, SeqInt), Superkingdom, Kingdom, Phylum, Class, Order, Family, Genus, species))
}

sha1_labels <- read_lines("unique_sha1_labels.txt")

blca_table <- make_blca_tax_table("uniques_trim.fasta.blca.out", min_confidence=0.5)
if(nrow(blca_table) == length(sha1_labels)) {
  blca_table$sha1 <- sha1_labels
  write_tsv(blca_table[,c(8,1:7)], "sha1_blca.tsv")
} else {
  stop("BLCA output and sha1 labels are of different lengths.")
}
