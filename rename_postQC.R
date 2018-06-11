
setwd("/groups/cbi/bryan/BCS_all/dada2_R")

library(stringr)
BCS_list <- list.dirs('..', recursive = FALSE)[grep("BCS", list.dirs('..', recursive = FALSE))]
BCS_fnFs <- lapply(BCS_list, function(x) sort(list.files(x, pattern="final_barcode.+_2.fastq", full.names = TRUE, recursive = FALSE))) # The forward read is actually labeled as read 2 as a result of running Flexbar by feeding the reverse read in first.
BCS_fnRs <- lapply(BCS_list, function(x) sort(list.files(x, pattern="final_barcode.+_1.fastq", full.names = TRUE, recursive = FALSE)))
libkey <- read.delim("../key.txt", header = TRUE)
fnFs <- as.character(unlist(BCS_fnFs))
fnRs <- as.character(unlist(BCS_fnRs))
file_roots <- substr(fnFs,1,nchar(fnFs)-8)
#make sure all unassigned barcode files moved out into a different folder.
adapter <- as.integer(str_match(file_roots,"(?<=_S)[0-9]+(?=_)"))
primer_tag <- as.integer(str_match(file_roots,"(?<=barcode_)[0-9]+"))
lib_number <- as.integer(str_match(file_roots,"(?<=BCS)[0-9]+"))

get_MLID <- function (lib_num, ad_num, prim_num) {
  return(as.character(libkey[libkey[,"Library"] == lib_num & libkey[,"Adapter.Order"] == ad_num & libkey[,"Primer.Tag.Number"] == prim_num,][,'MLID']))
}

#Check for any mismatched forward/reverse files
file_roots_R <- substr(fnRs,1,nchar(fnRs)-8)
if (any(!(file_roots_R %in% file_roots)) | any(!(file_roots %in% file_roots_R))) {
  print("Mismatched post-Flexbar forward and reverse reads. No corresponding pair found for files:")
  print(which(!(file_roots_R %in% file_roots)))
  print(which(!(file_roots %in% file_roots_R)))
  print("Continuing anyway...")
  }

fnFs <- fnFs[file_roots %in% file_roots_R]
fnRs <- fnRs[file_roots_R %in% file_roots]
read.df <- data.frame(fnFs, fnRs, file_roots, lib_number, adapter, primer_tag)
read.df$MLID <- apply(read.df[,c('lib_number', 'adapter','primer_tag')], 1, function (x) get_MLID(x[1],x[2],x[3]))
read.df$uniqID <- paste0('BCS',read.df$lib_number,'-',read.df$adapter,'-',read.df$primer_tag,'_',read.df$MLID)
  
file.copy( from = as.character(read.df$fnFs), to = paste0('../consolidated/',read.df$uniqID,'_1.fastq'))
file.copy( from = as.character(read.df$fnRs), to = paste0('../consolidated/',read.df$uniqID,'_2.fastq'))
