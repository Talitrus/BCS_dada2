library(dada2)
library(stringr)
setwd("/groups/cbi/bryan/BCS_all/dada2_R")
BCS_folder <- file.path("../consolidated")
seq_files <- list.files(path = "../consolidated", full.names = TRUE)
filt_paths <- file.path(path = BCS_folder, 'filtered' )
BCS_fnFs <- sort(list.files(BCS_folder, pattern="_1.fastq", full.names = TRUE))
BCS_fnRs <- sort(list.files(BCS_folder, pattern="_2.fastq", full.names = TRUE))
if(length(BCS_fnFs) != length(BCS_fnRs)) stop("Forward and reverse files do not match.")

BCS_fnFs_root <- substr(sort(list.files(BCS_folder, pattern="_1.fastq", full.names = FALSE)),1,nchar(sort(list.files(BCS_folder, pattern="_1.fastq", full.names = FALSE)))-8)

BCS_filtFs <- file.path(filt_paths, paste0(BCS_fnFs_root, '_F_filt.fastq.gz'))
BCS_filtRs <- file.path(filt_paths, paste0(BCS_fnFs_root, '_R_filt.fastq.gz'))
names(BCS_filtFs) <- BCS_fnFs_root
names(BCS_filtRs) <- BCS_fnFs_root

BCS_out <- filterAndTrim(BCS_fnFs, BCS_filtFs, BCS_fnRs, BCS_filtRs, rm.phix = TRUE, maxN = 0, maxEE = c(2,2), truncQ = 10, trimLeft = 26, multithread = TRUE) #set multithread to TRUE when running in normal R and not RStudio
BCS_out
lib_numbers <- as.factor(str_match(BCS_fnFs_root,"(?<=BCS)[0-9]+"))
lib_names <- unique(paste0("BCS",lib_numbers))
set.seed(100)



generate_RDS <- function(curr_lib_num, path = '') {
  errF <- learnErrors(BCS_filtFs[which(lib_numbers==curr_lib_num)], nreads = 1e6, multithread = TRUE) # change all multithread to TRUE on cluster
  errR <- learnErrors(BCS_filtRs[which(lib_numbers==curr_lib_num)], nreads = 1e6, multithread = TRUE)
  
  
  mergers <- vector("list", length(which(lib_numbers==curr_lib_num)))
  names(mergers) <- BCS_fnFs_root[which(lib_numbers==curr_lib_num)]
  
  
  #plotErrors(errF, nominalQ = TRUE)
  #plotErrors(errR, nominalQ = TRUE)
  
  for(sam in BCS_fnFs_root[which(lib_numbers==curr_lib_num)]) {
    cat("Processing:", sam, "\n")
    derepF <- derepFastq(BCS_filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(BCS_filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
  }
  
  rm(derepF); rm(derepR); rm(ddF); rm(ddR)
  
  # Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, paste0(path, "BCS",curr_lib_num,"_seqtab.rds"))
  
  rm(seqtab)
}

#BCS 1
generate_RDS(1)
generate_RDS(4)
generate_RDS(8)
generate_RDS(9)

