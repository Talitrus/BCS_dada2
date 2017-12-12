library(dada2)
setwd("/Users/bryan/Documents/Grad/STRI/prototyping/R")
BCS_folder <- file.path("../consolidated")
seq_files <- list.files(path = "../consolidated", full.names = TRUE)
filt_paths <- file.path(path = BCS_folder, 'filtered' )
BCS_fnFs <- sort(list.files(BCS_folder, pattern="_1.fastq", full.names = TRUE))
BCS_fnRs <- sort(list.files(BCS_folder, pattern="_2.fastq", full.names = TRUE))
if(length(BCS_fnFs) != length(BCS_fnRs)) stop("Forward and reverse files do not match.")

BCS_fnFs_root <- substr(sort(list.files(BCS_folder, pattern="_1.fastq", full.names = FALSE)),1,nchar(sort(list.files(BCS_folder, pattern="_1.fastq", full.names = FALSE)))-8)

BCS_filtFs <- file.path(filt_paths, paste0(BCS_fnFs_root, '_F_filt.fastq.gz'))
BCS_filtRs <- file.path(filt_paths, paste0(BCS_fnFs_root, '_R_filt.fastq.gz'))

BCS_out <- filterAndTrim(BCS_fnFs, BCS_filtFs, BCS_fnRs, BCS_filtRs, rm.phix = TRUE, maxN = 0, maxEE = c(2,2), truncQ = 10, trimLeft = 26, multithread = FALSE) #set multithread to TRUE when running in normal R and not RStudio
