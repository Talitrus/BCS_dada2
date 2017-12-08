library(dada2)

BCS_list <- list.dirs('..')[grep("BCS", list.dirs('..'))]
#We will create or use an existing directory called 'filtered' in each BCS folder
#Move Undetermined or unassigned sequences to a subdirectory.
filt_paths <- file.path(BCS_list, 'filtered')
#Maybe loop starting here
#BCS1_files = list.files(path = BCS_list[1], full.names = TRUE)
#BCS1_fnFs <- sort(list.files(BCS_list[1], pattern="_1.fastq", full.names = TRUE))
#BCS1_fnRs <- sort(list.files(BCS_list[1], pattern="_2.fastq", full.names = TRUE))
BCS_fnFs <- lapply(BCS_list, function(x) sort(list.files(x, pattern="_1.fastq", full.names = TRUE)))
BCS_fnRs <- lapply(BCS_list, function(x) sort(list.files(x, pattern="_2.fastq", full.names = TRUE)))

sample.names.1 <- c('a','b') #fix this with real sample names later
#plotQualityProfile(BCS1_fnFs[1:2]) #trim last ~2-5  nucleotides?
#plotQualityProfile(BCS1_fnRs[1:2]) #Could trim about ~5-10 also
BCS1_filtFs <- file.path(filt_paths[1], paste0(sample.names.1, '_F_filt.fastq.gz'))
BCS1_filtRs <- file.path(filt_paths[1], paste0(sample.names.1, '_R_filt.fastq.gz'))

BCS1_out <- filterAndTrim(unlist(BCS_fnFs[1])[c(34,19)], BCS1_filtFs, unlist(BCS_fnRs[1])[c(34,19)], BCS1_filtRs, rm.phix = TRUE, maxN = 0, maxEE = c(2,2), truncQ = 10, trimLeft = 32, multithread = TRUE) #set multithread to TRUE when running in normal R and not RStudio
