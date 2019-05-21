library(plotly)
library(dada2)
library(stringr)
setwd("/groups/cbi/Users/bnguyen/bocas/BCS_18S/scripts")
BCS_folder <- file.path("../consolidated")
seq_files <- list.files(path = "../consolidated", full.names = TRUE)
filt_paths <- file.path(path = BCS_folder, 'filtered' )
BCS_fnFs <- sort(list.files(BCS_folder, pattern="_1.fastq", full.names = TRUE))
BCS_fnRs <- sort(list.files(BCS_folder, pattern="_2.fastq", full.names = TRUE))
if(length(BCS_fnFs) != length(BCS_fnRs)) stop("Forward and reverse files do not match.")

BCS_fnFs_root <- substr(sort(list.files(BCS_folder, pattern="_1.fastq", full.names = FALSE)),1,nchar(sort(list.files(BCS_folder, pattern="_1.fastq", full.names = FALSE)))-11)

BCS_filtFs <- file.path(filt_paths, paste0(BCS_fnFs_root, '_F_filt.fastq.gz'))
BCS_filtRs <- file.path(filt_paths, paste0(BCS_fnFs_root, '_R_filt.fastq.gz'))
names(BCS_filtFs) <- BCS_fnFs_root
names(BCS_filtRs) <- BCS_fnFs_root

BCS_out <- filterAndTrim(BCS_fnFs, BCS_filtFs, BCS_fnRs, BCS_filtRs, rm.phix = TRUE, maxN = 0, truncQ = 2, maxEE = c(5, 10), truncLen = c(270,220), multithread = TRUE, verbose = TRUE) #set multithread to TRUE when running in normal R and not RStudio
BCS_out
lib_numbers <- as.factor(str_match(BCS_fnFs_root,"(?<=BCS)[0-9]+"))
lib_names <- unique(paste0("BCS",lib_numbers))


set.seed(100)



generate_RDS <- function(curr_lib_num, path = '', checkforRC = TRUE, refseq = "ATGCATGTCTAAGTTCACACTGTTTCACGGTGAAACCGCGAATGGCTCATTAAATCAGTCGAGGTTCCTTAGATGACACGATCCTACTTGGATAACTGTGGCAATTCTAGAGCTAATACATGCTTCGCAAGCTCCGACCCTCGTGGAAAGAGCGCTTTTATTAGTTCAAAACCAATCGTCGTTTCGCCGTTTCGGCGGCGGGCGGCGTCCCAATTGGTGACTCTGGATAACTTTGTGCTGATCGCATGGCCTTTTGTGCCGGCGACGCATCTTTCAAATGTCTGCCCTATCAAATGTCGATGGTACGTGATATGCCTACCATGTTTGTAACGGGTAACGGGGAATCAGGGTTCGATTCCGGAGAGGGAGCATGAGAAACGGCTACCACA") {
  errF <- learnErrors(BCS_filtFs[which(lib_numbers==curr_lib_num)], nbases = 1e8, multithread = TRUE) # change all multithread to TRUE on cluster
  errR <- learnErrors(BCS_filtRs[which(lib_numbers==curr_lib_num)], nbases = 1e8, multithread = TRUE)
  
  
  mergers <- vector("list", length(which(lib_numbers==curr_lib_num)))
  names(mergers) <- BCS_fnFs_root[which(lib_numbers==curr_lib_num)]
  
  
  f_err_plot <- plotErrors(errF, nominalQ = TRUE) + ggtitle(paste0("BCS", curr_lib_num, " forward"))
  r_err_plot <- plotErrors(errR, nominalQ = TRUE) + ggtitle(paste0("BCS", curr_lib_num, " reverse"))
  api_create(ggplotly(f_err_plot), filename = paste0("bocas/18S/dada2/BCS", curr_lib_num, "_F_errs"), sharing = "secret")
  api_create(ggplotly(r_err_plot), filename = paste0("bocas/18S/dada2/BCS", curr_lib_num, "_R_errs"), sharing = "secret")
  ggsave(paste0("../analysis/dada2/BCS", curr_lib_num, "_F_errs.png"), device = "png", plot = f_err_plot)
  ggsave(paste0("../analysis/dada2/BCS", curr_lib_num, "_R_errs.png"), device = "png", plot = r_err_plot)

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
  
  if(checkforRC) {
    score <- c()
    for(sam in names(mergers)) {
      seq.hd <- nwhamming(substr(refseq,1,70), substr(mergers[[sam]][,1],1,70), vec=TRUE, endsfree=F)
      seq.rc <- nwhamming(substr(dada2:::rc(refseq),1,70), substr(mergers[[sam]][,1],1,70), vec=TRUE, endsfree=F)
      seq.df <- seq.rc - seq.hd #positive = original orientation, negative = reverse orientation
      score <- c(score, seq.df)
      #api_create(ggplotly(orientation_hist), filename = paste0("bocas/18S/dada2/BCS", curr_lib_num, "_F_orientation_hist"), sharing = "secret")
      #ggsave(paste0("../analysis/dada2/BCS", curr_lib_num, "_orientation_hist.png"), device = "png", plot = orientation_hist)
      seq.r <- which(seq.df < 0)
      mergers[[sam]][seq.r,1] <- dada2:::rc(mergers[[sam]][seq.r,1])   
    }
    orientation_hist <- ggplot(data.frame(score = score), aes(x = score)) + geom_histogram(aes(y = ..density..)) + geom_density(alpha = 0.2, fill="#FF6666") + ggtitle(paste0("BCS", curr_lib_num, "sequence scores"))
    api_create(ggplotly(orientation_hist), filename = paste0("bocas/18S/dada2/BCS", curr_lib_num, "_F_orientation_hist"), sharing = "secret")
    #ggsave(paste0("../analysis/dada2/BCS", curr_lib_num, "_orientation_hist.png"), device = "png", plot = orientation_hist)
    } 
  # Construct sequence table
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, paste0(path, "BCS",curr_lib_num,"_seqtab.rds"))
  
  rm(seqtab)
}

#generate whatever
generate_RDS(2)
generate_RDS(5)
generate_RDS(7)
generate_RDS(21)
