#!/usr/bin/env Rscript

# pass an integer (number of cores) to the Rscript to use multiple cores. Otherwise it defaults to 1 core.
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	nslots <- 1
	} else if (length(args)==1 & is.integer(as.integer(args[1]))) {
	# default output file
	nslots <- as.integer(args[1])
	} else {
	stop("More than one argument supplied")
	}

# Load libraries ------------------------

library(phyloseq)
library(vegan)
library(ggplot2)
#library(DESeq2) #DESeq2 library not working on Hydra?
library(parallel)
library(plotly)
library(breakaway)
library(RDPutils)
library(igraph)
library(dplyr)
library(tibble)
library(readr)
library(stringr)


# Functions ----------

make_blca_tax_table <- function(in_file, min_confidence = 0) {
  class.table <- read.table(textConnection(gsub("[\t:]", ";", readLines(in_file))), sep = ";", fill = TRUE, na.strings = c("", "NA"), col.names = c("SeqID", "level1", "Superkingdom", "l1conf", "level2", "Phylum", "l2conf", "level3", "Class", "l3conf", "level4", "Order", "l4conf", "level5", "Family", "l5conf", "level6", "Genus", "l6conf", "level7", "species", "l7conf", "blank"))
  class.tib <- as.tibble(class.table[,1:22])
  rank.key <- c(3,6,9,12,15,18,21)
  names(rank.key) <- c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
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
  return(select(arrange(class.tib.mod, SeqInt), Superkingdom, Phylum, Class, Order, Family, Genus, species))
}

vegan_otu <- function(physeq) { #convert phyloseq OTU table into vegan OTU matrix
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


plotly_ps_network <- function(ps.object, distance = "bray", coords, def_threshold = 0.6625, segment_col = I("blue"), marker_col = "orange") {
  aval <- list()
  start_step <- round((def_threshold*100-50)*2/2.5)
  sam.df <- sample_data(ps.object)
  lowest_step = 1
  sam.df$longitude <- coords[sam.df$Site.Code,"longitude"]
  sam.df$latitude <- coords[sam.df$Site.Code,"latitude"]
  for(step in 1:40){
    dist.thres <- (step*2.5/2+50)/100
    igr <- make_network(ps.object, distance = "bray", max.dist = dist.thres)
    igr.edges <- as.tibble(as_edgelist(igr))
    #coral500.sam.df$LongID <- rownames(coral500.sam.df)
    
    #add coordinates to sample dataframe
    if(dim(igr.edges)[1] > 0) {
      igr.edges$slat <- unlist(sam.df[igr.edges$V1,'latitude'])
      igr.edges$slon <- unlist(sam.df[igr.edges$V1,'longitude'])
      igr.edges$elat <- unlist(sam.df[igr.edges$V2,'latitude'])
      igr.edges$elon <- unlist(sam.df[igr.edges$V2,'longitude'])
      aval[[step]] <-list(visible = FALSE,
                          name = step*2.5,
                          slat=igr.edges$slat,
                          slon=igr.edges$slon,
                          elat=igr.edges$elat,
                          elon=igr.edges$elon
      )
    }
    else {
      aval[[step]] <-list(visible = FALSE,
                          name = step*2.5,
                          slat=NULL,
                          slon=NULL,
                          elat=NULL,
                          elon=NULL
                          
      )
      lowest_step = lowest_step + 1
    }
  }
  aval[start_step][[1]]$visible = TRUE
  steps <- list()
  p <- plot_mapbox()
  for (i in lowest_step:40) {
    p <- add_segments(p,
                      x = aval[i][[1]]$slon, xend = aval[i][[1]]$elon,
                      y = aval[i][[1]]$slat, yend = aval[i][[1]]$elat,
                      visible = aval[i][[1]]$visible,
                      opacity = 0.3, name = "Connections",
                      hoverinfo = "none", showlegend = FALSE,
                      color = segment_col)
    step <- list(args = list('visible', c(rep(FALSE, length(lowest_step:40)), TRUE)),
                 method = 'restyle', label = (i*2.5/2+50)/100)
    step$args[[2]][i] = TRUE  
    steps[[i]] = step 
  }
  p <- p %>%
    layout(sliders = list(list(active = start_step,
                               currentvalue = list(prefix = "Threshold: "),
                               steps = steps))) %>%
    add_markers(data = sam.df, x = ~longitude, y = ~latitude, size = I(10), opacity = 0.5, text = ~Site.Code, marker = list(color = marker_col), name = "Sites", visible = TRUE) %>%
    layout(mapbox = list(
      zoom = 10,
      center = list(lat = 9.302979, lon = -82.231656),
      style = "mapbox://styles/nguyenbn/cjdnkuoo306vh2rnycjavmito"
    )
    )
  
  return(p)
}

plotly_ps_network_by_samplenames <- function(ps.object, distance = "bray", coords, def_threshold = 0.6625, segment_col = I("blue"), marker_col = "orange") {
  aval <- list()
  start_step <- round((def_threshold*100-50)*2/2.5)
  lowest_step = 1
  sam.df <- sample_data(ps.object)
  sam.df$longitude <- coords[rownames(sam.df),"longitude"]
  sam.df$latitude <- coords[rownames(sam.df),"latitude"]
  for(step in 1:40){
    dist.thres <- (step*2.5/2+50)/100
    igr <- make_network(ps.object, distance = "bray", max.dist = dist.thres)
    igr.edges <- as.tibble(as_edgelist(igr))
    #coral500.sam.df$LongID <- rownames(coral500.sam.df)
    
    #add coordinates to sample dataframe
    if(dim(igr.edges)[1] > 0) {
      igr.edges$slat <- unlist(sam.df[igr.edges$V1,'latitude'])
      igr.edges$slon <- unlist(sam.df[igr.edges$V1,'longitude'])
      igr.edges$elat <- unlist(sam.df[igr.edges$V2,'latitude'])
      igr.edges$elon <- unlist(sam.df[igr.edges$V2,'longitude'])
      aval[[step]] <-list(visible = FALSE,
                          name = step*2.5,
                          slat=igr.edges$slat,
                          slon=igr.edges$slon,
                          elat=igr.edges$elat,
                          elon=igr.edges$elon
      )
    }
    else {
      aval[[step]] <-list(visible = FALSE,
                          name = step*2.5,
                          slat=NULL,
                          slon=NULL,
                          elat=NULL,
                          elon=NULL
                          
      )
      lowest_step = lowest_step + 1
    }
  }
  aval[start_step][[1]]$visible = TRUE
  steps <- list()
  p <- plot_mapbox()
  for (i in lowest_step:40) {
    p <- add_segments(p,
                      x = aval[i][[1]]$slon, xend = aval[i][[1]]$elon,
                      y = aval[i][[1]]$slat, yend = aval[i][[1]]$elat,
                      visible = aval[i][[1]]$visible,
                      opacity = 0.3, name = "Connections",
                      hoverinfo = "none", showlegend = FALSE,
                      color = segment_col)
    step <- list(args = list('visible', c(rep(FALSE, length(lowest_step:40)), TRUE)),
                 method = 'restyle', label = (i*2.5/2+50)/100)
    step$args[[2]][i] = TRUE  
    steps[[i]] = step 
  }
  p <- p %>%
    layout(sliders = list(list(active = start_step,
                               currentvalue = list(prefix = "Threshold: "),
                               steps = steps))) %>%
    add_markers(data = sam.df, x = ~longitude, y = ~latitude, size = I(10), opacity = 0.5, text = rownames(sam.df), marker = list(color = marker_col), name = "Sites", visible = TRUE) %>%
    layout(mapbox = list(
      zoom = 10,
      center = list(lat = 9.302979, lon = -82.231656),
      style = "mapbox://styles/nguyenbn/cjdnkuoo306vh2rnycjavmito"
    )
    )
  
  return(p)
}

# Working code --------------------------

setwd("/data/genomics/leraynguyen/R/code/")
#tax <- make_tax_table(in_file="BCS_RDP_output.txt", min_confidence = 0.5)
tax.tib <- make_blca_tax_table(in_file = "uniques_b70.blca.out", min_confidence = 0.5)
tax <- as.matrix(tax.tib)
seqtab <- readRDS("seqtab_final.rds")
rownames(tax) <- colnames(seqtab) #make sure to check that these match by hand first. In the future, set the Sequence IDs for uniquesToFasta to the sequences themselves.


sample_meta_sheet <- read.delim("../key.txt")

row.names(sample_meta_sheet) = paste0("BCS",sample_meta_sheet$Library,'-',sample_meta_sheet$Adapter.Order,'-',sample_meta_sheet$Primer.Tag.Number,'_',sample_meta_sheet$MLID)
sample_meta_sheet$Site.Code2 <- paste0(sample_meta_sheet$Site,"-",sample_meta_sheet$Subsite)
sample_meta_sheet$rootnum <- str_match(sample_meta_sheet$Subsite, "(?<=R)([0-9]{1,2})")[,1]
root_site_table <- data.frame(subsite <- c(replicate(5, "RA"), replicate(5, "RB"), replicate(5, "RC")), row.names = 1:15)
sample_meta_sheet$Subsite[which(!is.na(sample_meta_sheet$rootnum))] <- root_site_table[sample_meta_sheet$rootnum[which(!is.na(sample_meta_sheet$rootnum))],1]
sample_meta_sheet$Site.Code <- paste0(sample_meta_sheet$Site,"-",sample_meta_sheet$Subsite)

meta_subset <- sample_meta_sheet[rownames(seqtab),]

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE), sample_data(meta_subset), tax_table(tax))
rm(seqtab)
rm(tax)
saveRDS(ps, file = "phyloseq.RDS")

# Subset BCS3 for François ---------------------------------------
#if(FALSE) { #delete to "uncomment" 1

#BCS3.ps <- subset_samples(ps, Library == 3)
#BCS3.ps.f <- filter_taxa(BCS3.ps, function (x) sum(x) > 0, TRUE)
#saveRDS(BCS3.ps.f, file = "BCS3FM.rds")

# Sequencing Depth plot ---------------------------------------------------


Reads <- sort(sample_sums(ps))
depth_p <- plot_ly( type = 'bar', y = Reads, x = names(Reads), color = ~Reads, name = "Reads") %>%
  layout(yaxis = list(title = "Reads"), xaxis = list(title = 'Sample', showticklabels = FALSE), title = " Post-QC Sequencing Depth")
api_create(depth_p, filename = "bocas/seq_depth", sharing = "secret")


# Diversity Estimation ----------------------------------------------------


frequency_count_list <- build_frequency_count_tables(t(vegan_otu(otu_table(ps))))
objBayesList <- mclapply(frequency_count_list, objective_bayes_negbin, answers = T, mc.cores = nslots)
saveRDS(objBayesList, file = "objBayes.rds")
#uncomment above to generate new diversity estimates
#objBayesList <- readRDS(file = "objBayes.rds")
bayes.df <- as.data.frame(t(matrix(unlist(objBayesList), nrow = length (unlist(objBayesList[1])))))
colnames(bayes.df) <- names(unlist(objBayesList[[1]]))
rownames(bayes.df) <- names(frequency_count_list)

bayes_combined <- cbind(bayes.df[1:19], sample_meta_sheet[rownames(bayes.df),])

est.N_p <-plot_ly(data = bayes_combined, y = ~est, boxpoints = "suspectedoutliers") %>% add_boxplot(x = ~Habitat) %>% 
  layout(yaxis = list(title = "Median estimated number of ASVs"), margin = list(b =70))
api_create(est.N_p, filename = "bocas/alpha/ObjBayesEst_all", sharing = "secret")

est.Nxtype_p <-plot_ly(data = bayes_combined, x = ~est, y = ~Habitat) %>% add_boxplot(color = ~Sample.Type, jitter = 0.3) %>%
  layout(boxmode = "group", yaxis=list(title=""), xaxis = list(title="Median estimated Amplicon Sequence Variant count"), margin = list(l = 120))
api_create(est.Nxtype_p, filename = "bocas/alpha/ObjBayesEstxtype_all", sharing = "secret")




#breakaway(frequency_count_list[[2]])
#bayesian_results <- objective_bayes_negbin(frequency_count_list[[1]], answers = T)
#shannon_better(frequency_count_list[[2]])

#obs_shannon_all_plot <- plot_richness(ps, x="Habitat", measures=c("Observed", "Shannon"), color="Sample.Type") + theme_bw()
#ggsave("obs_shannon_all.pdf", obs_shannon_all_plot, width = 11, height = 8.5)



# Species pool estimation -------------------------------------------------

coral.tab <- vegan_otu(otu_table(subset_samples(ps, (Sample.Type != "eDNA") & (Habitat == "Agaricia") & (Sample.Type != "Sediment"))))
seagrass.tab <- vegan_otu(otu_table(subset_samples(ps, (Sample.Type != "eDNA") & (Habitat == "Seagrass") & (Sample.Type != "Sediment"))))
mangrove.tab <- vegan_otu(otu_table(subset_samples(ps, (Sample.Type != "eDNA") & (Habitat == "Mangrove root") & (Sample.Type != "Sediment"))))
sediment.tab <- vegan_otu(otu_table(subset_samples(ps, (Sample.Type == "Sediment"))))


coral.pool <- poolaccum(coral.tab)
coral.pool.df <- as.data.frame(coral.pool$means)
coral.pool.df$Habitat <- rep("Coral", nrow(coral.pool.df))
seagrass.pool <- poolaccum(seagrass.tab)
seagrass.pool.df <- as.data.frame(seagrass.pool$means)
seagrass.pool.df$Habitat <- rep("Seagrass", nrow(seagrass.pool.df))
mangrove.pool <- poolaccum(mangrove.tab)
mangrove.pool.df <- as.data.frame(mangrove.pool$means)
mangrove.pool.df$Habitat <- rep("Mangrove", nrow(mangrove.pool.df))
sediment.pool <- poolaccum(sediment.tab)
sediment.pool.df <- as.data.frame(sediment.pool$means)
sediment.pool.df$Habitat <- rep("Sediment", nrow(sediment.pool.df))

asvpool <- rbind(coral.pool.df, seagrass.pool.df, mangrove.pool.df, sediment.pool.df)
asvpool <- group_by(asvpool, Habitat)
pool_S.p <- plot_ly(data = asvpool, x = ~N, y = ~S, type = "scatter", mode = "lines", color = ~Habitat) %>%
  layout(title = "Observed Diversity", xaxis = list(rangemode="tozero", title = "# Samples"), yaxis = list(rangemode="tozero", title = "Sequence Variants"))
api_create(pool_S.p, filename = "bocas/pools/observed", sharing = "secret")

pool_Chao.p <- plot_ly(data = asvpool, x = ~N, y = ~Chao, type = "scatter", mode = "lines", color = ~Habitat) %>%
  layout(title = "Chao estimator", xaxis = list(rangemode="tozero", title = "# Samples"), yaxis = list(rangemode="tozero", title = "Sequence Variants"))
api_create(pool_Chao.p, filename = "bocas/pools/Chao", sharing = "secret")

pool_combined.p <- subplot(pool_S.p, pool_Chao.p, shareX=TRUE, shareY=TRUE)
api_create(pool_combined.p, filename = "bocas/pools/S-Chao", sharing = "secret")



# Taxonomic composition ---------------------------------------------------
#} #delete to "uncomment" 1

ps.tr <- transform_sample_counts(ps, function(x) x / sum(x) ) #transformed to relative abundances
by.phylum.tr <- tax_glom(ps.tr, taxrank='Phylum')
by.phylum.tr.f <- filter_taxa(by.phylum.tr, function (x) mean(x) > 5e-3, TRUE)

Agaricia.tr <- filter_taxa(subset_samples(ps.tr, ((Habitat == "Agaricia") & (Sample.Type != "eDNA"))), function (x) sum(x) > 0, TRUE)

GPS.coords <- as.data.frame(read_tsv(file = "GPS.txt"))
row.names(GPS.coords) <- GPS.coords$name
coral.500 <- subset_samples(Agaricia.tr, Sample.Type == "500 um")
coral500.asvtab <- vegan_otu(otu_table(coral.500)) #use relative abundances as way of normalizing data
coral500.dist <- vegdist(coral500.asvtab)
coral500.dist.hist <- plot_ly(x = as.vector(coral500.dist), type = "histogram", cumulative = list(enabled=TRUE), histnorm = "probability") %>%
  layout(title = "Cumulative Bray-Curtis Distances for Agaricia 500 um fraction")
api_create(coral500.dist.hist, filename = "bocas/coral500hist", sharing = "secret")


# Network maps with sliders ----------------------------


coral500_network_plot <- plotly_ps_network(coral.500, coords = GPS.coords)
api_create(coral500_network_plot, filename = "bocas/network/coral500", sharing = "secret")
merge_test_no.sed <- merge_samples(subset_samples(Agaricia.tr, (Sample.Type == "Sessile") | (Sample.Type == "500 um") | (Sample.Type == "100 um")), "Site.Code")
BCS1.merged.plot <- plotly_ps_network_by_samplenames(merge_test_no.sed, coords = GPS.coords)
api_create(BCS1.merged.plot, filename = "bocas/network/Agaricia_merged", sharing = "secret")
Mangrove.tr <- filter_taxa(subset_samples(ps.tr, ((Habitat == "Mangrove root") & (Sample.Type != "eDNA") & ((Sample.Type == "Sessile") | (Sample.Type == "500 um") | (Sample.Type == "100 um")))), function (x) sum(x) > 0, TRUE)

Seagrass.tr <- filter_taxa(subset_samples(ps.tr, ((Habitat == "Seagrass") & (Sample.Type != "eDNA") & ((Sample.Type == "Sessile") | (Sample.Type == "500 um") | (Sample.Type == "100 um")))), function (x) sum(x) > 0, TRUE)

Mangrove.merged.tr <- merge_samples(Mangrove.tr, "Site.Code")
Seagrass.merged.tr <- merge_samples(Seagrass.tr, "Site.Code")

Mangrove.merged.plot <- plotly_ps_network_by_samplenames(Mangrove.merged.tr, coords = GPS.coords)
Seagrass.merged.plot <- plotly_ps_network_by_samplenames(Seagrass.merged.tr, coords = GPS.coords)

api_create(Mangrove.merged.plot, filename = "bocas/network/Mangrove_merged", sharing = "secret")
api_create(Seagrass.merged.plot, filename = "bocas/network/Seagrass_merged", sharing = "secret")

#These are glommed by phylum even though the names don't say so
BCS1.tr.f <- subset_samples(by.phylum.tr.f, Library == 1)
BCS4.tr.f <- subset_samples(by.phylum.tr.f, Library == 4)
BCS8.tr.f <- subset_samples(by.phylum.tr.f, Library == 8)
BCS9.tr.f <- subset_samples(by.phylum.tr.f, Library == 9)
BCS6.tr.f <- subset_samples(by.phylum.tr.f, Library == 6)
BCS11.tr.f <- subset_samples(by.phylum.tr.f, Library == 11)
BCS10.tr.f <- subset_samples(by.phylum.tr.f, Library == 10)
BCS1.phy_div_plot <- plot_bar(BCS1.tr.f,'MLID', 'Abundance', fill='Phylum', labs(y='Relative abundance', title='BCS1 COI Phylum-level Diversity'), facet_grid = "Sample.Type~.")
BCS4.phy_div_plot <- plot_bar(BCS4.tr.f,'MLID', 'Abundance', fill='Phylum', labs(y='Relative abundance', title='BCS4 COI Phylum-level Diversity'), facet_grid = "Sample.Type~.")
BCS8.phy_div_plot <- plot_bar(BCS8.tr.f,'MLID', 'Abundance', fill='Phylum', labs(y='Relative abundance', title='BCS8 COI Phylum-level Diversity'), facet_grid = "Sample.Type~.")
BCS9.phy_div_plot <- plot_bar(BCS9.tr.f,'MLID', 'Abundance', fill='Phylum', labs(y='Relative abundance', title='BCS9 COI Phylum-level Diversity'), facet_grid = "Sample.Type~.")
BCS6.phy_div_plot <- plot_bar(BCS6.tr.f,'MLID', 'Abundance', fill='Phylum', labs(y='Relative abundance', title='BCS6 COI Phylum-level Diversity'), facet_grid = "Sample.Type~.")
BCS11.phy_div_plot <- plot_bar(BCS11.tr.f,'MLID', 'Abundance', fill='Phylum', labs(y='Relative abundance', title='BCS11 COI Phylum-level Diversity'), facet_grid = "Sample.Type~.")
BCS10.phy_div_plot <- plot_bar(BCS10.tr.f,'MLID', 'Abundance', fill='Phylum', labs(y='Relative abundance', title='BCS10 COI Phylum-level Diversity'), facet_grid = "Sample.Type~.")
BCS1.phy_div_plotly <- ggplotly(BCS1.phy_div_plot)
BCS4.phy_div_plotly <- ggplotly(BCS1.phy_div_plot)
BCS8.phy_div_plotly <- ggplotly(BCS1.phy_div_plot)
BCS9.phy_div_plotly <- ggplotly(BCS1.phy_div_plot)
BCS6.phy_div_plotly <- ggplotly(BCS6.phy_div_plot)
BCS11.phy_div_plotly <- ggplotly(BCS11.phy_div_plot)
BCS10.phy_div_plotly <- ggplotly(BCS10.phy_div_plot)
api_create(BCS1.phy_div_plotly, filename = "bocas/blca/BCS1/phylum_tax", sharing = "secret")
api_create(BCS4.phy_div_plotly, filename = "bocas/blca/BCS4/phylum_tax", sharing = "secret")
api_create(BCS8.phy_div_plotly, filename = "bocas/blca/BCS8/phylum_tax", sharing = "secret")
api_create(BCS9.phy_div_plotly, filename = "bocas/blca/BCS9/phylum_tax", sharing = "secret")
api_create(BCS10.phy_div_plotly, filename = "bocas/blca/BCS10/phylum_tax", sharing = "secret")
api_create(BCS11.phy_div_plotly, filename = "bocas/blca/BCS11/phylum_tax", sharing = "secret")
api_create(BCS6.phy_div_plotly, filename = "bocas/blca/BCS6/phylum_tax", sharing = "secret")


# Class-level ------------------------------------


BCS1.ps <- subset_samples(ps, Library == 1)
BCS1.class <- tax_glom(BCS1.ps, taxrank='Class')
#BCS1.class.tr <- transform_sample_counts(BCS1.class, function(x) x / sum(x) )
#BCS1.class.tr.f <- filter_taxa(BCS1.class.tr, function (x) mean(x) > 5e-3, TRUE)
BCS1.cla_div_plot <- plot_bar(BCS1.class,'MLID', 'Abundance', fill='Class', labs(y='Relative abundance', title='BCS1 COI Class-level Diversity'))
BCS_class.ly <- ggplotly(BCS1.cla_div_plot)
api_create(BCS_class.ly, filename = "bocas/blca/BCS1/class_tax", sharing = "secret")

BCS1.fam <- tax_glom(BCS1.ps, taxrank='Family')
BCS1.fam_div_plot <- plot_bar(BCS1.fam,'MLID', 'Abundance', fill='Family', labs(y='Relative abundance', title='BCS1 COI Family-level Diversity'))
BCS1_fami.ly <- ggplotly(BCS1.fam_div_plot)
api_create(BCS1_fami.ly, filename = "bocas/blca/BCS1/family_tax", sharing = "secret")


# Vegan ----------------------------------------


vegan.asvtab <- vegan_otu(otu_table(ps.tr)) #use relative abundances as way of normalizing data
bc.ord <- metaMDS(vegan.asvtab, distance = 'bray', k = 3, trymax = 800)
bc.df <- data.frame(bc.ord$points, meta_subset)
bc.ord.plot <- plot_ly(type = 'scatter3d', mode = 'markers', data = bc.df, x = ~MDS1, y = ~MDS2, z = ~MDS3, color = ~as.factor(Sample.Type), colors = "Set2", symbol = ~Habitat, text =~Site.Code, marker = list(size = 8, opacity = 0.4), symbols = c('circle', 'circle-open','cross', 'square-open', 'diamond', 'square' , 'diamond-open', 'circle-open', 'diamond') ) %>%
  layout(title= "BCS Bray-Curtis NMDS")
api_create(bc.ord.plot, filename = "bocas/BCOrd", sharing = "secret")
saveRDS(bc.df, file = "bc.ordination.df.rds")

# DESeq2 variance stabilization ----------------------------

# DESeq2 not working on Hydra. I can't get the DESeq2 package installed properly.


#ps.deseq <- phyloseq_to_deseq2(ps, ~ 1)
#geoMeans = apply(counts(ps.deseq), 1, gm_mean)
# You must step through the size factor and dispersion estimates prior to calling the getVarianceStabilizedData() function.
#ps.deseq = estimateSizeFactors(ps.deseq, geoMeans = geoMeans)
#ps.deseq = estimateDispersions(ps.deseq)
#ps.vst = getVarianceStabilizedData(ps.deseq)
#saveRDS(ps.vst, file = "ps.vst.rds")
#uncomment above to generate variance stabilization data object

#load rds if existing
#ps.vst <- readRDS(file = "ps.vst.rds")
#ps.vst[ps.vst < 0] <- 0
#pst.vs.nc.nonneg <- phyloseq(otu_table(ps.vst, taxa_are_rows = TRUE), sample_data(meta_subset), tax_table(tax))




