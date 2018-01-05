library(phyloseq)
library(vegan)
library(ggplot2)
library(DESeq2)

setwd("/groups/cbi/bryan/BCS_all/dada2_R")
tax <- readRDS("tax_final.rds")
colnames(tax)[1] <- "Domain" #manual correction
seqtab <- readRDS("seqtab_final.rds")

vegan_otu <- function(physeq) { #convert phyloseq OTU table into vegan OTU matrix
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

samples.out <- rownames(seqtab)
MLIDs <- substr(samples.out, nchar(samples.out)-5, nchar(samples.out))

sample_meta_sheet <- read.delim("../key.txt")
sample_meta_sheet$Site.Code <- paste0(sample_meta_sheet$Site,"-",sample_meta_sheet$Subsite)
row.names(sample_meta_sheet) = paste0("BCS",sample_meta_sheet$Library,'-',sample_meta_sheet$Adapter.Order,'-',sample_meta_sheet$Primer.Tag.Number,'_',sample_meta_sheet$MLID)

meta_subset <- subset(sample_meta_sheet, Library == 1 | Library == 4 | Library == 8 | Library == 9)

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE), sample_data(meta_subset), tax_table(tax))
seq_depth <- sample_sums(ps)


obs_shannon_all_plot <- plot_richness(ps, x="Habitat", measures=c("Observed", "Shannon"), color="Sample.Type") + theme_bw()
ggsave("obs_shannon_all.pdf", obs_shannon_all_plot, width = 11, height = 8.5)
ps.tr <- transform_sample_counts(ps, function(x) x / sum(x) ) #transformed to relative abundances
by.phylum.tr <- tax_glom(ps.tr, taxrank='Phylum')
by.phylum.tr.f <- filter_taxa(by.phylum.tr, function (x) mean(x) > 5e-3, TRUE)

#phy_div_plot <- plot_bar(by.phylum.tr.f,'MLID', 'Abundance', fill='Phylum', labs(y='Relative abundance', title='BCS COI Phylum-level Diversity'))
#ggsave("phy_div_plot.png", phy_div_plot, width = 80, height = 5, limitsize = FALSE)


#despite the names, these are glommed by phylum
BCS1.tr.f <- subset_samples(by.phylum.tr.f, Library == 1)
BCS4.tr.f <- subset_samples(by.phylum.tr.f, Library == 4)
BCS8.tr.f <- subset_samples(by.phylum.tr.f, Library == 8)
BCS9.tr.f <- subset_samples(by.phylum.tr.f, Library == 9)
BCS1.phy_div_plot <- plot_bar(BCS1.tr.f,'MLID', 'Abundance', fill='Phylum', labs(y='Relative abundance', title='BCS1 COI Phylum-level Diversity'), facet_grid = "Sample.Type~.")
BCS4.phy_div_plot <- plot_bar(BCS4.tr.f,'MLID', 'Abundance', fill='Phylum', labs(y='Relative abundance', title='BCS4 COI Phylum-level Diversity'), facet_grid = "Sample.Type~.")
BCS8.phy_div_plot <- plot_bar(BCS8.tr.f,'MLID', 'Abundance', fill='Phylum', labs(y='Relative abundance', title='BCS8 COI Phylum-level Diversity'), facet_grid = "Sample.Type~.")
BCS9.phy_div_plot <- plot_bar(BCS9.tr.f,'MLID', 'Abundance', fill='Phylum', labs(y='Relative abundance', title='BCS9 COI Phylum-level Diversity'), facet_grid = "Sample.Type~.")


bray_dist <- phyloseq::distance(ps, method = "bray")
bray_MDS <- ordinate(ps, method = "MDS", distance = "bray")
bray_plot <- plot_ordination(ps, bray_MDS, color="Habitat",shape="Sample.Type")

# Vegan
vegan.asvtab <- vegan_otu(otu_table(ps))
pool <- specpool(vegan.asvtab, sample_data(ps)$Sample.Type)


# DESeq2 variance stabilization
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

nonnegative <- function (x) {
  return(max(x, 0))
}
Vnn <- Vectorize(nonnegative)


ps.deseq <- phyloseq_to_deseq2(ps, ~ 1)
geoMeans = apply(counts(ps.deseq), 1, gm_mean)
# You must step through the size factor and dispersion estimates prior to calling the getVarianceStabilizedData() function.
ps.deseq = estimateSizeFactors(ps.deseq, geoMeans = geoMeans)
ps.deseq = estimateDispersions(ps.deseq)
ps.vst = getVarianceStabilizedData(ps.deseq)
dim(ps.vst)
