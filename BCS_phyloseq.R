library(phyloseq)
library(vegan)
library(ggplot2)
library(DESeq2)
library(parallel)
library(plotly)
require(breakaway)



setwd("/groups/cbi/bryan/BCS_all/dada2_R/")
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

sample_meta_sheet <- read.delim("../key.txt")
sample_meta_sheet$Site.Code <- paste0(sample_meta_sheet$Site,"-",sample_meta_sheet$Subsite)
row.names(sample_meta_sheet) = paste0("BCS",sample_meta_sheet$Library,'-',sample_meta_sheet$Adapter.Order,'-',sample_meta_sheet$Primer.Tag.Number,'_',sample_meta_sheet$MLID)

#meta_subset <- subset(sample_meta_sheet, Library == 1 | Library == 4 | Library == 8 | Library == 9)

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE), sample_data(sample_meta_sheet), tax_table(tax))


seq_depth <- sample_sums(ps)

frequency_count_list <- build_frequency_count_tables(t(vegan_otu(otu_table(ps))))
objBayesList <- mclapply(frequency_count_list, objective_bayes_negbin, answers = T)
saveRDS(objBayesList, file = "objBayes.rds")
#uncomment above to generate new diversity estimates

#breakaway(frequency_count_list[[2]])
#bayesian_results <- objective_bayes_negbin(frequency_count_list[[1]], answers = T)
#shannon_better(frequency_count_list[[2]])

#obs_shannon_all_plot <- plot_richness(ps, x="Habitat", measures=c("Observed", "Shannon"), color="Sample.Type") + theme_bw()
#ggsave("obs_shannon_all.pdf", obs_shannon_all_plot, width = 11, height = 8.5)
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
BCS4.phy_div_plotly <- ggplotly(BCS4.phy_div_plot)
BCS8.phy_div_plotly <- ggplotly(BCS8.phy_div_plot)
BCS9.phy_div_plotly <- ggplotly(BCS9.phy_div_plot)
BCS6.phy_div_plotly <- ggplotly(BCS6.phy_div_plot)
BCS11.phy_div_plotly <- ggplotly(BCS11.phy_div_plot)
BCS10.phy_div_plotly <- ggplotly(BCS10.phy_div_plot)
api_create(BCS1.phy_div_plotly, filename = "bocas/BCS1/phylum_tax", sharing = "secret")
api_create(BCS4.phy_div_plotly, filename = "bocas/BCS4/phylum_tax", sharing = "secret")
api_create(BCS8.phy_div_plotly, filename = "bocas/BCS8/phylum_tax", sharing = "secret")
api_create(BCS9.phy_div_plotly, filename = "bocas/BCS9/phylum_tax", sharing = "secret")
api_create(BCS10.phy_div_plotly, filename = "bocas/BCS10/phylum_tax", sharing = "secret")
api_create(BCS11.phy_div_plotly, filename = "bocas/BCS11/phylum_tax", sharing = "secret")
api_create(BCS6.phy_div_plotly, filename = "bocas/BCS6/phylum_tax", sharing = "secret")

bray_dist <- phyloseq::distance(ps, method = "bray")
bray_MDS <- ordinate(ps, method = "MDS", distance = "bray")
bray_plot <- plot_ordination(ps, bray_MDS, color="Habitat",shape="Sample.Type")

# Vegan
vegan.asvtab <- vegan_otu(otu_table(ps))
bc.ord <- metaMDS(vegan.asvtab, distance = 'bray', k = 3, trymax = 800)
bc.df <- data.frame(bc.ord$points, sample_meta_sheet)
bc.ord.plot <- plot_ly(type = 'scatter3d', mode = 'markers', data = bc.df, x = ~MDS1, y = ~MDS2, z = ~MDS3, color = ~Habitat, shape = ~Sample.Type, text =~Site.Code)
api_create(bc.ord.plot, filename = "bocas/BCOrd", sharing = "secret")
#pool <- specpool(vegan.asvtab, sample_data(ps)$Sample.Type)
saveRDS(bc.df, file = "bc.ordination.df.rds")

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
saveRDS(ps.vst, file = "ps.vst.rds")
#uncomment above to generate variance stabilization data object

#load rds if existing
#ps.vst <- readRDS(file = "ps.vst.rds")
ps.vst[ps.vst < 0] <- 0
pst.vs.nc.nonneg <- phyloseq(otu_table(ps.vst, taxa_are_rows = TRUE), sample_data(sample_meta_sheet), tax_table(tax))



#objBayesList <- readRDS(file = "objBayes.rds")
bayes.df <- as.data.frame(t(matrix(unlist(objBayesList), nrow = length (unlist(objBayesList[1])))))
colnames(bayes.df) <- c(rownames(objBayesList[[1]]$results), rownames(objBayesList[[1]]$fits))
rownames(bayes.df) <- names(frequency_count_list)

bayes_combined <- cbind(bayes.df[1:19], sample_meta_sheet[rownames(bayes.df),])

median.C_p <-plot_ly(data = bayes_combined, y = ~median.C, boxpoints = "suspectedoutliers") %>% add_boxplot(x = ~Habitat) %>% 
  layout(yaxis = list(title = "Median estimated number of ASVs"))
api_create(median.C_p, filename = "bocas/medianC_all", sharing = "secret")

median.Cxtype_p <-plot_ly(data = bayes_combined, x = ~median.C, y = ~Habitat) %>% add_boxplot(color = ~Sample.Type, jitter = 0.3) %>%
  layout(boxmode = "group", yaxis=list(title=""), xaxis = list(title="Median estimated Exact Sequence Variant count"), margin = list(l = 90))
median.Cxtype_p
api_create(median.Cxtype_p, filename = "bocas/medianCxtype_all", sharing = "secret")
