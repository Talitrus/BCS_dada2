library(phyloseq)
library(vegan)
library(ggplot2)
library(DESeq2)
library(parallel)
library(plotly)
require(breakaway)
library(RDPutils)


setwd("/groups/cbi/bryan/BCS_all/dada2_R/")
tax <- make_tax_table(in_file="BCS_RDP_output.txt", confidence = 0.5)
seqtab <- readRDS("seqtab_final.rds")
rownames(tax) <- colnames(seqtab) #make sure to check that these match by hand first. In the future, set the Sequence IDs for uniquesToFasta to the sequences themselves.
vegan_otu <- function(physeq) { #convert phyloseq OTU table into vegan OTU matrix
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


sample_meta_sheet <- read.delim("../key.txt")
sample_meta_sheet$Site.Code <- paste0(sample_meta_sheet$Site,"-",sample_meta_sheet$Subsite)
row.names(sample_meta_sheet) = paste0("BCS",sample_meta_sheet$Library,'-',sample_meta_sheet$Adapter.Order,'-',sample_meta_sheet$Primer.Tag.Number,'_',sample_meta_sheet$MLID)

meta_subset <- sample_meta_sheet[rownames(seqtab),]

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE), sample_data(meta_subset), tax_table(tax))
rm(seqtab)
rm(tax)
saveRDS(ps, file = "phyloseq.RDS")

# Subset BCS3 for François ---------------------------------------
if(FALSE) { #delete to "uncomment" 1

BCS3.ps <- subset_samples(ps, Library == 3)
BCS3.ps.f <- filter_taxa(BCS3.ps, function (x) sum(x) > 0, TRUE)
saveRDS(BCS3.ps.f, file = "BCS3FM.rds")

# Sequencing Depth plot ---------------------------------------------------


Reads <- sort(sample_sums(ps))
depth_p <- plot_ly( type = 'bar', y = Reads, x = names(Reads), color = ~Reads, name = "Reads") %>%
  layout(yaxis = list(title = "Reads"), xaxis = list(title = 'Sample', showticklabels = FALSE), title = " Post-QC Sequencing Depth")
api_create(depth_p, filename = "bocas/seq_depth", sharing = "secret")


# Diversity Estimation ----------------------------------------------------


frequency_count_list <- build_frequency_count_tables(t(vegan_otu(otu_table(ps))))
#objBayesList <- mclapply(frequency_count_list, objective_bayes_negbin, answers = T)
#saveRDS(objBayesList, file = "objBayes.rds")
#uncomment above to generate new diversity estimates
objBayesList <- readRDS(file = "objBayes.rds")
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

} #delete to "uncomment" 1


#breakaway(frequency_count_list[[2]])
#bayesian_results <- objective_bayes_negbin(frequency_count_list[[1]], answers = T)
#shannon_better(frequency_count_list[[2]])

#obs_shannon_all_plot <- plot_richness(ps, x="Habitat", measures=c("Observed", "Shannon"), color="Sample.Type") + theme_bw()
#ggsave("obs_shannon_all.pdf", obs_shannon_all_plot, width = 11, height = 8.5)



# Species pool estimation -------------------------------------------------

coral.tab <- vegan_otu(otu_table(subset_samples(ps, (Library != 10) & (Habitat == "Agaricia") & (Sample.Type != "Sediment"))))
seagrass.tab <- vegan_otu(otu_table(subset_samples(ps, (Library != 10) & (Habitat == "Seagrass") & (Sample.Type != "Sediment"))))
mangrove.tab <- vegan_otu(otu_table(subset_samples(ps, (Library != 10) & (Habitat == "Mangrove root") & (Sample.Type != "Sediment"))))
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

ps.tr <- transform_sample_counts(ps, function(x) x / sum(x) ) #transformed to relative abundances
by.phylum.tr <- tax_glom(ps.tr, taxrank='Phylum')
by.phylum.tr.f <- filter_taxa(by.phylum.tr, function (x) mean(x) > 5e-3, TRUE)


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
api_create(BCS1.phy_div_plotly, filename = "bocas/midori/BCS1/phylum_tax", sharing = "secret")
api_create(BCS4.phy_div_plotly, filename = "bocas/midori/BCS4/phylum_tax", sharing = "secret")
api_create(BCS8.phy_div_plotly, filename = "bocas/midori/BCS8/phylum_tax", sharing = "secret")
api_create(BCS9.phy_div_plotly, filename = "bocas/midori/BCS9/phylum_tax", sharing = "secret")
api_create(BCS10.phy_div_plotly, filename = "bocas/midori/BCS10/phylum_tax", sharing = "secret")
api_create(BCS11.phy_div_plotly, filename = "bocas/midori/BCS11/phylum_tax", sharing = "secret")
api_create(BCS6.phy_div_plotly, filename = "bocas/midori/BCS6/phylum_tax", sharing = "secret")


# Class-level ------------------------------------


BCS1.ps <- subset_samples(ps, Library == 1)
BCS1.class <- tax_glom(BCS1.ps, taxrank='Class')
BCS1.class.tr <- transform_sample_counts(BCS1.class, function(x) x / sum(x) )
BCS1.class.tr.f <- filter_taxa(BCS1.class.tr, function (x) mean(x) > 5e-3, TRUE)
BCS1.cla_div_plot <- plot_bar(BCS1.class.tr.f,'MLID', 'Abundance', fill='Class', labs(y='Relative abundance', title='BCS1 COI Class-level Diversity'))
BCS_class.ly <- ggplotly(BCS1.cla_div_plot)
api_create(BCS_class.ly, filename = "bocas/midori/BCS1/class_tax", sharing = "secret")


# Vegan ----------------------------------------


#vegan.asvtab <- vegan_otu(otu_table(ps.tr)) #use relative abundances as way of normalizing data
#bc.ord <- metaMDS(vegan.asvtab, distance = 'bray', k = 3, trymax = 800)
#bc.df <- data.frame(bc.ord$points, meta_subset)
#bc.ord.plot <- plot_ly(type = 'scatter3d', mode = 'markers', data = bc.df, x = ~MDS1, y = ~MDS2, z = ~MDS3, color = ~as.factor(Sample.Type), colors = "Set2", symbol = ~Habitat, text =~Site.Code, marker = list(size = 8, opacity = 0.4), symbols = c('circle', 'circle-open','cross', 'square-open', 'diamond', 'square' , 'diamond-open', 'circle-open', 'diamond') ) %>%
#  layout(title= "BCS 1,3,4,6,8,9,10,11 Bray-Curtis NMDS")
#api_create(bc.ord.plot, filename = "bocas/BCOrd", sharing = "secret")
#saveRDS(bc.df, file = "bc.ordination.df.rds")

# DESeq2 variance stabilization ----------------------------
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}



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




