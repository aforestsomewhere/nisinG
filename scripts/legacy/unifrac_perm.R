#############################################
### Fig 2 UNIFRAC weighted         ###
#############################################

#flush previous plots/objects
dev.off(dev.list()["RStudioGD"])

####################
### Functions    ###
####################

metaphlanToPhyloseq = function(
  tax,
  metadat=NULL,
  simplenames=TRUE,
  roundtointeger=FALSE,
  split="|"){
  xnames = rownames(tax)
  shortnames = gsub(paste0(".+\\", split), "", xnames)
  if(simplenames){
    rownames(tax) = shortnames
  }
  if(roundtointeger){
    tax = round(tax * 1e4)
  }
  x2 = strsplit(xnames, split=split, fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] = x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  if(is.null(metadat)){
    res = phyloseq::phyloseq(taxmat, otutab)
  }else{
    res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
  }
  return(res)
}

####################
### read data    ###
####################

#Read cross reference of SGB to taxonomy
#remove 'SGB' and '_group' from first column to allow matching to the tree
mpa_legend <- read.table("data\\mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt", sep = '\t')
mpa_legend$V1 <- gsub("_group", "", mpa_legend$V1)
mpa_legend$V1 <- gsub("SGB", "", mpa_legend$V1)

#read metaphlan.tsv
meta <- read.csv("cleaned_data/metadata.csv", nrows=36, row.names=1)
mpa_table <- read.table("cleaned_data/mpa_table.tsv", comment.char = '#', sep = '\t', header = TRUE, row.names=NULL)
#filter to species level assignments
mpa_table <- mpa_table[grep('s__',mpa_table[,1]),]
mpa_table[,1] <- gsub(".+\\|s__", "", mpa_table[,1])
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]

#########################
### form Phyloseq obj ###
#########################

#generate base obj
mpa.phy = metaphlanToPhyloseq(mpa_table)
phy <- metaphlanToPhyloseq(
  phyloseq::otu_table(mpa.phy),
  metadat=phyloseq::sample_data(meta, errorIfNULL = TRUE),
  simplenames=TRUE,
  roundtointeger=FALSE,
  split="|")

#generate guide tree, importing full tree from MPA and filtering by taxa detected
read_tree("data\\mpa_vJan21_CHOCOPhlAnSGB_202103.nwk") -> mpa_tree
new_tree <- mpa_tree; new_tree$tip.label <- mpa_legend$V2[match(mpa_tree$tip.label,mpa_legend$V1)]
new_tree$tip.label <- gsub(".+\\|s__", "", new_tree$tip.label)
filt_tree <- ape::keep.tip(new_tree, intersect(rownames(phy@otu_table),new_tree$tip.label))
filt_mpa <- mpa_table[filt_tree$tip.label,] / 100.0
phy = phyloseq::merge_phyloseq(phy, filt_tree)

#########################
### Normalisation     ###
#########################


# #log10 transform
# phy.log10 <- phy
# phy.log10@otu_table <- log10(1 + phy@otu_table)
# 
# #Cumulative sum scaling (CSS) normalization of OTU abundance table.
# phy.cum <- phy
# phy.cum <- metagMisc::phyloseq_transform_css(phy.cum)
# 
# filt_mpa_table_log10 <- log10(1 + filt_mpa)
# filt_mpa_table_sinh_sqrt <- asinh(sqrt(filt_mpa))


#########################
### UNIFRAC Phyloseq  ###
#########################
#'There is no % of variance associated with each axis in nMDS in contrast with other 
#'Principal Component Methods like PCA, CA, PCoA (= MDS).
#'Contrary to PCA, PCoA, or CA, which are eigenvector-based methods, nMDS calculations 
#'do not maximize the variability associated with individual axes of the ordination

#Weighted, unnormalised
phy@sam_data$timepoint = factor(phy@sam_data$timepoint, levels = c("T_0", "T_6", "T_24"))
wUF.ordu = ordinate(phy, method="NMDS", distance="unifrac", weighted=TRUE)
p0 <- plot_ordination(phy, wUF.ordu, type="sites", color="treatment", shape="timepoint") +
  scale_color_manual(name="Treatment",
                     labels=c('Control', 
                              expression(paste(italic("F. nucleatum"))), 
                              expression(paste(italic("F. nucleatum/S. salivarius"))), 
                              expression(paste(italic('S. salivarius')))),
                     values = wes_palette("Darjeeling1", n = 4, type = "discrete")) +
  scale_fill_manual(name="Treatment",
                    labels=c('Control', 
                             expression(paste(italic("F. nucleatum"))), 
                             expression(paste(italic("F. nucleatum/S. salivarius"))), 
                             expression(paste(italic('S. salivarius')))),
                    values = wes_palette("Darjeeling1", n = 4, type = "discrete")) +
  labs(shape="Timepoint") +
  geom_point(size=4) +
  theme(axis.title.x= element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=18, face = "bold"),
        legend.text = element_text(size=16),
        # legend.position="bottom", legend.box = "horizontal",
        legend.position = c(0.3, 0.85), legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid")) +
  guides(col = guide_legend(order = 1),
         shape = guide_legend(order = 0),
         fill = "none") +
  ggtitle("Weighted UniFrac, nMDS, unnormalised, MPA4")
p0
ggsave("figures\\UNIFRAC_nMDS_unnorm.emf", width = 10, height = 10, dpi = 600, device = {function(filename, ...) devEMF::emf(file = filename, ...)})

#Unweighted, unnormalised
phy@sam_data$timepoint = factor(phy@sam_data$timepoint, levels = c("T_0", "T_6", "T_24"))
uUF.ordu = ordinate(phy, method="NMDS", distance="unifrac", weighted=FALSE)
p0 <- plot_ordination(phy, uUF.ordu, type="sites", color="timepoint") +
  scale_color_manual(name="Treatment",
                     labels=c('Control', 
                              expression(paste(italic("F. nucleatum"))), 
                              expression(paste(italic("F. nucleatum/S. salivarius"))), 
                              expression(paste(italic('S. salivarius')))),
                     values = wes_palette("Darjeeling1", n = 4, type = "discrete")) +
  scale_fill_manual(name="Treatment",
                    labels=c('Control', 
                             expression(paste(italic("F. nucleatum"))), 
                             expression(paste(italic("F. nucleatum/S. salivarius"))), 
                             expression(paste(italic('S. salivarius')))),
                    values = wes_palette("Darjeeling1", n = 4, type = "discrete")) +
  labs(shape="Timepoint") +
  theme(axis.title.x= element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=18, face = "bold"),
        legend.text = element_text(size=16),
        # legend.position="bottom", legend.box = "horizontal",
        legend.position = c(0.26, 0.85), legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid")) +
  guides(col = guide_legend(order = 1),
         shape = guide_legend(order = 0),
         fill = "none") +
  geom_mark_ellipse(aes(fill=timepoint),
                    alpha=0.1,
                    linetype = 2) +
  geom_point(aes(shape=treatment), size=4) +
  ggtitle("Unweighted UniFrac, nMDS, unnormalised, MPA4")
p0


Wunifrac.dist <- UniFrac(phy, 
                        weighted = TRUE, 
                        normalized = FALSE,  
                        parallel = FALSE, 
                        fast = TRUE)

Uunifrac.dist <- UniFrac(phy, 
                        weighted = FALSE, 
                        normalized = FALSE,  
                        parallel = FALSE, 
                        fast = TRUE)
Wpermanova <- adonis2(Wunifrac.dist ~ treatment*timepoint, data = meta)
Wpermanova <- adonis2(Wunifrac.dist ~ timepoint, data = meta)
Upermanova <- adonis2(Uunifrac.dist ~ treatment, data = meta)
Upermanova <- adonis2(Uunifrac.dist ~ timepoint, data = meta)
ggsave("figures\\UNIFRAC_nMDS_unnorm_unweighted.emf", width = 10, height = 10, dpi = 600, device = {function(filename, ...) devEMF::emf(file = filename, ...)})

#p_timepoint <-plot_grid(p0t,p1t,p2t,p3t)
#p_timepoint
