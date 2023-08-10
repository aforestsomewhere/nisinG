#' Script to assess beta diversity for NisinG samples
#' Removed T_0 due to large distances
#' Developed by A. Kate Falà (https://github.com/aforestsomewhere)
#' Computes Bray-Curtis, and UNIFRAC (weighted and unweighted)
#
###############################

##############################
### Packages and palettes  ###
##############################
library(tidyverse)
library(dplyr)
library(vegan)
library(cowplot)
library(devEMF)
library(webshot)
library(ggplot2)
library(rbiom)
library(phyloseq)
library(ggrepel)
library(microbiome)
library(pairwiseAdonis)
library(ggforce)
theme_set(theme_bw())

####################
### functions    ###
####################

sigcode = function(P){
  if(P<0){ 
    return("NaN") 
  }
  else if(P<0.001){ 
    return("***") 
  }
  else if(P<0.01){ 
    return("**") 
  }
  else if(P<=0.05){ 
    return("*") 
  }
  else if(P<0.1){ 
    return(".") 
  }
  else{ return(" ")}
}

####################
### read data    ###
####################
#filtering out T_0 samples
meta <- read.csv("cleaned_data/metadata.csv", nrows=36, row.names=1)
meta <- meta %>% filter(timepoint!="T_0")
mpa_table <- read.table("cleaned_data/mpa_table.tsv", comment.char = '#', sep = '\t', header = TRUE, row.names=NULL)

#filter to species level assignments
mpa_table <- mpa_table[grep('s__',mpa_table[,1]),]
mpa_table[,1] <- gsub(".+\\|s__", "", mpa_table[,1])
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]
#filter abundance data to timepoint T_24 only
samples <- rownames(meta)
mpa_table[,colnames(mpa_table) %in% samples] -> mpa_table

#generate matrix of abundance data
mpa_mat = t(as.matrix(mpa_table))
#order samples (rows) to properly merge with infO
mpa_mat <- mpa_mat[sort(rownames(mpa_mat)),]

#Phyloseq for UNIFRAC
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
#generate base Phyloseq obj
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



################
### ANOSIM   ###
################
#ANOSIM (ANalysis Of Similarities) is a non-parametric test of significant difference 
#between two or more groups, based on any distance measure (Clarke 1993). The distances 
#are converted to ranks. ANOSIM is normally used for taxa-in-samples data, where groups 
#of samples are to be compared.

dist_mat_bc = vegdist(mpa_mat, method="bray")
#UNIFRAC
#Weighted, unnormalised
phy@sam_data$timepoint = factor(phy@sam_data$timepoint, levels = c("T_6", "T_24"))
dist_mat_wUF <- UniFrac(phy, 
                        weighted = TRUE, 
                        normalized = FALSE,  
                        parallel = FALSE, 
                        fast = TRUE)
dist_mat_uUF <- UniFrac(phy, 
                        weighted = FALSE, 
                        normalized = FALSE,  
                        parallel =FALSE, 
                        fast = TRUE)
#initialise output DF
anosim_df_all <- data.frame(Data="MetaPhlAn4", 
                            Metric="Bray-Curtis dissimilarity",
                            Comparing=NA,
                            Test="ANOSIM",
                            Test_Statistic=NA,
                            P_value=NA,
                            Significance_codes="NA")
#timepoint
calc_anosim <- function(distmat, group, dataset, metric, groupname){
  res <- anosim(distmat,group, permutations = 9999)
  P <- res$signif
  F<-res$statistic
  S<-sigcode(P)
  anosim_df <- c(dataset, metric, groupname, "anosim",F,P,S)
  return(anosim_df)
}
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_bc,group=meta$timepoint, dataset="MetaPhlAn4", metric="Bray-Curtis dissimilarity", groupname="Timepoint"))
anosim_df_all[-1,]-> anosim_df_all #remove first shaper row (empty)
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_bc,group=meta$treatment, dataset="MetaPhlAn4", metric="Bray-Curtis dissimilarity", groupname="Treatment"))
#Weighted UNIFRAC
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_wUF,group=meta$timepoint, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_wUF,group=meta$treatment, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Treatment"))
#Unweighted UNIFRAC
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_uUF,group=meta$timepoint, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_uUF,group=meta$treatment, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Treatment"))

#adonis
calc_adonis <- function(distmat, group, dataset, metric, groupname){
  res <- adonis2(distmat~group, permutations = 9999)
  P<- res["group","Pr(>F)"]
  F<-res["group","F"]
  S<-sigcode(P)
  adonis_df <- c(dataset, metric, groupname, "adonis2",F,P,S)
  return(adonis_df)
}
#BC
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_bc,group=meta$timepoint, dataset="MetaPhlAn4", metric="Bray-Curtis dissimilarity", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_bc,group=meta$treatment, dataset="MetaPhlAn4", metric="Bray-Curtis dissimilarity", groupname="Treatment"))
#Weighted UNIFRAC
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_wUF,group=meta$timepoint, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_wUF,group=meta$treatment, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Treatment"))
#Unweighted UNIFRAC
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_uUF,group=meta$timepoint, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_uUF,group=meta$treatment, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Treatment"))

#beta dispersion
#checkk beta dispersions within groups - TREATMENT
#adonis
calc_disp <- function(distmat, group, dataset, metric, groupname){
  res <- betadisper(distmat, group)
  disp <- permutest(res, pairwise=TRUE, permutations=1000)
  R<-disp$tab$F[1]
  P<- disp$tab$'Pr(>F)'[1]
  S<-sigcode(P)
  adonis_df <- c(dataset, metric, groupname, "beta_dispersion",R,P,S)
  return(adonis_df)
}

#calc_disp(dist_mat_uUF,group=meta$timepoint, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Timepoint")
#BC
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_bc,group=meta$timepoint, dataset="MetaPhlAn4", metric="Bray-Curtis dissimilarity", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_bc,group=meta$treatment, dataset="MetaPhlAn4", metric="Bray-Curtis dissimilarity", groupname="Treatment"))
#Weighted UNIFRAC
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_wUF,group=meta$timepoint, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_wUF,group=meta$treatment, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Treatment"))
#Unweighted UNIFRAC
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_uUF,group=meta$timepoint, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_uUF,group=meta$treatment, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Treatment"))


write.csv(anosim_df_all, file="reports\\beta_testing_T6_T24.csv")
#############
# Plots     #
#############

#unweighted
#phy@sam_data$timepoint = factor(phy@sam_data$timepoint, levels = c("T_0", "T_6", "T_24"))
uUF.ordu = ordinate(phy, method="NMDS", distance="unifrac", weighted=FALSE)

#base plot for legend
leggo <- plot_ordination(phy, uUF.ordu, type="sites", color="timepoint") +
  scale_color_manual(name="Timepoint",
                     labels=c("T_6",
                              "T_24"),
                     values = wes_palette("Darjeeling1", n = 3, type = "discrete")) +
  scale_fill_manual(name="Timepoint",
                    labels=c("T_6",
                             "T_24"),
                    values = wes_palette("Darjeeling1", n = 3, type = "discrete")) +
  labs(shape="Treatment") +
  theme(axis.title.x= element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=18, face = "bold"),
        legend.text = element_text(size=16),
        legend.position="bottom", legend.box = "vertical",
        #legend.position = c(0.26, 0.85), legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid")) +
  guides(shape = guide_legend(nrow=4, byrow=TRUE),
         fill = "none") +
  geom_mark_ellipse(aes(fill=timepoint, label=timepoint),
                    alpha=0.1,
                    linetype = 2) +
  geom_point(aes(shape=treatment),size=3) +
  scale_shape_manual(name="Treatment",
                     labels=c('Control', 
                              expression(paste(italic("F. nucleatum"))), 
                              expression(paste(italic("F. nucleatum/S. salivarius"))), 
                              expression(paste(italic('S. salivarius')))),
                     values = c(15,16,17,18)) +
  ggtitle("Unweighted UniFrac, nMDS, unnormalised, MPA4")
plegend <- get_legend(leggo)
#unweighted UNIFRAC
p0 <- plot_ordination(phy, uUF.ordu, type="sites", color="timepoint") +
  scale_color_manual(name="Timepoint",
                     labels=c("T_6",
                              "T_24"),
                     values = wes_palette("Darjeeling1", n = 3, type = "discrete")) +
  scale_fill_manual(name="Timepoint",
                    labels=c("T_6",
                             "T_24"),
                    values = wes_palette("Darjeeling1", n = 3, type = "discrete")) +
  labs(shape="Treatment") +
  theme(axis.title.x= element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=16),
        legend.position="bottom", legend.box = "horizontal",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid")) +
  guides(col = "none",
         shape = "none",
         fill = "none") +
  geom_mark_ellipse(aes(fill=timepoint, label=timepoint),
                    alpha=0.1,
                    linetype = 2) +
  geom_point(aes(shape=treatment),size=3, alpha=1) +
  scale_shape_manual(name="Treatment",
                     labels=c('Control', 
                              expression(paste(italic("F. nucleatum"))), 
                              expression(paste(italic("F. nucleatum/S. salivarius"))), 
                              expression(paste(italic('S. salivarius')))),
                     values = c(15,16,17,18)) +
  ggtitle("nMDS of unweighted UNIFRAC distances")
p0 <- p0 + expand_limits(x = c(-.4,.4), y = c(-.2, .2))
p0
#weighted
wUF.ordu = ordinate(phy, method="NMDS", distance="unifrac", weighted=TRUE)
p1 <- plot_ordination(phy, wUF.ordu, type="sites", color="timepoint") +
  scale_color_manual(name="Timepoint",
                     labels=c("T_6",
                              "T_24"),
                     values = wes_palette("Darjeeling1", n = 3, type = "discrete")) +
  scale_fill_manual(name="Timepoint",
                    labels=c("T_6",
                             "T_24"),
                    values = wes_palette("Darjeeling1", n = 3, type = "discrete")) +
  labs(shape="Treatment") +
  theme(axis.title.x= element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=16),
        legend.position="bottom", legend.box = "horizontal",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid")) +
  guides(col = "none",
         shape = "none",
         fill = "none") +
  geom_mark_ellipse(aes(fill=timepoint, label=timepoint),
                    alpha=0.1,
                    linetype = 2) +
  geom_point(aes(shape=treatment),size=3, alpha=1) +
  scale_shape_manual(name="Treatment",
                     labels=c('Control', 
                              expression(paste(italic("F. nucleatum"))), 
                              expression(paste(italic("F. nucleatum/S. salivarius"))), 
                              expression(paste(italic('S. salivarius')))),
                     values = c(15,16,17,18)) +
  ggtitle("nMDS of weighted UNIFRAC distances")
p1 <- p1 + expand_limits(x = c(-.21,.25), y = c(-.16, .16))
p1
#Bray Curtis

#distance matrix, BC, classic (metric) multidimensional scaling (PCA)
dist_mat = vegdist(mpa_mat, method="bray")
cmd_res = cmdscale(dist_mat, 
                   k = (nrow(mpa_mat) - 1),
                   eig = TRUE)

# calculate the proportion of variance in the data which is explained by 
#the first two PCoA axes
PC1 <- round(100*(cmd_res$eig[1]/(sum(cmd_res$eig))), digits = 2)
PC2 <- round(100*(cmd_res$eig[2]/(sum(cmd_res$eig))), digits = 2)
#extract details of features
pcoa_df = tibble(PC1 = cmd_res$points[,1], 
                 PC2 = cmd_res$points[,2])

#combine pcoa with sample data
meta$timepoint <- factor(meta$timepoint,levels = c("T_6", "T_24"))
meta$treatment <- factor(meta$treatment)
pcoa_meta = bind_cols(pcoa_df, meta)
p2 = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=timepoint)) + 
  scale_color_manual(name="Timepoint",
                     labels=c("T_6",
                              "T_24"),
                     values = wes_palette("Darjeeling1", n = 3, type = "discrete")) +
  scale_fill_manual(name="Timepoint",
                    labels=c("T_6",
                             "T_24"),
                    values = wes_palette("Darjeeling1", n = 3, type = "discrete")) +
  labs(shape="Treatment") +
  theme(axis.title.x= element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=16),
        legend.position="bottom", legend.box = "horizontal",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid")) +
  guides(col = "none",
         shape = "none",
         fill = "none") +
  ylab(paste("PC2 (",PC2,"%)"))+
  xlab(paste("PC1 (",PC1,"%)"))+
  geom_mark_ellipse(aes(fill=timepoint, label=timepoint),
                    alpha=0.1,
                    linetype = 2) +
  geom_point(aes(shape=treatment),size=3, alpha=1) +
  scale_shape_manual(name="Treatment",
                     labels=c('Control', 
                              expression(paste(italic("F. nucleatum"))), 
                              expression(paste(italic("F. nucleatum/S. salivarius"))), 
                              expression(paste(italic('S. salivarius')))),
                     values = c(15,16,17,18)) +
  ggtitle("PCoA of Bray-Curtis distances")
p2 <- p2 + expand_limits(x=c(-.2,.2),y = c(-.2, .2))
plots <- plot_grid(p0, p1, p2, plegend)
plots
ggsave("figures/beta_T6_T24.pdf", width=29, height=21, dpi=300, units="cm")

#Figure X Taxonomic β-diversity of sampled timepoints