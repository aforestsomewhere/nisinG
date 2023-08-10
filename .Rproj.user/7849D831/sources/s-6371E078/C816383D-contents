#final beta

#' Script to assess beta diversity for NisinG samples
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
library(ggforce)
theme_set(theme_classic())
palette_time <- c("#9A8822", "#F8AFA8", "#74A089")
palette_groups <- c("#332211","#ffcc44", "#44aacc", "#bb2211")
palette_collapsed <- c("#FF0000", "#00A08A")

####################
### functions    ###
####################

#new vectorised
sigcode <- function(P) {
  sapply(P, function(p) {
    if (p <= 0) {
      return("NaN")
    } else if (p < 0.001) {
      return("***")
    } else if (p < 0.01) {
      return("**")
    } else if (p <= 0.05) {
      return("*")
    } else if (p < 0.1) {
      return(".")
    } else {
      return(" ")
    }
  })
}
#Functions to adjust X and Y axes for PCoA computed through Phyloseq
new_ylab <- function(plot) {
  y_axis <- get_plot_component(plot, "ylab-l")
  ylab <- y_axis$children[[1]]$label
  # Find the position of the opening and closing square brackets
  start_pos <- regexpr("\\[", ylab)
  end_pos <- regexpr("\\]", ylab)
  # Extract the substring between the square brackets
  percent <- substring(ylab, start_pos+1, end_pos-1)
  result <- paste("PC2 (",percent,")")
  return(result)
}
new_xlab <- function(plot) {
  x_axis <- get_plot_component(plot, "xlab-b")
  xlab <- x_axis$children[[1]]$label
  # Find the position of the opening and closing square brackets
  start_pos <- regexpr("\\[", xlab)
  end_pos <- regexpr("\\]", xlab)
  # Extract the substring between the square brackets
  percent <- substring(xlab, start_pos + 1, end_pos - 1)
  result <- paste("PC1 (",percent,")")
  return(result)
}

####################
### read data    ###
####################
meta <- read.csv("cleaned_data/metadata.csv", nrows=36, row.names=1)
meta <- meta %>%
  mutate(strep_exposure = case_when(
    treatment == "FS" ~ "Y", 
    treatment == "S" ~ "Y",
    treatment == "F" ~ "N",
    treatment == "B" ~ "N"))


mpa_table <- read.table("cleaned_data/mpa_table.tsv", comment.char = '#', sep = '\t', header = TRUE, row.names=NULL)
#filter to species level assignments
mpa_table <- mpa_table[grep('s__',mpa_table[,1]),]
mpa_table[,1] <- gsub(".+\\|s__", "", mpa_table[,1])
#remove duplicates by filtering to t__ level
#mpa_table <- mpa_table %>% filter(!grepl('t__',row.names))
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]
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
#filter out dupes

keep_rows <- grepl("t__", row.names(mpa_table), perl = TRUE)
mpa_table_filtered <- mpa_table[!keep_rows, ]
#generate matrix of abundance data
mpa_mat = t(as.matrix(mpa_table_filtered))
#order samples (rows) to properly merge with infO
mpa_mat <- mpa_mat[sort(rownames(mpa_mat)),]
dist_mat_ai = vegdist(mpa_mat, method="robust.aitchison")
dist_mat_bc = vegdist(mpa_mat, methd = "bray")
#UNIFRAC
#Weighted, unnormalised
phy@sam_data$timepoint = factor(phy@sam_data$timepoint, levels = c("T_0", "T_6", "T_24"))
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
                            Metric="Robust Aitchison Distances",
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
#Bray
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_ai,group=meta$timepoint, dataset="MetaPhlAn4", metric="Bray-Curtis Distances", groupname="Timepoint"))
anosim_df_all[-1,]-> anosim_df_all #remove first shaper row (empty)
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_ai,group=meta$treatment, dataset="MetaPhlAn4", metric="Bray-Curtis Distances", groupname="Treatment"))
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_ai,group=meta$strep_exposure, dataset="MetaPhlAn4", metric="Bray-Curtis Distances", groupname="Any exposure"))
#Robust Aitchisons
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_ai,group=meta$timepoint, dataset="MetaPhlAn4", metric="Robust Aitchison Distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_ai,group=meta$treatment, dataset="MetaPhlAn4", metric="Robust Aitchison Distances", groupname="Treatment"))
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_ai,group=meta$strep_exposure, dataset="MetaPhlAn4", metric="Robust Aitchison Distances", groupname="Any exposure"))
#Weighted UNIFRAC
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_wUF,group=meta$timepoint, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_wUF,group=meta$treatment, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Treatment"))
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_wUF,group=meta$strep_exposure, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Any exposure"))
#Unweighted UNIFRAC
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_uUF,group=meta$timepoint, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_uUF,group=meta$treatment, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Treatment"))
anosim_df_all <-rbind(anosim_df_all, calc_anosim(dist_mat_uUF,group=meta$strep_exposure, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Any exposure"))

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
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_bc,group=meta$strep_exposure, dataset="MetaPhlAn4", metric="Bray-Curtis dissimilarity", groupname="Any exposure"))
#Robust Aitchisons
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_ai,group=meta$timepoint, dataset="MetaPhlAn4", metric="Robust Aitchison distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_ai,group=meta$treatment, dataset="MetaPhlAn4", metric="Robust Aitchison distances", groupname="Treatment"))
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_ai,group=meta$strep_exposure, dataset="MetaPhlAn4", metric="Robust Aitchison distances", groupname="Any exposure"))
#Weighted UNIFRAC
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_wUF,group=meta$timepoint, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_wUF,group=meta$treatment, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Treatment"))
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_wUF,group=meta$strep_exposure, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Any exposure"))
#Unweighted UNIFRAC
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_uUF,group=meta$timepoint, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_uUF,group=meta$treatment, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Treatment"))
anosim_df_all <-rbind(anosim_df_all, calc_adonis(dist_mat_uUF,group=meta$strep_exposure, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Any exposure"))

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
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_bc,group=meta$strep_exposure, dataset="MetaPhlAn4", metric="Bray-Curtis dissimilarity", groupname="Any exposure"))
#Robust Aitchisons
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_ai,group=meta$timepoint, dataset="MetaPhlAn4", metric="Robust Aitchison distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_ai,group=meta$treatment, dataset="MetaPhlAn4", metric="Robust Aitchison distances", groupname="Treatment"))
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_ai,group=meta$strep_exposure, dataset="MetaPhlAn4", metric="Robust Aitchison distances", groupname="Any exposure"))
#Weighted UNIFRAC
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_wUF,group=meta$timepoint, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_wUF,group=meta$treatment, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Treatment"))
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_wUF,group=meta$strep_exposure, dataset="MetaPhlAn4", metric="Weighted UNIFRAC distances", groupname="Any exposure"))
#Unweighted UNIFRAC
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_uUF,group=meta$timepoint, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Timepoint"))
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_uUF,group=meta$treatment, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Treatment"))
anosim_df_all <-rbind(anosim_df_all, calc_disp(dist_mat_uUF,group=meta$strep_exposure, dataset="MetaPhlAn4", metric="Unweighted UNIFRAC distances", groupname="Any exposure"))

beta_stats_all <- anosim_df_all
write.csv(beta_stats_all, "reports/beta_alltimepoints_testing.csv")

###############################
# Plots, all timepoints       #
###############################
###############################
# Robust Aitchisons           #
###############################

cmd_res = cmdscale(dist_mat_ai, 
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
pcoa_meta = bind_cols(pcoa_df, meta)
#timepoint

p2 = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=timepoint)) + 
  scale_color_manual(name="Timepoint",
                     labels=c("T_0", 
                              "T_6",
                              "T_24"),
                     values = palette_time) +
  scale_fill_manual(name="Timepoint",
                    labels=c("T_0", 
                             "T_6",
                             "T_24"),
                    values = palette_time) +
  labs(shape="Treatment") +
  theme(axis.title.x= element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.key.size = unit(0.5, 'mm'),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size=5),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.position="bottom", legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  ylab(paste("PC2 (",PC2,"%)", sep=""))+
  xlab(paste("PC1 (",PC1,"%)", sep=""))+
  geom_mark_ellipse(aes(fill=timepoint),
                    alpha=0.05,
                    linetype = 2) +
  geom_point(aes(shape=treatment),size=1, alpha=1) +
  scale_shape_manual(name="Treatment",
                     labels=c("Control", #b
                              expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                              expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                              expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                     values = c(15,16,17,18)) +
  guides(fill= "none",
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=palette_time)),
         shape = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE)) +
  ggtitle("Robust Aitchisons: All timepoints ~ timepoint")
p2a <- p2 + expand_limits(x = c(-16,11), y = c(-16,8))
p2leg <- get_legend(p2a)
p2b <- p2a + theme(legend.position='none')
p2c <- ggdraw(plot_grid(p2b, p2leg,
                        nrow=2, rel_heights=c(10, 2),
                        rel_widths=c(17, 2),
                        axis = 'l', align = 'h'))
p2c

#ggsave("figures/beta_aitchisonsp2.pdf", width=29, height=21, dpi=300, units="cm")

#group
p3 = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=treatment)) + 
  scale_color_manual(name="Treatment",
                     labels=c("Control", #b
                              expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                              expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                              expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                     values = palette_groups) + 
  scale_fill_manual(name="Treatment",
                    labels=c("Control", #b
                             expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                             expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                             expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                    values = palette_groups) +
  labs(shape="Treatment") +
  theme(axis.title.x= element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.key.size = unit(0.5, 'mm'),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size=5),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.position="bottom", legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  ylab(paste("PC2 (",PC2,"%)",sep=""))+
  xlab(paste("PC1 (",PC1,"%)",sep=""))+
  geom_mark_ellipse(aes(fill=treatment),
                    alpha=0.05,
                    linetype = 2) +
  geom_point(aes(shape=timepoint),size=1, alpha=1) +
  scale_shape_manual(name="Timepoint",
                     labels=c("T_0", 
                              "T_6",
                              "T_24"),
                     values =  c(15,16,17,18)) +
  guides(fill= "none",
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=palette_groups)),
         shape = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE)) +
  ggtitle("Robust Aitchisons: All timepoints ~ treatment")
p3
p3a <- p3 + expand_limits(x = c(-17,9), y = c(-15,10))
p3a
p3leg <- get_legend(p3a)
p3b <- p3a + theme(legend.position='none')
p3c <- ggdraw(plot_grid(p3b, p3leg,
                        nrow=2, rel_heights=c(10, 2),
                        rel_widths=c(17, 2),
                        axis = 'l', align = 'h'))

##########################
#any strep exposure
p4 = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=strep_exposure)) + 
  scale_color_manual(name="Any exposure",
                     labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                              expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                     values = palette_collapsed) + 
  scale_fill_manual(name="Any exposure",
                    labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                             expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                    values = palette_collapsed) +
  labs(shape="Timepoint") +
  theme(axis.title.x= element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.key.size = unit(0.5, 'mm'),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size=5),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.position="bottom", legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  ylab(paste("PC2 (",PC2,"%)",sep=""))+
  xlab(paste("PC1 (",PC1,"%)",sep=""))+
  geom_mark_ellipse(aes(fill=strep_exposure),
                    alpha=0.05,
                    linetype = 2) +
  geom_point(aes(shape=timepoint),size=1, alpha=1) +
  scale_shape_manual(name="Timepoint",
                     labels=c("T_0", 
                              "T_6",
                              "T_24"),
                     values =  c(15,16,17,18)) +
  guides(fill= "none",
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=palette_collapsed)),
         shape = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE)) +
  ggtitle("Robust Aitchisons: All timepoints ~ exposure")
p4
p4a <- p4 + expand_limits(x = c(-15,11), y = c(-12,8))
p4a
p4leg <- get_legend(p4a)
p4b <- p4a + theme(legend.position='none')
p4c <- ggdraw(plot_grid(p4b, p4leg,
                        nrow=2, rel_heights=c(10, 2),
                        rel_widths=c(17, 2),
                        axis = 'l', align = 'h'))

p0 <- ggdraw(plot_grid(p2c,p3c,p4c,
                       align = "h", ncol = 1,
                       labels = "AUTO"))
p0
ggsave("figures/beta_all_time_aitchisons.png", width=13, height=29, dpi=300, units="cm")
ggsave("figures/beta_all_time_aitchisons.pdf", width=13, height=29, dpi=300, units="cm")

###############################
# Bray Curtis                 #
###############################

cmd_res = cmdscale(dist_mat_bc, 
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
pcoa_meta = bind_cols(pcoa_df, meta)
#timepoint

p2 = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=timepoint)) + 
  scale_color_manual(name="Timepoint",
                     labels=c("T_0", 
                              "T_6",
                              "T_24"),
                     values = palette_time) +
  scale_fill_manual(name="Timepoint",
                    labels=c("T_0", 
                             "T_6",
                             "T_24"),
                    values = palette_time) +
  labs(shape="Treatment") +
  theme(axis.title.x= element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.key.size = unit(0.5, 'mm'),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size=5),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.position="bottom", legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  ylab(paste("PC2 (",PC2,"%)", sep=""))+
  xlab(paste("PC1 (",PC1,"%)", sep=""))+
  geom_mark_ellipse(aes(fill=timepoint),
                    alpha=0.05,
                    linetype = 2) +
  geom_point(aes(shape=treatment),size=1, alpha=1) +
  scale_shape_manual(name="Treatment",
                     labels=c("Control", #b
                              expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                              expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                              expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                     values = c(15,16,17,18)) +
  guides(fill= "none",
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=palette_time)),
         shape = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE)) +
  ggtitle("Bray-Curtis: All timepoints ~ timepoint")
p2
p2a <- p2 + expand_limits(x = c(-.3,.3), y = c(-.2,.2))
p2leg <- get_legend(p2a)
p2b <- p2a + theme(legend.position='none')
p2c <- ggdraw(plot_grid(p2b, p2leg,
                        nrow=2, rel_heights=c(10, 2),
                        rel_widths=c(17, 2),
                        axis = 'l', align = 'h'))
p2c

#ggsave("figures/beta_aitchisonsp2.pdf", width=29, height=21, dpi=300, units="cm")

#group
p3 = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=treatment)) + 
  scale_color_manual(name="Treatment",
                     labels=c("Control", #b
                              expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                              expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                              expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                     values = palette_groups) + 
  scale_fill_manual(name="Treatment",
                    labels=c("Control", #b
                             expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                             expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                             expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                    values = palette_groups) +
  labs(shape="Treatment") +
  theme(axis.title.x= element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.key.size = unit(0.5, 'mm'),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size=5),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.position="bottom", legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  ylab(paste("PC2 (",PC2,"%)",sep=""))+
  xlab(paste("PC1 (",PC1,"%)",sep=""))+
  geom_mark_ellipse(aes(fill=treatment),
                    alpha=0.05,
                    linetype = 2) +
  geom_point(aes(shape=timepoint),size=1, alpha=1) +
  scale_shape_manual(name="Timepoint",
                     labels=c("T_0", 
                              "T_6",
                              "T_24"),
                     values =  c(15,16,17,18)) +
  guides(fill= "none",
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=palette_groups)),
         shape = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE)) +
  ggtitle("Bray-Curtis: All timepoints ~ treatment")
p3
p3a <- p3 + expand_limits(x = c(-.3,.3), y = c(-.2,.2))
p3a
p3leg <- get_legend(p3a)
p3b <- p3a + theme(legend.position='none')
p3c <- ggdraw(plot_grid(p3b, p3leg,
                        nrow=2, rel_heights=c(10, 2),
                        rel_widths=c(17, 2),
                        axis = 'l', align = 'h'))

##########################
#any strep exposure
p4 = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=strep_exposure)) + 
  scale_color_manual(name="Any exposure",
                     labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                              expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                     values = palette_collapsed) + 
  scale_fill_manual(name="Any exposure",
                    labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                             expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                    values = palette_collapsed) +
  labs(shape="Timepoint") +
  theme(axis.title.x= element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.key.size = unit(0.5, 'mm'),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size=5),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.position="bottom", legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  ylab(paste("PC2 (",PC2,"%)",sep=""))+
  xlab(paste("PC1 (",PC1,"%)",sep=""))+
  geom_mark_ellipse(aes(fill=strep_exposure),
                    alpha=0.05,
                    linetype = 2) +
  geom_point(aes(shape=timepoint),size=1, alpha=1) +
  scale_shape_manual(name="Timepoint",
                     labels=c("T_0", 
                              "T_6",
                              "T_24"),
                     values =  c(15,16,17,18)) +
  guides(fill= "none",
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=palette_collapsed)),
         shape = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE)) +
  ggtitle("Bray-Curtis: All timepoints ~ exposure")
p4
p4a <- p4 + expand_limits(x = c(-.3,.3), y = c(-.2,.2))
p4a
p4leg <- get_legend(p4a)
p4b <- p4a + theme(legend.position='none')
p4c <- ggdraw(plot_grid(p4b, p4leg,
                        nrow=2, rel_heights=c(10, 2),
                        rel_widths=c(17, 2),
                        axis = 'l', align = 'h'))

p0 <- ggdraw(plot_grid(p2c,p3c,p4c,
                       align = "h", ncol = 1,
                       labels = "AUTO"))
p0
ggsave("figures/beta_all_time_bc.png", width=13, height=29, dpi=300, units="cm")
ggsave("figures/beta_all_time_bc.pdf", width=13, height=29, dpi=300, units="cm")

###############################
# Unweighted UNIFRAC            #
###############################

uUF.ordu = ordinate(phy, method="NMDS", distance="unifrac", weighted=FALSE)
#timepoint
p2 <- plot_ordination(phy, uUF.ordu, type="sites", color="timepoint") +
  scale_color_manual(name="Timepoint",
                     labels=c("T_0", 
                              "T_6",
                              "T_24"),
                     values = palette_time) +
  scale_fill_manual(name="Timepoint",
                    labels=c("T_0", 
                             "T_6",
                             "T_24"),
                    values = palette_time) +
  labs(shape="Treatment") +
  theme(axis.title.x= element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.key.size = unit(0.5, 'mm'),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size=5),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.position="bottom", legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  #ylab(paste("PC2 (",PC2,"%)", sep=""))+
 # xlab(paste("PC1 (",PC1,"%)", sep=""))+
  geom_mark_ellipse(aes(fill=timepoint),
                    alpha=0.05,
                    linetype = 2) +
  geom_point(aes(shape=treatment),size=1, alpha=1) +
  scale_shape_manual(name="Treatment",
                     labels=c("Control", #b
                              expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                              expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                              expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                     values = c(15,16,17,18)) +
  guides(fill= "none",
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=palette_time)),
         shape = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE)) +
  annotate("text", x=-Inf,y=-Inf,hjust=-0.1,vjust=-.98,label=(paste("Stress=",round(uUF.ordu$stress, digits = 6))))+
  ggtitle("Unweighted UNIFRAC: All timepoints ~ timepoint")

p2c <- p2 + expand_limits(x = c(-.4,.4), y = c(-.18, .13))
#group
p3 = plot_ordination(phy, uUF.ordu, type="sites", color="treatment") +
  scale_color_manual(name="Treatment",
                     labels=c("Control", #b
                              expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                              expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                              expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                     values = palette_groups) + 
  scale_fill_manual(name="Treatment",
                    labels=c("Control", #b
                             expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                             expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                             expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                    values = palette_groups) +
  labs(shape="Treatment") +
  theme(axis.title.x= element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.key.size = unit(0.5, 'mm'),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size=5),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.position="bottom", legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  #ylab(paste("PC2 (",PC2,"%)",sep=""))+
  #xlab(paste("PC1 (",PC1,"%)",sep=""))+
  geom_mark_ellipse(aes(fill=treatment),
                    alpha=0.05,
                    linetype = 2) +
  geom_point(aes(shape=timepoint),size=1, alpha=1) +
  scale_shape_manual(name="Timepoint",
                     labels=c("T_0", 
                              "T_6",
                              "T_24"),
                     values =  c(15,16,17,18)) +
  guides(fill= "none",
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=palette_groups)),
         shape = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE)) +
  annotate("text", x=-Inf,y=-Inf,hjust=-0.1,vjust=-.98,label=(paste("Stress=",round(uUF.ordu$stress, digits = 6))))+
  ggtitle("Unweighted UNIFRAC: All timepoints ~ treatment")
p3a <- p3 + expand_limits(x = c(-.3,.3), y = c(-.18,.14))
p3a
p3leg <- get_legend(p3a)
p3b <- p3a + theme(legend.position='none')
p3c <- ggdraw(plot_grid(p3b, p3leg,
                        nrow=2, rel_heights=c(10, 2),
                        rel_widths=c(17, 2),
                        axis = 'l', align = 'h'))
##########################
#any strep exposure
p4 = plot_ordination(phy, uUF.ordu, type="sites", color="strep_exposure") +
  scale_color_manual(name="Any exposure",
                     labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                              expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                     values = palette_collapsed) + 
  scale_fill_manual(name="Any exposure",
                    labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                             expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                    values = palette_collapsed) +
  labs(shape="Timepoint") +
  theme(axis.title.x= element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.key.size = unit(0.5, 'mm'),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size=5),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.position="bottom", legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
 # ylab(paste("PC2 (",PC2,"%)",sep=""))+
  #xlab(paste("PC1 (",PC1,"%)",sep=""))+
  geom_mark_ellipse(aes(fill=strep_exposure),
                    alpha=0.05,
                    linetype = 2) +
  geom_point(aes(shape=timepoint),size=1, alpha=1) +
  scale_shape_manual(name="Timepoint",
                     labels=c("T_0", 
                              "T_6",
                              "T_24"),
                     values =  c(15,16,17,18)) +
  guides(fill= "none",
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=palette_collapsed)),
         shape = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE)) +
  annotate("text", x=-Inf,y=-Inf,hjust=-0.1,vjust=-.98,label=(paste("Stress=",round(uUF.ordu$stress, digits = 6))))+
  ggtitle("Unweighted UNIFRAC: All timepoints ~ exposure")
p4
p4a <- p4 + expand_limits(x = c(-.3,.3), y = c(-.2,.2))
p4a
p4leg <- get_legend(p4a)
p4b <- p4a + theme(legend.position='none')
p4c <- ggdraw(plot_grid(p4b, p4leg,
                        nrow=2, rel_heights=c(10, 2),
                        rel_widths=c(17, 2),
                        axis = 'l', align = 'h'))

p0 <- ggdraw(plot_grid(p2c,p3c,p4c,
                       align = "h", ncol = 1,
                       labels = "AUTO"))
p0
ggsave("figures/beta_all_time_uUF.png", width=13, height=29, dpi=300, units="cm")
ggsave("figures/beta_all_time_uUF.pdf", width=13, height=29, dpi=300, units="cm")

###############################
# Unweighted UNIFRAC            #
###############################

wUF.ordu = ordinate(phy, method="NMDS", distance="unifrac", weighted=TRUE)
#timepoint
p2 <- plot_ordination(phy, wUF.ordu, type="sites", color="timepoint") +
  scale_color_manual(name="Timepoint",
                     labels=c("T_0", 
                              "T_6",
                              "T_24"),
                     values = palette_time) +
  scale_fill_manual(name="Timepoint",
                    labels=c("T_0", 
                             "T_6",
                             "T_24"),
                    values = palette_time) +
  labs(shape="Treatment") +
  theme(axis.title.x= element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.key.size = unit(0.5, 'mm'),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size=5),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.position="bottom", legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  #ylab(paste("PC2 (",PC2,"%)", sep=""))+
  #xlab(paste("PC1 (",PC1,"%)", sep=""))+
  geom_mark_ellipse(aes(fill=timepoint),
                    alpha=0.05,
                    linetype = 2) +
  geom_point(aes(shape=treatment),size=1, alpha=1) +
  scale_shape_manual(name="Treatment",
                     labels=c("Control", #b
                              expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                              expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                              expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                     values = c(15,16,17,18)) +
  guides(fill= "none",
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=palette_time)),
         shape = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE)) +
  annotate("text", x=-Inf,y=-Inf,hjust=-0.1,vjust=-.98,label=(paste("Stress=",round(wUF.ordu$stress, digits = 6))))+
  ggtitle("Weighted UNIFRAC: All timepoints ~ timepoint")

p2c <- p2 + expand_limits(x = c(-.3,.23), y = c(-.18, .13))
p2c
#group
p3 = plot_ordination(phy, wUF.ordu, type="sites", color="treatment") +
  scale_color_manual(name="Treatment",
                     labels=c("Control", #b
                              expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                              expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                              expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                     values = palette_groups) + 
  scale_fill_manual(name="Treatment",
                    labels=c("Control", #b
                             expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                             expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                             expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                    values = palette_groups) +
  labs(shape="Treatment") +
  theme(axis.title.x= element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.key.size = unit(0.5, 'mm'),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size=5),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.position="bottom", legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  #ylab(paste("PC2 (",PC2,"%)",sep=""))+
  #xlab(paste("PC1 (",PC1,"%)",sep=""))+
  geom_mark_ellipse(aes(fill=treatment),
                    alpha=0.05,
                    linetype = 2) +
  geom_point(aes(shape=timepoint),size=1, alpha=1) +
  scale_shape_manual(name="Timepoint",
                     labels=c("T_0", 
                              "T_6",
                              "T_24"),
                     values =  c(15,16,17,18)) +
  guides(fill= "none",
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=palette_groups)),
         shape = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE)) +
  annotate("text", x=-Inf,y=-Inf,hjust=-0.1,vjust=-.98,label=(paste("Stress=",round(wUF.ordu$stress, digits = 6))))+
  ggtitle("Weighted UNIFRAC: All timepoints ~ treatment")
p3
p3a <- p3 + expand_limits(x = c(-.25,.25), y = c(-.14,.10))
p3a
p3leg <- get_legend(p3a)
p3b <- p3a + theme(legend.position='none')
p3c <- ggdraw(plot_grid(p3b, p3leg,
                        nrow=2, rel_heights=c(10, 2),
                        rel_widths=c(17, 2),
                        axis = 'l', align = 'h'))
##########################
#any strep exposure
p4 = plot_ordination(phy, wUF.ordu, type="sites", color="strep_exposure") +
  scale_color_manual(name="Any exposure",
                     labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                              expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                     values = palette_collapsed) + 
  scale_fill_manual(name="Any exposure",
                    labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                             expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                    values = palette_collapsed) +
  labs(shape="Timepoint") +
  theme(axis.title.x= element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.key.size = unit(0.5, 'mm'),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size=5),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.position="bottom", legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  #ylab(paste("PC2 (",PC2,"%)",sep=""))+
  #xlab(paste("PC1 (",PC1,"%)",sep=""))+
  geom_mark_ellipse(aes(fill=strep_exposure),
                    alpha=0.05,
                    linetype = 2) +
  geom_point(aes(shape=timepoint),size=1, alpha=1) +
  scale_shape_manual(name="Timepoint",
                     labels=c("T_0", 
                              "T_6",
                              "T_24"),
                     values =  c(15,16,17,18)) +
  guides(fill= "none",
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=palette_collapsed)),
         shape = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE)) +
  annotate("text", x=-Inf,y=-Inf,hjust=-0.1,vjust=-.98,label=(paste("Stress=",round(wUF.ordu$stress, digits = 6))))+
  ggtitle("Weighted UNIFRAC: All timepoints ~ exposure")
p4
p4a <- p4 + expand_limits(x = c(-.25,.25), y = c(-.14,.10))
p4a
p4leg <- get_legend(p4a)
p4b <- p4a + theme(legend.position='none')
p4c <- ggdraw(plot_grid(p4b, p4leg,
                        nrow=2, rel_heights=c(10, 2),
                        rel_widths=c(17, 2),
                        axis = 'l', align = 'h'))

p0 <- ggdraw(plot_grid(p2c,p3c,p4c,
                       align = "h", ncol = 1,
                       labels = "AUTO"))
p0
ggsave("figures/beta_all_time_wUF.png", width=13, height=29, dpi=300, units="cm")
ggsave("figures/beta_all_time_wUF.pdf", width=13, height=29, dpi=300, units="cm")

