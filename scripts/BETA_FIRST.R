#' Script to calculate beta diversity on NisinG samples
#' Developed by A. Kate Falà (https://github.com/aforestsomewhere)
#' all groups 

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
#rm()
#dev.off()
#read sample metadata
meta <- read.csv("cleaned_data/metadata.csv", nrows=36, row.names=1)
#read MPA4 data
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

####################
### normalisation ##
####################
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9499858/
#arcsin-sqrt
#stabilises variance where values are between 0 and 1 - need to convert first
#create new vector where each value is divided by max value
pre_asin_sqrt <- function(x, max_x){
  y <- x / max_x
  return(y)
}
asin_sqrt <- function(x){
  y <- asin(sqrt(x))
  return(y)
}
mpa_mat_arc <- pre_asin_sqrt(mpa_mat, max(mpa_mat))
mpa_mat_arc <- asin_sqrt(mpa_mat_arc)

################
### ANOSIM   ###
################
#ANOSIM (ANalysis Of Similarities) is a non-parametric test of significant difference 
#between two or more groups, based on any distance measure (Clarke 1993). The distances 
#are converted to ranks. ANOSIM is normally used for taxa-in-samples data, where groups 
#of samples are to be compared.

dist_mat = vegdist(mpa_mat_arc, method="bray")
anosim_df_all <- data.frame(Data="MetaPhlAn4", 
                            Metric="Bray-Curtis dissimilarity",
                            Comparing=NA,
                            Test="ANOSIM",
                            Test_Statistic=NA,
                            P_value=NA,
                            Significance_codes="NA") #initialise output DF
#timepoint
anosim_time <- anosim(dist_mat,meta$timepoint, permutations = 9999) #TIMEPOINT
P<- anosim_time$signif
F<-anosim_time$statistic
S<-sigcode(P)
anosim_df <- c("MetaPhlAn4", "Bray-Curtis dissimilarity", "Timepoint", "ANOSIM",F,P,S)
anosim_df_all <-rbind(anosim_df_all, anosim_df)
anosim_df_all[-1,]-> anosim_df_all #remove first shaper row (empty)
#treatment
anosim_treat <- anosim(dist_mat,meta$treatment, permutations = 9999) #TREATMENT
P<- anosim_treat$signif
F<-anosim_treat$statistic
S<-sigcode(P)
name<-"anosim_treatment_BC"
anosim_df <- c("MetaPhlAn4", "Bray-Curtis dissimilarity", "Treatment", "ANOSIM",F,P,S)
anosim_df_all <-rbind(anosim_df_all, anosim_df)
#exp_group
anosim_treat <- anosim(dist_mat,meta$exp_group, permutations = 9999) #TREATMENT
P<- anosim_treat$signif
F<-anosim_treat$statistic
S<-sigcode(P)
name<-"anosim_exp_group_BC"
anosim_df <- c("MetaPhlAn4", "Bray-Curtis dissimilarity", "Exp_group", "ANOSIM",F,P,S)
anosim_df_all <-rbind(anosim_df_all, anosim_df)

#plot
############################################
### Fig 1 BC arcsin MPA4 sp level         ###
#############################################
#distance matrix, BC, classic (metric) multidimensional scaling (PCA)
dist_mat = vegdist(mpa_mat_arc, method="bray")
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
meta$timepoint <- factor(meta$timepoint,levels = c("T_0","T_6", "T_24"))
meta$treatment <- factor(meta$treatment)
pcoa_meta = bind_cols(pcoa_df, meta)

#base plot to extract legend
leggo = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=treatment, fill=treatment)) + 
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
  scale_shape(name="Timepoint")+
  labs(shape="Timepoint") +
  geom_point(aes(shape = timepoint), size=4) +
  ylab(paste("PC2 (",PC2,"%)"))+
  xlab(paste("PC1 (",PC1,"%)"))+
  theme(axis.title.x= element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=18, face = "bold"),
        legend.text = element_text(size=16),
        # legend.position="bottom", legend.box = "horizontal",
        legend.position = c(0.3, 0.85), legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid")) +
  xlim(c(-0.3,0.3))+
  ylim(c(-0.2,0.2))+
  geom_mark_ellipse(aes(fill=treatment, color=treatment),
                    alpha=0.1,
                    linetype = 2) +
  guides(col = guide_legend(order = 1),
         shape = guide_legend(order = 0),
         fill = "none") +
  ggtitle("PCoA, BC, arcsin-transformed MPA4, sp. level assignment")
plegend <- get_legend(leggo)
#####################################################################################
p2 = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=treatment, fill=treatment)) + 
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
  scale_shape(name="Timepoint")+
  labs(shape="Timepoint") +
  geom_point(aes(shape = timepoint), size=4) +
  ylab(paste("PC2 (",PC2,"%)"))+
  xlab(paste("PC1 (",PC1,"%)"))+
  theme(axis.title.x= element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=18, face = "bold"),
        legend.text = element_text(size=16),
        # legend.position="bottom", legend.box = "horizontal",
        legend.position = c(0.3, 0.85), legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid")) +
  xlim(c(-0.3,0.3))+
  ylim(c(-0.2,0.2))+
  geom_mark_ellipse(aes(fill=treatment, color=treatment),
                    alpha=0.1,
                    linetype = 2) +
  guides(col = "none",
         shape = "none",
         fill = "none") +
  ggtitle("PCoA, BC, arcsin-transformed MPA4, sp. level assignment")
p2 


#aitchison
#Now the same using Aitchison distance. This metric corresponds to Euclidean 
#distances between CLR transformed sample abundance vectors.
# Does clr transformation. Pseudocount is added, because data contains zeros.
tse <- mia::transformCounts(mpa_mat,  assay_name = abund_values, method = "clr", pseudocount = 1)

dist_mat <-  vegdist(mpa_mat, method="aitchison",pseudocount = 1)
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
meta$timepoint <- factor(meta$timepoint,levels = c("T_0","T_6", "T_24"))

pcoa_meta = bind_cols(pcoa_df, meta)
p3 = ggplot(pcoa_meta,
            aes(x = PC1, y = PC2, color = treatment, shape=timepoint)) + 
  scale_color_manual(values = wes_palette("Darjeeling1", n = 4, type = "discrete"))+
  geom_point(size=3)+
  ggtitle ("Logit-transformed, Euclidean distances")+
  ylab(paste("PC2 (",PC2,"%)"))+
  xlab(paste("PC1 (",PC1,"%)"))
p3

plots <- plot_grid(p0, p1, p2)
plots
