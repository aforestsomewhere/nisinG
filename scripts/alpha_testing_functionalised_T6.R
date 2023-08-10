#' Script to calculate alpha diversity on NisinG samples
#' Developed by A. Kate Falà (https://github.com/aforestsomewhere)
#' collapsing any S. salivarius-treated; 

##############################
### Packages and palettes  ###
##############################
library(tidyverse)
library(dplyr)
library(vegan)
library(cowplot)
theme_set(theme_cowplot(font_size=20))
library(devEMF)
library(ggforce)
library(webshot)
library(ggplot2)
library(rstatix)
library(ggsignif)
library(showtext)
library(ggpubr)
pal1 <- c("#072e51","#ffdf66","#332211", "#ffcc44", "#44aacc", "#bb2211")
#install.packages("showtext")
#library('devtools')
#theme_name = "theme_starwars" # Pick which theme you want
#theme_url = paste0("https://raw.githubusercontent.com/MatthewBJane/theme_park/main/", theme_name ,".R")
#source_url(theme_url)

################
# FUNCTIONS    #
################


####################
### read data    ###
####################
#rm()
#dev.off()
#read sample metadata
meta <- read.csv("cleaned_data/metadata.csv", nrows=36, row.names=1)
meta <- meta %>% mutate(treatment = str_replace(treatment, "FS", "S")) %>%
  mutate(treatment = str_replace(treatment, "F", "B")) %>%
  filter(timepoint=="T_6")
#read MPA4 data
mpa_table <- read.table("cleaned_data/mpa_table.tsv", comment.char = '#', sep = '\t', header = TRUE, row.names=NULL)
#filter to species level assignments
mpa_table <- mpa_table[grep('s__',mpa_table[,1]),]
mpa_table[,1] <- gsub(".+\\|s__", "", mpa_table[,1])
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]
#filter abundance data to timepoint T_6 only
samples <- rownames(meta)
mpa_table[,colnames(mpa_table) %in% samples] -> mpa_table
#generate matrix of abundance data
mpa_mat = t(as.matrix(mpa_table))
#order samples (rows) to properly merge with infO
mpa_mat <- mpa_mat[sort(rownames(mpa_mat)),] 

####################################
### calculate alpha diversity    ###
####################################
calc_vegan_table <- function(mpa, metadata, time){
  richness_vec <-(specnumber(mpa))
  shannon_vec <- vegan::diversity(mpa, index = "shannon")
  simp_vec <- vegan::diversity(mpa, index = "simpson")
  inv_simp_vec <- vegan::diversity(mpa, index = "invsimpson")
  alpha_values <- cbind.data.frame(comb_id = rownames(meta), 
                                   richness = richness_vec, 
                                   shannon = shannon_vec, 
                                   simpson = simp_vec, 
                                   inv_simpson = inv_simp_vec, 
                                   meta = meta)
  names(alpha_values) <- c("sample_name", "richness","shannon","simpson", "inv_simpson", "sample_id", "treatment", "replicate", "timepoint", "treatment_replicate", "exp_group")
  alpha_values <- alpha_values %>% 
    select(sample_name, sample_id, treatment, replicate, treatment_replicate, timepoint, exp_group, richness, shannon, simpson, inv_simpson) %>%
    filter(timepoint==time)
  return(alpha_values)
}
alpha_values <- calc_vegan_table(mpa_mat, meta, "T_6")

####################################
# violin plot functions            #
####################################
#function to generate significance codes from p values
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

#function to perform wilcox test
calc_wilcox <- function(alphas, metric, group){
  res.stat <- alphas %>% 
    rstatix::wilcox_test(as.formula(paste(metric, group, sep = "~"))) %>%
    rstatix::adjust_pvalue(method = "BH") %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_xy_position(x=group, fun="max")
  return(res.stat)
}
#function to generate base violin plot
plot_base <- function(df, metric, group){
  p <- ggviolin(df, x=group, y=metric,
                fill=group, xlab="", ylab=metric, alpha=0.6)
  return(p)
}
#function to generate plot to extract the legend
plot_leg <- function(base_plot){
  p1 <- base_plot + 
    scale_fill_manual(name="Treatment group",
                      labels=c(expression(paste("No ", italic("S. salivarius"))),expression(paste("Any ", italic("S. salivarius")))),
                      values = pal1)+
    geom_point() +
    theme(legend.box = "horizontal", 
          legend.position = "bottom",
          legend.title = element_text(face="bold"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  return(p1)
}
#function to combine the base plot with the significance testing
plot_layer <- function(base_plot, stats){
  p1 <- base_plot + 
    stat_pvalue_manual(stats, label = "p.adj = {p.adj}", tip.length = 0.01)+
    scale_fill_manual(name="Treatment",
                      labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                      values = pal1)+
    geom_point() +
    theme(legend.position='none', 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"))
}
#function to generate the report of statistical testing
compute_stats <-function(wilcox_results){
  alpha_res <- c("MetaPhlAn4", wilcox_results$.y., paste(wilcox_results$group1," v.s ", wilcox_results$group2), "Wilcoxon rank sum",wilcox_results$statistic,wilcox_results$p.adj,"BH",wilcox_results$p.adj.signif)
}

#######################
# Main
#######################

#initialise output DF
alpha_all <- data.frame(Data="MetaPhlAn4", 
                        Metric="Metric",
                        Comparing="Comparing",
                        Test="Wilcox",
                        Test_Statistic="F",
                        P_value="P",
                        p_adjust_method="BH",
                        Significance_codes="NA")
#richness
richness <- plot_layer(plot_base(alpha_values, metric="richness", group = "treatment"),
                  calc_wilcox(alphas = alpha_values, metric = "richness", group="treatment"))
alpha_all <-rbind(alpha_all, compute_stats(calc_wilcox(alphas = alpha_values, metric = "richness", group="treatment")))
#shannon
shannon <- plot_layer(plot_base(alpha_values, metric="shannon", group = "treatment"),
                   calc_wilcox(alphas = alpha_values, metric = "shannon", group="treatment"))
alpha_all <-rbind(alpha_all, compute_stats(calc_wilcox(alphas = alpha_values, metric = "shannon", group="treatment")))
#simpson
simpson <- plot_layer(plot_base(alpha_values, metric="simpson", group = "treatment"),
                      calc_wilcox(alphas = alpha_values, metric = "simpson", group="treatment"))
alpha_all <-rbind(alpha_all, compute_stats(calc_wilcox(alphas = alpha_values, metric = "simpson", group="treatment")))
#inv_simpson
invsimpson <- plot_layer(plot_base(alpha_values, metric="inv_simpson", group = "treatment"),
                      calc_wilcox(alphas = alpha_values, metric = "inv_simpson", group="treatment"))
alpha_all <-rbind(alpha_all, compute_stats(calc_wilcox(alphas = alpha_values, metric = "inv_simpson", group="treatment")))

#######################
# generate final plot #
#######################
plots <- plot_grid(richness,shannon,simpson,invsimpson, nrow=1, align='v',  
                   label_size = 14)
#make title
title <- ggdraw() + draw_label("Alpha diversity between groups at T_6", fontface='bold')
plots2 <- plot_grid(title, plots, ncol=1, rel_heights=c(0.1, 1)) 
#combine wth legend pulled from base plot
leggo <- get_legend(plot_leg(richness))
legend <- get_legend(plot_base(alpha_values, metric="richness", group = "treatment"))


p_all <- ggdraw(plot_grid(plots2, leggo, 
                          nrow=2, align='v', rel_heights=c(1, 0.2),
                          rel_widths=c(1, 2))) + 
  theme(plot.background = element_rect(fill="#FFFFFF", color = NA),
        axis.title.y = element_text(size=12, face="bold", colour = "black"))
p_all
#ggsave("figures\\alpha_violins_T_6.emf", width = 4, height = 4, dpi = 300, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
#ggsave("figures\\alpha_violins_T_6.png", p_all, dpi=600)
ggsave("figures\\alpha_violins_T_6.pdf", p_all, height = 8, width = 12, units = "in", dpi = 150)
#outputstats
alpha_all<- alpha_all[-1,]
write.csv(alpha_all, "reports\\alpha_tested_T_6.csv", row.names = FALSE)


