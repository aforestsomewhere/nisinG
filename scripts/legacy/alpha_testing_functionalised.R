#' Script to calculate alpha diversity on NisinG samples
#' Developed by A. Kate Falà (https://github.com/aforestsomewhere)
#' collapsing any S. salivarius-treated; 

####################
### Packages     ###
####################
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
#install.packages("showtext")
#library('devtools')
#theme_name = "theme_starwars" # Pick which theme you want
#theme_url = paste0("https://raw.githubusercontent.com/MatthewBJane/theme_park/main/", theme_name ,".R")
#source_url(theme_url)

################
# FUNCTIONS    #
################
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
rm()
dev.off()
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
pal1 <- c("#072e51",  "#ffdf66","#332211", "#ffcc44", "#44aacc", "#bb2211")
########################
#   Richness boxplot   #
########################
#initialise output DF
alpha_all <- data.frame(Data="MetaPhlAn4", 
                        Metric="Metric",
                        Comparing="Comparing",
                        Test="Wilcox",
                        Test_Statistic="F",
                        P_value="P",
                        p_adjust_method="BH",
                        Significance_codes="NA")
#function to perform wilcox test
calc_wilcox <- function(alphas, metric, group){
  res.stat <- alphas %>% 
    rstatix::wilcox_test(as.formula(paste(metric, group, sep = "~"))) %>%
    rstatix::adjust_pvalue(method = "BH") %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_xy_position(x = "metric", fun = "max")
  return(res.stat)
}
#calc_wilcox(alphas = alpha_values, metric = "richness", group="treatment")
plot_base <- function(df, metric, group){
  p <- ggviolin(df, x=group, y=metric,
                fill=group, xlab="", ylab=metric, alpha=0.6)
  return(p)
}
#layered plot function
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

#try together for richness

#ptest <- plot_layer(vio, res.stat)

#ptest
all <- plot_layer(plot_base(alpha_values, metric="richness", group = "treatment"),
                  calc_wilcox(alphas = alpha_values, metric = "richness", group="treatment"))
all

#richness
#base plot
vio <- plot_base(alpha_values, metric="richness", group = "treatment")
#stat test
res.stat <- calc_wilcox(alphas = alpha_values, metric = "richness", group="treatment")

alpha_res <- c("MetaPhlAn4", res.stat$.y., paste(res.stat$group1," v.s ", res.stat$group2), "Wilcoxon rank sum",res.stat$statistic,res.stat$p.adj,"BH",res.stat$p.adj.signif)
alpha_all <-rbind(alpha_all, alpha_res)

p1<- vio + stat_pvalue_manual(res.stat, label = "p.adj = {p.adj}", tip.length = 0.01)+
  scale_fill_manual(name="Treatment",
                    labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = pal1)+
  geom_point() +
  theme(legend.position='none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"))
p1





#richness_vec <-(specnumber(mpa_mat))
#shannon_vec <- vegan::diversity(mpa_mat, index = "shannon")
#simp_vec <- vegan::diversity(mpa_mat, index = "simpson")
#inv_simp_vec <- vegan::diversity(mpa_mat, index = "invsimpson")
# alpha_values <- cbind.data.frame(comb_id = rownames(meta), 
#                                  richness = richness_vec, 
#                                  shannon = shannon_vec, 
#                                  simpson = simp_vec, 
#                                  inv_simpson = inv_simp_vec, 
#                                  meta = meta)
# rm(richness_vec, shannon_vec, simp_vec, inv_simp_vec)

#alpha_df$meta.treatment = factor(alpha_df$meta.treatment)
#Rationalise column names
# names(alpha_values) <- c("sample_name", "richness","shannon","simpson", "inv_simpson", "sample_id", "treatment", "replicate", "timepoint", "treatment_replicate", "exp_group")
# alpha_values <- alpha_values %>% 
#   select(sample_name, sample_id, treatment, replicate, treatment_replicate, timepoint, exp_group, richness, shannon, simpson, inv_simpson) %>%
#   filter(timepoint=="T_6")
#write.csv(alpha_df, "reports/alpha_all.csv", row.names = FALSE)

################################################
### test differences in alpha diversity      ###
### due to sample size n=3, no testing for   ###
### individual groups / timepoints           ###
### violin plots to show points              ###
################################################
theme_set(theme_bw())
#first plto to extract common legend
pal1 <- c("#072e51",  "#ffdf66","#332211", "#ffcc44", "#44aacc", "#bb2211")
legend <- ggviolin(alpha_values, 
                   x = "treatment", 
                   y = "richness",
                   fill = "treatment", 
                   xlab="", 
                   ylab="Observed Species", 
                   alpha=0.6)  +
  scale_fill_manual(name="Treatment group",
                    labels=c('No exposure', 
                             expression(paste(italic("Any S. salivarius")))), 
                    values = pal1) + 
  xlab("") + ylab("") +
  theme(legend.box = "horizontal", 
        legend.position = "bottom",
        axis.ticks.x = element_blank())
legend <- get_legend(legend)

########################
#   Richness boxplot   #
########################
#base plot
vio <- ggviolin(alpha_values, x = "treatment", y = "richness", 
                fill = "treatment", xlab="", ylab="Observed Species", alpha=0.6)
#perform stat testing
res.stat <- alpha_values %>% 
  rstatix::wilcox_test(data = ., richness ~ treatment) %>% 
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "treatment", fun = "max")
#compile results to report
alpha_all <- data.frame(Data="MetaPhlAn4", 
                        Metric="Metric",
                        Comparing="Comparing",
                        Test="Wilcox",
                        Test_Statistic="F",
                        P_value="P",
                        p_adjust_method="BH",
                        Significance_codes="NA") #initialise output DF
alpha_res <- c("MetaPhlAn4", res.stat$.y., paste(res.stat$group1," v.s ", res.stat$group2), "Wilcoxon rank sum",res.stat$statistic,res.stat$p.adj,"BH",res.stat$p.adj.signif)
alpha_all <-rbind(alpha_all, alpha_res)
p1<- vio + stat_pvalue_manual(res.stat, label = "p.adj = {p.adj}", tip.length = 0.01)+
  scale_fill_manual(name="Treatment",
                    labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = pal1)+
  geom_point() +
  theme(legend.position='none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"))
p1

########################
#   Shannon boxplot    #
########################

#base plot
vio <- ggviolin(alpha_values, x = "treatment", y = "shannon", 
                fill = "treatment", xlab="", ylab="Shannon's Index", alpha=0.6)
#perform stat testing
res.stat <- alpha_values %>% 
  rstatix::wilcox_test(data = ., shannon ~ treatment) %>% 
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "treatment", fun = "max")

#compile results to report
alpha_df <- c("MetaPhlAn4", res.stat$.y., paste(res.stat$group1," v.s ",res.stat$group2), "Wilcoxon rank sum",res.stat$statistic,res.stat$p.adj,"BH",res.stat$p.adj.signif)
alpha_all <-rbind(alpha_all, alpha_df)
p2<- vio + stat_pvalue_manual(res.stat, label = "p.adj = {p.adj}", tip.length = 0.01)+
  scale_fill_manual(name="Treatment",
                    labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = pal1)+  
  geom_point() + theme(legend.position='none', 
                       axis.text.x = element_blank(), 
                       axis.ticks.x = element_blank(),
                       axis.title.y = element_text(size=12, face="bold", colour = "black"))
p2

########################
#   Simpsons boxplot   #
########################

#base plot
vio <- ggviolin(alpha_values, x = "treatment", y = "simpson", 
                fill = "treatment", xlab="", ylab="Simpson's Index", alpha=0.6)

#perform stat testing
res.stat <- alpha_values %>% 
  rstatix::wilcox_test(data = ., simpson ~ treatment) %>% 
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "treatment", fun = "max")

#compile results to report
alpha_df <- c("MetaPhlAn4", res.stat$.y., paste(res.stat$group1," v.s ",res.stat$group2), "Wilcoxon rank sum",res.stat$statistic,res.stat$p.adj,"BH",res.stat$p.adj.signif)
alpha_all <-rbind(alpha_all, alpha_df)
p3<- vio + stat_pvalue_manual(res.stat, label = "p.adj = {p.adj}", tip.length = 0.01)+
  scale_fill_manual(name="Treatment",
                    labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = pal1)+  
  geom_point() + theme(legend.position='none', 
                       axis.text.x = element_blank(), 
                       axis.ticks.x = element_blank(),
                       axis.title.y = element_text(size=12, face="bold", colour = "black"))
p3

########################
#   Simpsons boxplot   #
########################

#base plot
vio <- ggviolin(alpha_values, x = "treatment", y = "inv_simpson", 
                fill = "treatment", xlab="", ylab="Inverse Simpson's Index", alpha=0.6)

#perform stat testing
res.stat <- alpha_values %>% 
  rstatix::wilcox_test(data = ., inv_simpson ~ treatment) %>% 
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "treatment", fun = "max")

#compile results to report
alpha_df <- c("MetaPhlAn4", res.stat$.y., paste(res.stat$group1," v.s ",res.stat$group2), "Wilcoxon rank sum",res.stat$statistic,res.stat$p.adj,"BH",res.stat$p.adj.signif)
alpha_all <-rbind(alpha_all, alpha_df)
p4<- vio + stat_pvalue_manual(res.stat, label = "p.adj = {p.adj}", tip.length = 0.01)+
  scale_fill_manual(name="Treatment",
                    labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = pal1)+  
  geom_point() + theme(legend.position='none', 
                       axis.text.x = element_blank(), 
                       axis.ticks.x = element_blank(),
                       axis.title.y = element_text(size=12, face="bold", colour = "black"))
p4

#######################
# generate final plot #
#######################
plots <- plot_grid(p1,p2,p3,p4, nrow=1, align='v',  
                   label_size = 14)
#make title
title <- ggdraw() + draw_label("Alpha diversity between groups at T_6", fontface='bold')
plots2 <- plot_grid(title, plots, ncol=1, rel_heights=c(0.1, 1)) 
#combine
p_all <- ggdraw(plot_grid(plots2, legend, 
                          nrow=2, align='v', rel_heights=c(1, 0.2),
                          rel_widths=c(1, 2))) + 
  theme(plot.background = element_rect(fill="#FFFFFF", color = NA),
        axis.title.y = element_text(size=12, face="bold", colour = "black"))
#ggsave("figures\\alpha_violins_T_6.emf", width = 4, height = 4, dpi = 300, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
#ggsave("figures\\alpha_violins_T_6.png", p_all, dpi=600)
ggsave("figures\\alpha_violins_T_6.pdf", p_all, height = 8, width = 12, units = "in", dpi = 150)
#outputstats
alpha_all<- alpha_all[-1,]
write.csv(alpha_all, "reports\\alpha_tested_T_6.csv", row.names = FALSE)


















legend <-ggplot(data = alpha_df, aes(x=treatment, y=richness, fill=treatment))+
  geom_violin() +geom_point(color="black", alpha=0.8)  +
  scale_fill_manual(name="Treatment group",
                    labels=c('Control', 
                             expression(paste(italic("F. nucleatum"))), 
                             expression(paste(italic("F. nucleatum/S. salivarius"))), 
                             expression(paste(italic('S. salivarius')))),
                    values = pal1) + 
  xlab("") + ylab("Observed species") +
  theme(legend.box = "horizontal", 
        legend.position = "bottom",
        legend.title = element_text(face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + facet_wrap(~timepoint)
legend <- get_legend(legend)
#Richness
p1 <-ggplot(data = alpha_df, aes(x=treatment, y=richness, fill=treatment))+
  geom_violin() +geom_point(color="black", alpha=0.8)  +
  scale_fill_manual(name="Treatment Group",
                    labels=c('Control', 
                             expression(paste(italic("F. nucleatum"))), 
                             expression(paste(italic("F. nucleatum/S. salivarius"))), 
                             expression(paste(italic('S. salivarius')))),
                    values = pal1) + 
  xlab("") + ylab("Observed species") +
  theme(legend.box = "vertical", 
        legend.position = "right",
        legend.title = element_text(face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        strip.text.x = element_text(size = 12, face="bold", colour = "black" ),
        strip.text.y = element_text(size = 12, face="bold", colour = "black")) + facet_wrap(~timepoint) +
  guides(fill="none", alpha="none", color="none")
p1
#Shannon
p2 <-ggplot(data = alpha_df, aes(x=treatment, y=shannon, fill=treatment))+
  geom_violin() +geom_point(color="black", alpha=0.8)  +
  scale_fill_manual(name="Treatment Group",
                    labels=c('Control', 
                             expression(paste(italic("F. nucleatum"))), 
                             expression(paste(italic("F. nucleatum/S. salivarius"))), 
                             expression(paste(italic('S. salivarius')))),
                    values = pal1) + 
  xlab("") + ylab("Shannon's Index") +
  theme(legend.box = "vertical", 
        legend.position = "right",
        legend.title = element_text(face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        strip.text.x = element_text(size = 12, face="bold", colour = "black" ),
        strip.text.y = element_text(size = 12, face="bold", colour = "black")) + facet_wrap(~timepoint) +
  guides(fill="none", alpha="none", color="none")
p2
#Simpson
p3 <-ggplot(data = alpha_df, aes(x=treatment, y=simpson, fill=treatment))+
  geom_violin() +geom_point(color="black", alpha=0.8)  +
  scale_fill_manual(name="Treatment Group",
                    labels=c('Control', 
                             expression(paste(italic("F. nucleatum"))), 
                             expression(paste(italic("F. nucleatum/S. salivarius"))), 
                             expression(paste(italic('S. salivarius')))),
                    values = pal1) + 
  xlab("") + ylab("Simpson's Index") +
  theme(legend.box = "vertical", 
        legend.position = "right",
        legend.title = element_text(face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        strip.text.x = element_text(size = 12, face="bold", colour = "black" ),
        strip.text.y = element_text(size = 12, face="bold", colour = "black")) + facet_wrap(~timepoint) +
  guides(fill="none", alpha="none", color="none")
p3
#Inverse Simpson
p4 <-ggplot(data = alpha_df, aes(x=treatment, y=inv_simpson, fill=treatment))+
  geom_violin() +geom_point(color="black", alpha=0.8)  +
  scale_fill_manual(name="Treatment Group",
                    labels=c('Control', 
                             expression(paste(italic("F. nucleatum"))), 
                             expression(paste(italic("F. nucleatum/S. salivarius"))), 
                             expression(paste(italic('S. salivarius')))),
                    values = pal1) + 
  xlab("") + ylab("Inverse Simpson's Index") +
  theme(legend.box = "vertical", 
        legend.position = "right",
        legend.title = element_text(face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        strip.text.x = element_text(size = 12, face="bold", colour = "black" ),
        strip.text.y = element_text(size = 12, face="bold", colour = "black")) + facet_wrap(~timepoint) +
  guides(fill="none", alpha="none", color="none")
p4

all_alpha <- plot_grid(p1,p2,p3,p4, nrow=2)

all_alpha_leg <- plot_grid(all_alpha, legend, 
                           nrow=2, align='v', rel_heights=c(1, 0.2),
                           rel_widths=c(1, 2)) + theme(plot.background = element_rect(fill="#FFFFFF", color = NA))
all_alpha_leg
ggsave("figures/alpha_all.pdf", width=10, height=10, dpi=300, units="cm")
ggsave("figures/alpha_all.emf", width = 10, height = 10, dpi = 600, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
