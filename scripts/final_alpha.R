#final alpha

#' Script to calculate alpha diversity on NisinG samples
#' Developed by A. Kate Falà (https://github.com/aforestsomewhere)
#' TESTING FOR DIFFERENCE BY TREATMENT OR TIMEPOINT

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
#pal1 <- c("#072e51","#ffdf66","#332211", "#ffcc44", "#44aacc", "#bb2211")
palette_time <- c("#9A8822", "#F8AFA8", "#74A089")
palette_groups <- c("#332211","#ffcc44", "#44aacc", "#bb2211")
palette_collapsed <- c("#FF0000", "#00A08A")

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
meta <- meta %>%
  mutate(strep_exposure = case_when(
    treatment == "FS" ~ "Y", 
    treatment == "S" ~ "Y",
    treatment == "F" ~ "N",
    treatment == "B" ~ "N"))
#read MPA4 data
mpa_table <- read.table("cleaned_data/mpa_table.tsv", comment.char = '#', sep = '\t', header = TRUE, row.names=NULL)
#filter to species level assignments
mpa_table <- mpa_table[grep('s__',mpa_table[,1]),]
mpa_table[,1] <- gsub(".+\\|s__", "", mpa_table[,1])
mpa_table <- mpa_table %>% filter(!grepl("t__",row.names))
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]
#filter abundance data to timepoint T_0 only
samples <- rownames(meta)
mpa_table[,colnames(mpa_table) %in% samples] -> mpa_table
#generate matrix of abundance data
mpa_mat = t(as.matrix(mpa_table))
#order samples (rows) to properly merge with infO
mpa_mat <- mpa_mat[sort(rownames(mpa_mat)),] 

####################################
### calculate alpha diversity    ###
####################################
calc_vegan_table <- function(mpa, metadata){
  richness_vec <-(specnumber(mpa))
  shannon_vec <- vegan::diversity(mpa, index = "shannon")
  #simp_vec <- vegan::diversity(mpa, index = "simpson")
  inv_simp_vec <- vegan::diversity(mpa, index = "invsimpson")
  alpha_values <- cbind.data.frame(comb_id = rownames(meta), 
                                   richness = richness_vec, 
                                   shannon = shannon_vec, 
                                   inv_simpson = inv_simp_vec,
                                   meta=meta)
  #alpha_values <- merge(alpha_values, meta, by=row.names)
  names(alpha_values) <- c("sample_name", "richness","shannon","inv_simpson", "sample_id", "treatment", "replicate", "timepoint", "treatment_replicate", "exp_group", "strep_exposure")
  alpha_values <- alpha_values %>% 
    dplyr::select(-sample_name)
  return(alpha_values)
}
alpha_values <- calc_vegan_table(mpa_mat, meta)

#number of observed species
alpha_values %>%
  summarise(mean=mean(richness), sd=sd(richness), median=median(richness))

######################################
### Kruskal Wallis testing overall ###
######################################

#kruskal testing for differences between groups
#stat testing (uncollapsed)
indices <- c("richness","shannon","inv_simpson")
#can only test groups where n>5 - treatment, timepoint
results_time <- data.frame()
for (i in (1:length(indices))){
  result <- alpha_values %>%
    do(tidy(kruskal.test(x = .[[indices[i]]], g = .$timepoint)))
  result <- result %>% 
    mutate(metric=indices[i]) %>%
    mutate(group="timepoint") %>%
    mutate(sig_code=sigcode(p.value))
  results_time <- bind_rows(results_time, result)
}
results_treat <- data.frame()
for (i in (1:length(indices))){
  result <- alpha_values %>%
    do(tidy(kruskal.test(x = .[[indices[i]]], g = .$treatment)))
  result <- result %>% 
    mutate(metric=indices[i]) %>%
    mutate(group="treatment") %>%
    mutate(sig_code=sigcode(p.value))
  results_treat <- bind_rows(results_treat, result)
}
results_strep <- data.frame()
for (i in (1:length(indices))){
  result <- alpha_values %>%
    do(tidy(kruskal.test(x = .[[indices[i]]], g = .$strep_exposure)))
  result <- result %>% 
    mutate(metric=indices[i]) %>%
    mutate(group="strep_exposure") %>%
    mutate(sig_code=sigcode(p.value))
  results_strep <- bind_rows(results_strep, result)
}
kruskal_results <- rbind(results_time,results_treat, results_strep)

############################################
### Follow up test sig. from KW          ###
### Only timepoint -> Wilcoxon           ###
############################################
#only timepoint was significant - follow up with wilcoxon
alpha_time <- data.frame()
results4 <- data.frame()
#wilcoxon for each metric, pairwise between groups
for (i in 1:length(indices)){
  alpha_time <-rbind(alpha_time, calc_wilcox(alphas = alpha_values, metric = indices[i], group="timepoint"))
  results4 <- bind_rows(results4, alpha_time)
}
results_wilcox_alpha_time <- results4  %>%
  rename("diversity_metric" = 1) %>%
  select(1:8) %>%
  mutate(between_groups="timepoint") %>%
  mutate(sig_code=sigcode(as.numeric(p.adj)))

############################################
### Sig. differences between timepoints  ###
### Group within timepoints and test     ###
### Collapse for n=6                     ###
### Kruskal Wallis                       ###
############################################

results_exposure <- data.frame()
#testing for differences between strep_exposure Y/N
for (i in (1:length(indices))){
  result <- alpha_values %>%
    group_by(timepoint) %>%
    do(tidy(kruskal.test(x = .[[indices[i]]], g = .$strep_exposure)))
  result <- result %>% 
    mutate(metric=indices[i]) %>%
    mutate(sig_code=sigcode(p.value)) %>%
    mutate(group = "strep_exposure")
  results_exposure <- bind_rows(results_exposure, result)
}
#not significant differences

####################################
# violin plot functions            #
####################################

#function to perform wilcox test
calc_wilcox <- function(alphas, metric, group){
  res.stat <- alphas %>% 
    rstatix::wilcox_test(as.formula(paste(metric, group, sep = "~"))) %>%
    rstatix::adjust_pvalue(method = "BH") %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_xy_position(x=group, fun="max")
  stdevvar <- sd(alphas[,as.character(metric)])
  res.stat <- res.stat %>%
    mutate(y.position = y.position + (2*stdevvar)) #very manual hack to move p values up
  return(res.stat)
}

#calc_wilcox(alphas = alpha_values, metric = "richness", group="timepoint")
#function to generate base violin plot
plot_base <- function(df, metric, group){
  p <- ggviolin(df, x=group, y=metric,
                fill=group, xlab="", ylab=metric, alpha=0.6)
  return(p)
}
#function to generate plot to extract the legend
plot_leg <- function(base_plot,legendtitle, pal){
  p1 <- base_plot + 
    scale_fill_manual(name=legendtitle,
                      #labels=c(expression(paste("No ", italic("S. salivarius"))),expression(paste("Any ", italic("S. salivarius")))),
                      values = pal)+
    geom_point() +
    theme(legend.box = "horizontal", 
          legend.position = "bottom",
          legend.title = element_text(face="bold"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  return(p1)
}
#function to combine the base plot with the significance testing
plot_layer <- function(base_plot, stats, legendtitle, pal){
  p1 <- base_plot + 
    stat_pvalue_manual(stats, label = "p.adj = {p.adj}", tip.length = 0.01)+
    scale_fill_manual(name=legendtitle,
                      #labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                      values = pal)+
    geom_point() +
    theme(legend.position='none', 
          #axis.text.x = element_blank(), 
          #axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"))
}
#function to generate the report of statistical testing
# compute_stats <-function(wilcox_results){
#   alpha_res <- c("MetaPhlAn4", wilcox_results$.y., paste(wilcox_results$group1," v.s ", wilcox_results$group2), "Wilcoxon rank sum",wilcox_results$statistic,wilcox_results$p.adj,"BH",wilcox_results$p.adj.signif)
# }

#######################
# Main
#######################
#grouped by timepoint - wilcoxon - three timepoints, three possible comparisons...
#initialise output DF
alpha_all <- data.frame()
# alpha_all <- data.frame(Data="MetaPhlAn4", 
#                         Metric="Metric",
#                         Comparing="Comparing",
#                         Test="Wilcox",
#                         Test_Statistic="F",
#                         P_value="P",
#                         p_adjust_method="BH",
#                         Significance_codes="NA")
#richness x timepoint
richness <- plot_base(alpha_values, metric="richness", group = "timepoint") +
  scale_fill_manual(name="Timepoint",
                    #labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = palette_time) + 
  theme(axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position='none') +
  geom_point()#only plot_base - no sig differences by KW testing
richness
leg_time <- get_legend(plot_leg(plot_base(alpha_values, metric="richness", group = "timepoint"), legendtitle = "Timepoint", pal=palette_time))
#alpha_all <-rbind(alpha_all, (calc_wilcox(alphas = alpha_values, metric = "richness", group="timepoint")))
#shannon x timepoint
shannon <- plot_layer(plot_base(alpha_values, metric="shannon", group = "timepoint"),
                      calc_wilcox(alphas = alpha_values, metric = "shannon", group="timepoint"), 
                      legendtitle = "Timepoint", pal=palette_time)
shannon
alpha_all <-rbind(alpha_all, (calc_wilcox(alphas = alpha_values, metric = "shannon", group="timepoint")))
#inverse simpson x timepoint
invsimpson <- plot_layer(plot_base(alpha_values, metric="inv_simpson", group = "timepoint"),
                         calc_wilcox(alphas = alpha_values, metric = "inv_simpson", group="timepoint"),
                         legendtitle = "Timepoint", pal=palette_time)
invsimpson
alpha_all <-rbind(alpha_all, (calc_wilcox(alphas = alpha_values, metric = "inv_simpson", group="timepoint")))

#plot together
plots <- plot_grid(richness,shannon,invsimpson, nrow=3, align='v',  
                   label_size = 14)
#make title
title <- ggdraw() + draw_label("Alpha diversity between timepoints", fontface='bold')
plots2 <- plot_grid(title, plots, ncol=1, rel_heights=c(0.1, 1)) 
#combine wth legend pulled from base plot
#leggo <- get_legend(plot_leg(richness))
#legend <- get_legend(plot_base(alpha_values, metric="richness", group = "treatment"))
p_all <- ggdraw(plot_grid(plots2, leg_time, 
                          nrow=2, align='v', rel_heights=c(1, 0.1),
                          rel_widths=c(1, 2))) + 
  theme(plot.background = element_rect(fill="#FFFFFF", color = NA),
        axis.title.y = element_text(size=12, face="bold", colour = "black"))
p_all
ggsave("figures\\alpha_violins_between_times.pdf", p_all, height = 40, width = 12, units = "cm", dpi = 300)
#outputstats
alpha_all <- as.data.frame(alpha_all)
alpha_all <- apply(alpha_all,2,as.character)
write.csv(alpha_all, "reports\\alpha_tested_between_times_allgroups.csv")

#######################
# All treatment - no stats as KW insig
#######################
richness <- plot_base(alpha_values, metric="richness", group = "treatment") +
  scale_fill_manual(name="Treatment",
                    #labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = palette_groups) + 
  theme(axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position='none') +
  geom_point()#only plot_base - no sig differences by KW testing
richness
leg_treat <- get_legend(plot_leg(plot_base(alpha_values, metric="richness", group = "treatment"), legendtitle = "Treatment", pal=palette_groups))


#check <- compute_stats(calc_wilcox(alphas = alpha_values, metric = "richness", group="timepoint"))
#shannon
shannon <- plot_base(alpha_values, metric="shannon", group = "treatment") +
  scale_fill_manual(name="Treatment",
                    #labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = palette_groups) + 
  theme(axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position='none') +
  geom_point()#only plot_base - no sig differences by KW testing
shannon
#alpha_all <-rbind(alpha_all, compute_stats(calc_wilcox(alphas = alpha_values, metric = "shannon", group="treatment")),
#simpson
# simpson <- plot_layer(plot_base(alpha_values, metric="simpson", group = "treatment"),
#                       calc_wilcox(alphas = alpha_values, metric = "simpson", group="treatment"))
# alpha_all <-rbind(alpha_all, compute_stats(calc_wilcox(alphas = alpha_values, metric = "simpson", group="treatment")))
#inv_simpson
invsimpson <- plot_base(alpha_values, metric="inv_simpson", group = "treatment") +
  scale_fill_manual(name="Treatment",
                    #labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = palette_groups) + 
  theme(axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position='none') +
  geom_point()#only plot_base - no sig differences by KW testing
invsimpson


#plot together
plots <- plot_grid(richness,shannon,invsimpson, nrow=3, align='v',  
                   label_size = 14)
#make title
title <- ggdraw() + draw_label("Alpha diversity between treatment groups", fontface='bold')
plots2 <- plot_grid(title, plots, ncol=1, rel_heights=c(0.1, 1)) 
#combine wth legend pulled from base plot
#leggo <- get_legend(plot_leg(richness))
#legend <- get_legend(plot_base(alpha_values, metric="richness", group = "treatment"))
p_all <- ggdraw(plot_grid(plots2, leg_treat, 
                          nrow=2, align='v', rel_heights=c(1, 0.1),
                          rel_widths=c(1, 2))) + 
  theme(plot.background = element_rect(fill="#FFFFFF", color = NA),
        axis.title.y = element_text(size=12, face="bold", colour = "black"))
p_all
ggsave("figures\\alpha_violins_between_groups.pdf", p_all, height = 40, width = 12, units = "cm", dpi = 300)

################################################
# All strep exposed Y/N - no stats as KW insig #
################################################
richness <- plot_base(alpha_values, metric="richness", group = "strep_exposure") +
  scale_fill_manual(name="Any exposure",
                    #labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = palette_collapsed) + 
  theme(axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position='none') +
  geom_point()#only plot_base - no sig differences by KW testing
richness
leg_strep <- get_legend(plot_leg(plot_base(alpha_values, metric="richness", group = "strep_exposure"), legendtitle = "Any exposure", pal=palette_collapsed))

#shannon
shannon <- plot_base(alpha_values, metric="shannon", group = "strep_exposure") +
  scale_fill_manual(name="Any exposure",
                    values = palette_collapsed) + 
  theme(axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position='none') +
  geom_point()#only plot_base - no sig differences by KW testing
shannon

#inv_simpson
invsimpson <- plot_base(alpha_values, metric="inv_simpson", group = "strep_exposure") +
  scale_fill_manual(name="Any exposure",
                    #labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = palette_collapsed) + 
  theme(axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position='none') +
  geom_point()#only plot_base - no sig differences by KW testing
invsimpson


#plot together
plots <- plot_grid(richness,shannon,invsimpson, nrow=3, align='v',  
                   label_size = 14)
#make title
title <- ggdraw() + draw_label("Alpha diversity between any exposure", fontface='bold')
plots2 <- plot_grid(title, plots, ncol=1, rel_heights=c(0.1, 1)) 
#combine wth legend pulled from base plot
#leggo <- get_legend(plot_leg(richness))
#legend <- get_legend(plot_base(alpha_values, metric="richness", group = "treatment"))
p_all <- ggdraw(plot_grid(plots2, leg_strep, 
                          nrow=2, align='v', rel_heights=c(1, 0.1),
                          rel_widths=c(1, 2))) + 
  theme(plot.background = element_rect(fill="#FFFFFF", color = NA),
        axis.title.y = element_text(size=12, face="bold", colour = "black"))
p_all
ggsave("figures\\alpha_violins_between_strep_exposure.pdf", p_all, height = 40, width = 12, units = "cm", dpi = 300)




#######################
# generate final plot #
#######################
# plots <- plot_grid(richness,shannon,simpson,invsimpson, nrow=1, align='v',  
#                    label_size = 14)
# #make title
# title <- ggdraw() + draw_label("Alpha diversity between groups at T_0", fontface='bold')
# plots2 <- plot_grid(title, plots, ncol=1, rel_heights=c(0.1, 1)) 
# #combine wth legend pulled from base plot
# leggo <- get_legend(plot_leg(richness))
# legend <- get_legend(plot_base(alpha_values, metric="richness", group = "treatment"))
# 
# 
# p_all <- ggdraw(plot_grid(plots2, leggo, 
#                           nrow=2, align='v', rel_heights=c(1, 0.2),
#                           rel_widths=c(1, 2))) + 
#   theme(plot.background = element_rect(fill="#FFFFFF", color = NA),
#         axis.title.y = element_text(size=12, face="bold", colour = "black"))
# p_all
# #ggsave("figures\\alpha_violins_T_0.emf", width = 4, height = 4, dpi = 300, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
# #ggsave("figures\\alpha_violins_T_0.png", p_all, dpi=600)
# ggsave("figures\\alpha_violins_T_0.pdf", p_all, height = 8, width = 12, units = "in", dpi = 150)
# #outputstats
# alpha_all<- alpha_all[-1,]
# write.csv(alpha_all, "reports\\alpha_tested_T_0.csv", row.names = FALSE)

#######################
#generate output file #
#######################
list_of_datasets <- list("Alpha diversity metrics" = alpha_values,
                         "Kruskal Wallis, overall groups" = kruskal_results,
                         "Wilcoxon, timepoint" = results_wilcox_alpha_time,
                         "Kruskal Wallis, collapsed" = results_exposure)
openxlsx::write.xlsx(list_of_datasets, file = "reports/alpha_diversity.xlsx")


