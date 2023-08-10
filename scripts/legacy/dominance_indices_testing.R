#' Script to calculate dominance indices for NisinG samples
#' Developed by A. Kate Falà (https://github.com/aforestsomewhere)
#' Compute raw dominance indices
#' For stat testing need to collapse B with F and FS with s
#
###############################
#'Dominant taxa
#'First done with species level assignment,
#'Then with full taxonomic ordination
#'Dominance indices are in general negatively correlated with diversity, and 
#'sometimes used in ecological literature. High dominance is obtained when one 
#'or few species have a high share of the total species abundance in the community.
#'MPA4
#'13/03/23
#'NisinG
#'Katie F
###############################
palette_time <- c("#9A8822", "#F8AFA8", "#74A089")
palette_groups <- c("#332211","#ffcc44", "#44aacc", "#bb2211")
palette_collapsed <- c("#FF0000", "#00A08A")
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
### Packages     ###
####################
library(dplyr)
library(rstatix)
library(tidyr)
library(tidyverse)
library(vegan)
library(phyloseq)
library(FSA)
library(openxlsx)
theme_set(theme_classic())

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
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]

#########################
### form Phyloseq obj ###
#########################
#generate base obj
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


# mpa.phy = metaphlanToPhyloseq(mpa_table)
# phy <- metaphlanToPhyloseq(
#   phyloseq::otu_table(mpa.phy),
#   metadat=phyloseq::sample_data(meta, errorIfNULL = TRUE),
#   simplenames=TRUE,
#   roundtointeger=FALSE)

#dominance
#make relative first
#microbiome::aggregate_rare(phy, level = 'Kingdom',
#  detection = 0.1/100, prevalence = 5/100) -> phy_agg
microbiome::transform(phy, 'compositional') -> phy2
dom <- cbind(phy2@sam_data, microbiome::dominant(phy2, level = NULL))
dom$sample_name <- rownames(dom)
microbiome::dominance(phy2, index = "all", rank = 1, relative = TRUE, aggregate = TRUE) -> dom2
dom2$sample_name <- rownames(dom2)
dom_big <- dplyr::left_join(dom2, dom, by="sample_name")
names(dom_big)[names(dom_big) == 'microbiome::dominant(phy2, level = NULL)'] <- "dominant_species"

#stat testing (uncollapsed)
subbed <- dom_big[,1:7]
indices <- names(subbed)
#can only test groups where n>5 - treatment, timepoint
#Kruskal - don't need to adjust (not pairwise - Dunn/Wilcoxon yes)
results <- data.frame()
for (i in (1:length(indices))){
    result <- dom_big %>%
      do(tidy(kruskal.test(x = .[[indices[i]]], g = .$treatment)))
    result <- result %>% 
      mutate(metric=indices[i]) %>%
      mutate(group="treatment") %>%
      mutate(sig_code=sigcode(p.value))
    results <- bind_rows(results, result)
}
results2 <- data.frame()
for (i in (1:length(indices))){
  result <- dom_big %>%
    do(tidy(kruskal.test(x = .[[indices[i]]], g = .$timepoint)))
  result <- result %>% 
    mutate(metric=indices[i]) %>%
    mutate(group="timepoint") %>%
    mutate(sig_code=sigcode(p.value))
  results2 <- bind_rows(results2, result)
}
results3 <- data.frame()
for (i in (1:length(indices))){
  result <- dom_big %>%
    do(tidy(kruskal.test(x = .[[indices[i]]], g = .$strep_exposure)))
  result <- result %>% 
    mutate(metric=indices[i]) %>%
    mutate(group="strep_exposure") %>%
    mutate(sig_code=sigcode(p.value))
  results3 <- bind_rows(results3, result)
}
results_uncollapsed <- rbind(results, results2, results3)
#only timepoint was significant
c <- dunnTest(dbp ~ timepoint, data = dom_big, method="bh")
dunn_result <- c$res %>%
  mutate(method="Dunn's Test") %>%
  mutate(correction="BH") %>%
  rowwise() %>%
  mutate(sig_code=sigcode(as.numeric(P.adj)))

#####################
# Plots for KW/Dunn #
#####################
#timepoint sig - follow up
#function to perform wilcox test
calc_wilcox <- function(dom, metric, group){
  res.stat <- dom %>% 
    rstatix::wilcox_test(as.formula(paste(metric, group, sep = "~"))) %>%
    rstatix::adjust_pvalue(method = "BH") %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_xy_position(x=group, fun="max")
  stdevvar <- sd(dom[,as.character(metric)])
  res.stat <- res.stat %>%
    mutate(y.position = y.position + (2*stdevvar)) #very manual hack to move p values up
  return(res.stat)
}
#calc_wilcox(dom = dom_big, metric = "dbp", group="timepoint")
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
    theme(legend.box = "vertical", 
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

#######################
# Main
#######################
alpha_all <- data.frame()
#get list of dominance indices
#stat testing
subbed <- dom_big[,1:7]
metric_list <-  names(subbed)
rm(subbed)

#timepoint sig different - pairwise wilcoxon to ascertain which different, violin plots to show
#make plot for each dominance index (all sig dif)
myplots <- vector("list", 7)
#get legend for group plot
legend <- get_legend(plot_leg(plot_base(df=dom_big, metric = metric_list[i], group="timepoint"), 
                              legendtitle="Timepoint", pal=palette_time) +
                       theme(legend.box = "vertical", 
                             legend.position = "bottom",
                             legend.title = element_text(face="bold")))
for (i in 1:length(metric_list)){
  myplots[[i]] <- local({
    i <- i
    p1 <- plot_layer(plot_base(df=dom_big, metric = metric_list[i], group="timepoint"),
                     calc_wilcox(dom=dom_big, metric = metric_list[i], group="timepoint"),
                               legendtitle="Timepoint", pal=palette_time)
    print(p1)
  })
}
final_plot <- plot_grid(ncol=3, 
                        plotlist = myplots, legend, 
                        labels = "AUTO")
final_plot
ggsave("figures\\dom_violins_between_timepoints.pdf", final_plot, height = 29, width = 21, units = "cm", dpi = 300)

#treatment not sig different - violin plots to show
myplots <- vector("list", 7)
#get legend for group plot
newlegend <- get_legend(plot_leg(plot_base(df=dom_big, metric = metric_list[i], group="treatment"), 
                              legendtitle="Treatment", pal=palette_groups) +
                       theme(legend.direction = "vertical", 
                             legend.position = "bottom",
                             legend.title = element_text(face="bold")))
for (i in 1:length(metric_list)){
  myplots[[i]] <- local({
    i <- i
    p1 <- plot_base(df=dom_big, metric = metric_list[i], group="treatment") + 
      scale_fill_manual(name="Treatment group", values = palette_groups) +
      theme(axis.title.y = element_text(size=12, face="bold", colour = "black"),
            legend.position='none') + geom_point()
    print(p1)
  })
}
final_plot <- plot_grid(ncol=3, 
                        plotlist = myplots, newlegend, 
                        labels = "AUTO")
final_plot
ggsave("figures\\dom_violins_between_treatments.pdf", final_plot, height = 29, width = 21, units = "cm", dpi = 300)

#Strep exposure not significant # violin plots to show
myplots <- vector("list", 7)
#get legend for group plot
legend <- get_legend(plot_leg(plot_base(df=dom_big, metric = metric_list[i], group="strep_exposure"), 
                                 legendtitle="Any exposure", pal=palette_collapsed) +
                          theme(legend.direction = "vertical", 
                                legend.position = "bottom",
                                legend.title = element_text(face="bold")))
for (i in 1:length(metric_list)){
  myplots[[i]] <- local({
    i <- i
    p1 <- plot_base(df=dom_big, metric = metric_list[i], group="strep_exposure") + 
      scale_fill_manual(name="Any exposure", values = palette_collapsed) +
      theme(axis.title.y = element_text(size=12, face="bold", colour = "black"),
            legend.position='none') + geom_point()
    print(p1)
  })
}
final_plot <- plot_grid(ncol=3, 
                        plotlist = myplots, legend, 
                        labels = "AUTO")
final_plot
ggsave("figures\\dom_violins_between_exposure.pdf", final_plot, height = 29, width = 21, units = "cm", dpi = 300)

##########################
# Final test - collapsed to exposure, but grouped within timepoints
#grouping_var <- c("strep_exposure", "dominant_species")
results_collapsed <- data.frame()
#testing for differences between strep_exposure Y/N
for (i in (1:length(metric_list))){
    result <- dom_big %>%
      group_by(timepoint) %>%
      do(tidy(kruskal.test(x = .[[metric_list[i]]], g = .$strep_exposure)))
    result <- result %>% 
      mutate(metric=metric_list[i]) %>%
      mutate(sig_code=sigcode(p.value)) %>%
      mutate(group = "strep_exposure")
    results_collapsed <- bind_rows(results_collapsed, result)
}
# Final test - collapsed to dominant species, where more than 2 available, but only T_24
for (i in 1:length(metric_list)) {
  result <- dom_big %>%
    #filter(dominant_species != "Bifidobacterium_bifidum") %>%
    #group_by(timepoint) %>%
    summarise(p.value = ifelse(n_distinct(dominant_species) >= 2,
                               kruskal.test(get(metric_list[i]), dominant_species)$p.value,
                               NA)) %>%
    mutate(metric = metric_list[i],
           sig_code = case_when(!is.na(p.value) & p.value < 0.05 ~ "*", 
                                !is.na(p.value) & p.value >= 0.05 ~ "",
                                TRUE ~ NA_character_),
           group = "dominant_species",
           method = "Kruskal-Wallis rank sum test")
  results_collapsed <- bind_rows(results_collapsed, result)
}


# for (i in (1:length(metric_list))){
#   result <- dom_big %>%
#     #filter(dominant_species!="Bifidobacterium_bifidum") %>%
#     group_by(timepoint) %>%
#     if(nlevels(dominant_species) >=2){
#       do(tidy(kruskal.test(x = .[[metric_list[i]]], g = .$dominant_species)))
#     }
#   result <- result %>% 
#     mutate(metric=metric_list[i]) %>%
#     mutate(sig_code=sigcode(p.value)) %>%
#     mutate(group = "dominant_species")
#   results_collapsed <- bind_rows(results_collapsed, result)
# }

#Follow up testing: only T24, only DBP, absolute, relative, both strep_exposure and dominant species
# Final test - collapsed to dominant species, where more than 2 available, but only T_24
dom_24 <- dom_big %>% filter(timepoint=="T_24")
indices <- c("dbp", "absolute", "relative")
groups <- c("strep_exposure", "dominant_species")

#wilcoxon?
#function to perform wilcox test
#calc_wilcox(dom = dom_24, metric = "dbp", group="strep_exposure")
leg_strep <- plot_leg(plot_base(df=dom_24_spec, metric="dbp", group = "strep_exposure"),
                     legendtitle = "Any exposure", pal =palette_collapsed)
dbp_strep <- plot_layer(plot_base(df=dom_24, metric="dbp", group = "strep_exposure"),
                       calc_wilcox(dom = dom_24, metric = "dbp", group="strep_exposure"),
                       legendtitle = "Any exposure", pal = palette_collapsed)
absolute_strep <- plot_layer(plot_base(df=dom_24, metric="absolute", group = "strep_exposure"),
                        calc_wilcox(dom = dom_24, metric = "absolute", group="strep_exposure"),
                        legendtitle = "Any exposure", pal = palette_collapsed)
rel_strep <- plot_layer(plot_base(df=dom_24, metric="relative", group = "strep_exposure"),
                             calc_wilcox(dom = dom_24, metric = "relative", group="strep_exposure"),
                             legendtitle = "Any exposure", pal = palette_collapsed)
#dominant species
#drop the sample where B bifidum is dominant
dom_24_spec <- dom_24 %>% filter(dominant_species!="Bifidobacterium_bifidum")
leg_spec <- plot_leg(plot_base(df=dom_24_spec, metric="dbp", group = "dominant_species"),
                     legendtitle = "Any exposure", pal = wes_palette("BottleRocket2", 2))
dbp_spec <- print(plot_layer(plot_base(df=dom_24_spec, metric="dbp", group = "dominant_species"),
                        calc_wilcox(dom = dom_24_spec, metric = "dbp", group="dominant_species"),
                        legendtitle = "Any exposure", pal = wes_palette("BottleRocket2", 2)))
absolute_spec <- print(plot_layer(plot_base(df=dom_24_spec, metric="absolute", group = "dominant_species"),
                             calc_wilcox(dom = dom_24_spec, metric = "absolute", group="dominant_species"),
                             legendtitle = "Any exposure", pal = wes_palette("BottleRocket2", 2)))
rel_spec <- print(plot_layer(plot_base(df=dom_24_spec, metric="relative", group = "dominant_species"),
                        calc_wilcox(dom = dom_24_spec, metric = "relative", group="dominant_species"),
                        legendtitle = "Any exposure", pal = wes_palette("BottleRocket2", 2)))
final_plot <- plot_grid(plotlist = dbp_spec, absolute_spec, rel_spec)

plotlistw <- c("dbp_strep", "absolute_strep", "rel_strep","dbp_spec", "absolute_spec", "rel_spec")
final_plot <- plot_grid(plotlist = plotlistw)


final_plot <- plot_grid(plotlist = dbp_strep, absolute_strep, rel_strep,dbp_spec, absolute_spec, rel_spec)
final_plot
ggsave("figures\\dom_violins_between_exposure.pdf", final_plot, height = 29, width = 21, units = "cm", dpi = 300)


#BottleRocket1, BottleRocket2, Rushmore1, Royal1, Royal2, Zissou1, Darjeeling1, Darjeeling2, Chevalier1 , 
#FantasticFox1 , Moonrise1, Moonrise2, Moonrise3, Cavalcanti1, GrandBudapest1, GrandBudapest2, IsleofDogs1, IsleofDogs2


#dominant species approaches significance level - follow up
for (i in (1:length(indices))){
  c2 <- dunnTest(indices[i] ~ strep_exposure, data = dom_big, method="bh")
  dunn_result2 <- c2$res %>%
    mutate(method="Dunn's Test") %>%
    mutate(correction="BH") %>%
    rowwise() %>%
    mutate(sig_code=sigcode(as.numeric(P.adj)))
  results_dunn <- bind_rows(results_dunn, dunn_result2)
}





#write.csv(results_uncollapsed, file="reports\\dominance_uncollapsed_KW.csv")
#write.csv(dunn_result, file="reports\\dominance_uncollapsed_dunn.csv")
write.csv(dom_big, file="reports\\dominance_all_values.csv")

########################################################
### #collapse because need n=5 for Kruskal Wallis    ###
########################################################
# dom_big <- dom_big %>% 
#   mutate(treatment = str_replace(treatment, "FS", "S")) %>%
#    mutate(treatment = str_replace(treatment, "F", "B"))
# 
# #stat testing
# subbed <- dom_big[,1:7]
# indices <- names(subbed)
# rm(subbed)
#summary stats
# dom_big %>% group_by(exp_group) %>%
#   summarise(n=n(), mean = mean(dbp), sd = sd(dbp))
# #timepoint
# dom_big %>% group_by(timepoint) %>%
#   summarise(n=n(), mean = mean(dbp), sd = sd(dbp))

#######################
#generate output file #
#######################
list_of_datasets <- list("Dominance Indices" = dom_big,
                         "Kruskal Wallis, all groups" = results_uncollapsed,
                         "Dunn Testing, Timepoint" = dunn_result,
                         "Kruskal Wallis, collapsed" = results_collapsed,
                         "Dunn Testing, dominant species" = dunn_result2)
openxlsx::write.xlsx(list_of_datasets, file = "reports/dominance.xlsx")

#write.csv(results_collapsed, file="reports\\dominance_collapsed_KW.csv")
# #boxplot(p.value~parameter, data=results)
# #group by timepoint as it has evident influence
# grouping_var <- c("treatment", "dominant_species")
# results <- data.frame()
# for (i in (1:length(indices))){
#   for (j in (1:length(grouping_var))){
#     result <- dom_big %>%
#       group_by(timepoint) %>%
#       do(tidy(kruskal.test(x = .[[indices[i]]], g = .[[grouping_var[j]]])))
#     result <- result %>% 
#       mutate(metric=indices[i]) %>%
#       mutate(group=grouping_var[j])
#     results <- bind_rows(results, result)
#   }
# }


###################################################
results <- data.frame()
for (j in (1:length(grouping_var))){
  result <- dom_big %>% kruskal_test(as.numeric(dbp) ~ grouping_var[j])
  #result <- result %>% 
    #mutate(metric=indices[i]) %>%
    #mutate(group=group_var)
  results <- bind_rows(results, result)
}
dom_big$
grouping_var[2]
length(indices)
grouping_var <- c("treatment", "exp_group", "dominant_species")
results <- data.frame()
for (i in (1:length(indices))){
  for (j in (1:length(grouping_var))){
    result <- dom_big %>%
      group_by(timepoint) %>%
      do(tidy(kruskal.test(x = .[[indices[i]]], g = .[[grouping_var[j]]])))
    result <- result %>% 
      mutate(metric=indices[i]) %>%
      mutate(group=grouping_var[j])
    results <- bind_rows(results, result)
  }
}
dom_big %>%
  summarise(fit = list(kruskal.test(simpson ~ timepoint) %>% tidy)) %>% 
  unnest_wider(fit)

boxplot(dbp~exp_group, data=dom_big)
write.csv(dom_big, file="reports\\dominance.csv")
write.csv(results, file="reports\\dominance_stats.csv")

###############
# Plots
###############

dom <- dom_big %>% 
  dplyr::select(c(dbp, sample_name,exp_group, dominant_species, timepoint, treatment, replicate))

#multiple up dominance indices
dom2 <- dom %>%
  group_by(exp_group, dominant_species, timepoint, replicate, treatment) %>%
  mutate(dom_sum = sum(dbp)) %>% ungroup()

dom2$timepoint = factor(dom$timepoint, levels=c("T_0", "T_6", "T_24"))
ggplot(dom2, aes(x=treatment, y=dominant_species, color=treatment))+
  geom_point(data = dom2,  name="Berger–Parker index",
             aes(x = treatment, y = dominant_species,
                 size=dom_sum*10)) +
  xlab("") + ylab("") +
  scale_color_manual(name="Treatment",
                     labels=c('Control', 
                              expression(paste(italic("F. nucleatum"))), 
                              expression(paste(italic("F. nucleatum/S. salivarius"))), 
                              expression(paste(italic('S. salivarius')))),
                     values = wes_palette("Darjeeling1", n = 4, type = "discrete"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face="italic"))+
  scale_size_continuous(name="Berger–Parker index (summed)")+
  # range = c(0.0,7.0),
  #breaks = c(0.2, 2, 6),
  # labels = c("0.3", "0.4", "1.0"))+
  facet_wrap(~timepoint+replicate)
ggsave("figures\\dominance.emf", width = 10, height = 10, dpi = 300, device = {function(filename, ...) devEMF::emf(file = filename, ...)})





#plot alpha diversity against dominance indices
alpha <-microbiome::alpha(phy2, index = "all")
#join table of dominance with table of alpha diversities by sample
rownames(alpha)->alpha$sample_name
df <- dplyr::left_join(dom_big, alpha, by="sample_name")
df$timepoint <- factor(df$timepoint, levels = c("T_0", "T_6", "T_24"))
#plot alpha against dominance?
ggplot(df, aes(x=treatment, y=diversity_inverse_simpson), color=timepoint) +
  geom_point(aes(color=factor(timepoint)))
#theme(axis.text.x = element_blank())
ggbetweenstats(data=df, 
               x=treatment, 
               y= diversity_inverse_simpson,
               type = "nonparametric")

p0<-grouped_ggbetweenstats(data=df, 
                           x=treatment, 
                           y= diversity_inverse_simpson,
                           type = "nonparametric",
                           grouping.var=timepoint)
p0<-grouped_ggbetweenstats(data=df, 
                           x=timepoint, 
                           y= diversity_inverse_simpson,
                           type = "nonparametric",
                           grouping.var=treatment)
p0
df0 <- extract_stats(p0)
df0p <- df0$pairwise_comparisons_data
df0p <- df0$subtitle_data

#scale_color_viridis_d("Sex", breaks = c(-1,1), labels = c("female", "male")) 


###stat test

model <- lm(gini ~ treatment, data = dom_big)
summary(model)$coef

##################################################################
#Heatmap and line plots

#pull out relative abundances of dominant taxa ()
dominant_taxa <- unique(dom_big$dominant_species)
mpa_table_t <- t(mpa_table)
mpa_table_t[,colnames(mpa_table_t) %in% dominant_taxa] -> mpa2

cbind(mpa2, info) -> df
df$timepoint=factor(df$timepoint, levels=c("T_0", "T_6", "T_24"))
df <- df %>% dplyr::select(Bacteroides_uniformis, Phocaeicola_vulgatus, Streptococcus_salivarius, Bifidobacterium_adolescentis, treatment, timepoint)
df$sample <- rownames(df)
df2 <- pivot_longer(df, 
                    cols=c(Bacteroides_uniformis, Phocaeicola_vulgatus, Streptococcus_salivarius, Bifidobacterium_adolescentis),
                    names_to = c("species"),
                    values_to = c("abundance"))

#levels(df$treatment)= c("B"= "Control", "F" = expression(paste(italic("F. nucleatum"))), "FS" = expression(paste(italic("F.nucleatum/S. salivarius"))), "S" = expression(paste(italic("S. salivarius"))))
treatment_names <- c('B' = "Control", 'F' = "F. nucleatum",'FS' = "F.nucleatum/S. salivarius",'S' = "S. salivarius")
#scatter/lineplot
ggplot(data=df2, aes(x=timepoint, y=abundance, color=species, group=species, alpha=5/10)) +
  geom_point() + geom_smooth(se=FALSE) + facet_wrap(~treatment,labeller=as_labeller(treatment_names)) +
  scale_color_viridis(name="Species", discrete=TRUE) +
  guides(alpha="none") + xlab("Relative abundance") + ylab("Timepoint") +
  theme(strip.text = element_text(size=12, face = "italic"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12, face= "italic"),
        legend.position = "bottom",
        legend.box = "horizontal")
ggsave("figures\\dominance_geom_line.emf", width = 10, height = 10, dpi = 600, device = {function(filename, ...) devEMF::emf(file = filename, ...)})

####repeat line plot including fuso
#pull out relative abundances of dominant taxa ()
#dominant_taxa_f <- c(dominant_taxa, "Fusobacterium_nucleatum")
#mpa_table_t <- t(mpa_table)
#mpa_table_t[,colnames(mpa_table_t) %in% dominant_taxa_f] -> mpa2_f

#cbind(mpa2_f, info) -> df
#df$timepoint=factor(df$timepoint, levels=c("T_0", "T_6", "T_24"))
#df <- df %>% dplyr::select(Bacteroides_uniformis, Phocaeicola_vulgatus, Streptococcus_salivarius, Bifidobacterium_adolescentis, Fusobacterium_nucleatum, treatment, timepoint)
#df$sample <- rownames(df)
#df2 <- pivot_longer(df, 
#                    cols=c(Bacteroides_uniformis, Phocaeicola_vulgatus, Streptococcus_salivarius, Bifidobacterium_adolescentis,Fusobacterium_nucleatum),
#                    names_to = c("species"),
#                   values_to = c("abundance"))

#scatter/lineplot
#ggplot(data=df2, aes(x=timepoint, y=abundance, color=species, group=species, alpha=5/10)) +
#  geom_point() + geom_smooth(se=FALSE) + facet_wrap(~treatment) +
#  scale_color_viridis(discrete=TRUE) +
#  guides(alpha="none")

#HEATMAP
mat <- as.matrix(t(mpa2))
colours <- list('treatment' = c('B' = '#9A8822', 'F' = '#F8AFA8', 'FS'='#FDDDA0', 'S'='#74A089'),
                'timepoint' = c('T_0' = '#899DA4', 'T_6' = '#C93312', 'T_24' = '#DC863B'))

dom_map <- ComplexHeatmap::Heatmap(mat,
                                   name="Relative abundance",
                                   col=viridis(3),
                                   show_column_names = FALSE,
                                   row_names_gp = gpar(fontsize = 17, fontface = "italic"),
                                   width = unit (150, "mm"), 
                                   bottom_annotation = ComplexHeatmap::HeatmapAnnotation(treatment=info[colnames(mat),]$treatment,
                                                                                         timepoint=info[colnames(mat),]$timepoint,
                                                                                         col=colours, 
                                                                                         height = unit(5, "mm"),
                                                                                         gp = gpar(fontsize = 12, fontface = "bold")
                                   ),
                                   heatmap_legend_param = list(title_gp = gpar(fontsize = 12, fontface = "bold"),
                                                               labels_gp = gpar(fontsize = 12),  
                                                               grid_height=unit(5, "mm"),
                                                               grid_width=unit(20,"mm"),
                                                               legend_direction = "horizontal"))
png(file="figures\\dominant_heatmap.png",width=17,height=12,units="in",res=1200)
draw(dom_map, heatmap_legend_side ="bottom", annotation_legend_side = "bottom")
dev.off()