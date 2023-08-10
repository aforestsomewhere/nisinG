#' Script to calculate alpha diversity on NisinG samples
#' Developed by A. Kate Falà (https://github.com/aforestsomewhere)
#' all timepoints and treatments

####################
### Packages     ###
####################
library(tidyverse)
library(dplyr)
library(vegan)
library(cowplot)
library(devEMF)
library(ggforce)
library(webshot)
library(ggplot2)
library(rstatix)
library(ggsignif)
theme_set(theme_classic())

####################
### read data    ###
####################
#read MPA4 data
mpa_table <- read.table("cleaned_data/mpa_table.tsv", comment.char = '#', sep = '\t', header = TRUE, row.names=NULL)
#filter to species level assignments
mpa_table <- mpa_table[grep('s__',mpa_table[,1]),]
mpa_table[,1] <- gsub(".+\\|s__", "", mpa_table[,1])
#remove duplicates by filtering to t__ level
mpa_table <- mpa_table %>% filter(grepl('t__',row.names))
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]
#generate matrix of abundance data
mpa_mat = t(as.matrix(mpa_table))
mpa_mat <- mpa_mat[sort(rownames(mpa_mat)),] #order samples (rows) to properly merge with infO
mpa_mat_t <- t(mpa_mat)
#read sample metadata
meta <- read.csv("cleaned_data/metadata.csv", nrows=36, row.names=1)

####################################
### calculate alpha diversity    ###
####################################
richness_vec <-(specnumber(mpa_mat))
shannon_vec <- vegan::diversity(mpa_mat, index = "shannon")
simp_vec <- vegan::diversity(mpa_mat, index = "simpson")
inv_simp_vec <- vegan::diversity(mpa_mat, index = "invsimpson")
alpha_df <- cbind.data.frame(comb_id = rownames(meta), richness = richness_vec, shannon = shannon_vec, simpson = simp_vec, inv_simpson = inv_simp_vec, meta = meta)
#alpha_df$meta.treatment = factor(alpha_df$meta.treatment)
#Rationalise column names
names(alpha_df) <- c("sample_name", "richness","shannon","simpson", "inv_simpson", "sample_id", "treatment", "replicate", "timepoint", "treatment_replicate", "exp_group")
alpha_df <- alpha_df %>% 
  select(sample_name, sample_id, treatment, replicate, treatment_replicate, timepoint, exp_group, richness, shannon, simpson, inv_simpson) %>%
  mutate(timepoint=factor(timepoint, levels=c("T_0", "T_6", "T_24")))
write.csv(alpha_df, "reports/alpha_all.csv", row.names = FALSE)

################################################
### visualise differences in alpha diversity ###
### due to sample size n=3, no testing for   ###
### individual groups / timepoints           ###
### violin plots to show points              ###
################################################
theme_set(theme_classic())
#first plto to extract common legend
pal1 <- c("#332211", "#ffcc44", "#44aacc", "#bb2211")
legend <-ggplot(data = alpha_df, aes(x=treatment, y=richness, fill=treatment))+
  geom_violin(trim=FALSE) +geom_point(color="black", alpha=0.8)  +
  scale_fill_manual(name="Treatment group",
                    labels=c("Control", #b
                             expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                             expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                             expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                    values = pal1) + 
  xlab("") + ylab("Observed species") +
  theme(legend.direction = "horizontal",
        legend.text = element_text(size=10, hjust=0),
        legend.position = "bottom",
        legend.title = element_text(face="bold", size=14, hjust=0),
        axis.ticks.x = element_blank(),
         axis.text.x = element_blank()) + facet_wrap(~timepoint) +
  guides(fill = guide_legend(title.position = "top"))
  #                              nrow=2, byrow=TRUE))
legend <- get_legend(legend)
#Richness
p1 <-ggplot(data = alpha_df, aes(x=treatment, y=richness, fill=treatment))+
  geom_violin(trim=FALSE) +geom_point(color="black", alpha=0.8)  +
  scale_fill_manual(name="Treatment Group",
                    labels=c("Control", #b
                             expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                             expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                             expression(paste(italic("S. salivarius"), " DPC6487"))), #s
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
  geom_violin(trim=FALSE) +geom_point(color="black", alpha=0.8)  +
  scale_fill_manual(name="Treatment Group",
                    labels=c("Control", #b
                             expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                             expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                             expression(paste(italic("S. salivarius"), " DPC6487"))), #s
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
# #Simpson
# p3 <-ggplot(data = alpha_df, aes(x=treatment, y=simpson, fill=treatment))+
#   geom_violin() +geom_point(color="black", alpha=0.8)  +
#   scale_fill_manual(name="Treatment Group",
#                     labels=c("Control", #b
#                              expression(paste(italic("F. nucleatum"), " DSM15643")), #f
#                              expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
#                              expression(paste(italic("S. salivarius"), " DPC6487"))), #s
#                     values = pal1) + 
#   xlab("") + ylab("Simpson's Index") +
#   theme(legend.box = "vertical", 
#         legend.position = "right",
#         legend.title = element_text(face="bold"),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.title.y = element_text(size=12, face="bold", colour = "black"),
#         strip.text.x = element_text(size = 12, face="bold", colour = "black" ),
#         strip.text.y = element_text(size = 12, face="bold", colour = "black")) + facet_wrap(~timepoint) +
#   guides(fill="none", alpha="none", color="none")
# p3
#Inverse Simpson
p4 <-ggplot(data = alpha_df, aes(x=treatment, y=inv_simpson, fill=treatment))+
  geom_violin(trim=FALSE) +geom_point(color="black", alpha=0.8)  +
  scale_fill_manual(name="Treatment Group",
                    labels=c("Control", #b
                             expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                             expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                             expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                    values = pal1) + 
  xlab("") + ylab("Inverse Simpson's Index") +
  theme(legend.box = "horizontal", 
        legend.position = "right",
        legend.title = element_text(face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        strip.text.x = element_text(size = 12, face="bold", colour = "black" ),
        strip.text.y = element_text(size = 12, face="bold", colour = "black")) + facet_wrap(~timepoint) +
  guides(fill="none", alpha="none", color="none")
p4

all_alpha <- plot_grid(p1,p2,p4, nrow=3)
all_alpha
all_alpha_leg <- plot_grid(all_alpha, legend, 
                nrow=2, align='v', rel_heights=c(1, 0.2),
                rel_widths=c(1, 2)) + theme(plot.background = element_rect(fill="#FFFFFF", color = NA))
all_alpha_leg
ggsave("figures/alpha_all.pdf", width=21, height=29, dpi=300, units="cm")
ggsave("figures/alpha_all.emf", width = 10, height = 10, dpi = 600, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
write.csv(alpha_df, "reports\\alpha_all_metrics.csv", row.names = FALSE)
