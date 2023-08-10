#superfocus t-24 plots
#SUPERFOCUS DIVERSITY lvl 1
library(ggplot2)
library(dplyr)
library(vegan)
theme_set(theme_classic())
palette_groups <- c("#332211","#ffcc44", "#44aacc", "#bb2211")
palette_time <- c("#9A8822", "#F8AFA8", "#74A089")
palette_collapsed <- c("#FF0000", "#00A08A")
###################
# Read and clean  #
###################

sf <- read.csv("data/superfocus_1.csv",skip=4, row.names = 1, check.names = FALSE)
#read sample metadata
meta <- read.csv("cleaned_data/metadata.csv", nrows=36, row.names=1)
meta <- meta %>% filter(timepoint=="T_24") %>%
  mutate(strep_exposure = case_when(
    treatment == "FS" ~ "Y", 
    treatment == "S" ~ "Y",
    treatment == "F" ~ "N",
    treatment == "B" ~ "N"))
samplenames <- rownames(meta)
sf2 <- sf %>%
  dplyr::select(contains("%")) %>% #filter to relative abundance entries
  dplyr::select(matches(paste(samplenames, collapse = "|"))) %>% #select nisin samples
  dplyr::select(seq(2, ncol(.), by = 2)) %>% #select every other (forward/reverse reads)
  rename_with(~ gsub("_microbial_2.fastq %", "", .))
sf2 <- sf2[rowSums(sf2) >= 0.01,]
sf2 <- t(sf2)
rm(sf)

#PLOT1 SF L1
dist_mat_sf_filt = vegdist(sf2, method="robust.aitchison")
cmd_res = cmdscale(dist_mat_sf_filt, 
                   k = (nrow(sf2) - 1),
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
#GET LEGEND
p1 = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=strep_exposure)) + 
  geom_mark_ellipse(aes(fill=strep_exposure),
                    alpha=0.05,
                    linetype = 2) +
  scale_color_manual(name="Any exposure",
                     labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                              expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                     values = palette_collapsed) + 
  scale_fill_manual(name="Any exposure",
                    labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                             expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                    values = palette_collapsed) +
  labs(shape="Treatment") +
  theme(legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size=8),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.position="bottom", legend.box = "horizontal",
        legend.background = element_rect(colour="black",linewidth = 0.5, linetype="solid"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  ylab(paste("PC2 (",PC2,"%)",sep=""))+
  xlab(paste("PC1 (",PC1,"%)",sep=""))+
  geom_point(aes(shape=treatment),size=1, alpha=1) +
  scale_shape_manual(name="Treatment",
                     labels=c("Control", #b
                              expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                              expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                              expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                     values =  c(15,16,17,18)) +
  guides(fill= "none",
         color = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE),
         shape = guide_legend(title.position = "top",
                              nrow=4, byrow=TRUE,
                              override.aes=list(col=palette_collapsed))) +  
  ggtitle("SUPER-FOCUS Level 1 Pathways")
legend1 <- get_legend(p1)
p1 <- p1 + expand_limits(x = c(-.7, .7),y = c(-0.4, 0.5)) + theme(legend.position = "none")
p1
rm(cmd_res, pcoa_df, pcoa_meta, meta, sf2, dist_mat_sf_filt, PC1, PC2)

#plot 2 level 2
sf <- read.csv("data/superfocus_2.csv",skip=4, row.names = 1, check.names = FALSE)
#read sample metadata
meta <- read.csv("cleaned_data/metadata.csv", nrows=36, row.names=1)
meta <- meta %>% filter(timepoint=="T_24") %>%
  mutate(strep_exposure = case_when(
    treatment == "FS" ~ "Y", 
    treatment == "S" ~ "Y",
    treatment == "F" ~ "N",
    treatment == "B" ~ "N"))
samplenames <- rownames(meta)
sf2 <- sf %>%
  dplyr::select(contains("%")) %>% #filter to relative abundance entries
  dplyr::select(matches(paste(samplenames, collapse = "|"))) %>% #select nisin samples
  dplyr::select(seq(2, ncol(.), by = 2)) %>% #select every other (forward/reverse reads)
  rename_with(~ gsub("_microbial_2.fastq %", "", .))
sf2 <- sf2[rowSums(sf2) >= 0.01,]
sf2 <- t(sf2)
rm(sf)
#plot
dist_mat_sf_filt = vegdist(sf2, method="robust.aitchison")
cmd_res = cmdscale(dist_mat_sf_filt, 
                   k = (nrow(sf2) - 1),
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
p2 = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=strep_exposure)) + 
  geom_mark_ellipse(aes(fill=strep_exposure),
                    alpha=0.05,
                    linetype = 2) +
  scale_color_manual(name="Any exposure",
                     labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                              expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                     values = palette_collapsed) + 
  scale_fill_manual(name="Any exposure",
                    labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                             expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                    values = palette_collapsed) +
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
  geom_point(aes(shape=treatment),size=1, alpha=1) +
  scale_shape_manual(name="Treatment",
                     labels=c("Control", #b
                              expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                              expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                              expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                     values =  c(15,16,17,18)) +
  guides(fill= "none",
         color = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE),
         shape = guide_legend(title.position = "top",
                              nrow=4, byrow=TRUE,
                              override.aes=list(col=palette_collapsed))) +  
  ggtitle("SUPER-FOCUS Level 2 Pathways")
p2 <- p2 + expand_limits(x = c(-4, 4),y = c(-5, 5)) + theme(legend.position = "none")
p2
rm(cmd_res, pcoa_df, pcoa_meta, meta, sf2, dist_mat_sf_filt, PC1, PC2)

#plot 3 level 3
sf <- read.csv("data/superfocus_3.csv",skip=4, row.names = 1, check.names = FALSE)
#read sample metadata
meta <- read.csv("cleaned_data/metadata.csv", nrows=36, row.names=1)
meta <- meta %>% filter(timepoint=="T_24") %>%
  mutate(strep_exposure = case_when(
    treatment == "FS" ~ "Y", 
    treatment == "S" ~ "Y",
    treatment == "F" ~ "N",
    treatment == "B" ~ "N"))
samplenames <- rownames(meta)
sf2 <- sf %>%
  dplyr::select(contains("%")) %>% #filter to relative abundance entries
  dplyr::select(matches(paste(samplenames, collapse = "|"))) %>% #select nisin samples
  dplyr::select(seq(2, ncol(.), by = 2)) %>% #select every other (forward/reverse reads)
  rename_with(~ gsub("_microbial_2.fastq %", "", .))
sf2 <- sf2[rowSums(sf2) >= 0.01,]
sf2 <- t(sf2)
rm(sf)
#plot
dist_mat_sf_filt = vegdist(sf2, method="robust.aitchison")
cmd_res = cmdscale(dist_mat_sf_filt, 
                   k = (nrow(sf2) - 1),
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
p3 = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=strep_exposure)) + 
  geom_mark_ellipse(aes(fill=strep_exposure),
                    alpha=0.05,
                    linetype = 2) +
  scale_color_manual(name="Any exposure",
                     labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                              expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                     values = palette_collapsed) + 
  scale_fill_manual(name="Any exposure",
                    labels=c(expression(paste("No ", italic("S. salivarius"), " DPC6487")),
                             expression(paste("Any ", italic("S. salivarius"), " DPC6487"))),
                    values = palette_collapsed) +
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
  geom_point(aes(shape=treatment),size=1, alpha=1) +
  scale_shape_manual(name="Treatment",
                     labels=c("Control", #b
                              expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                              expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                              expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                     values =  c(15,16,17,18)) +
  guides(fill= "none",
         color = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE),
         shape = guide_legend(title.position = "top",
                              nrow=4, byrow=TRUE,
                              override.aes=list(col=palette_collapsed))) +  
  ggtitle("SUPER-FOCUS Level 3 Pathways")
p3 <- p3 + expand_limits(x = c(-10, 6),y = c(-10, 10)) + theme(legend.position = "none")
p3
rm(cmd_res, pcoa_df, pcoa_meta, meta, sf2, dist_mat_sf_filt, PC1, PC2)

#plot together
plots <- plot_grid(p1, p2, p3, legend1, nrow=2, labels="AUTO")

title <- ggdraw() + draw_label("Functional Beta-diversity (T_24)", fontface='bold')
plots2 <- plot_grid(title, plots, ncol=1, rel_heights=c(0.1, 1)) 
plots2
ggsave("figures\\func_sf_all.pdf", plots2, height = 21, width = 24, units = "cm", dpi = 300)
ggsave("figures\\func_sf_all.png", plots2, height = 21, width = 24, units = "cm", dpi = 300)
