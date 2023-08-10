#humann t_24 plots strat and unstrat
#humann unstratified
library(ggplot2)
library(dplyr)
library(vegan)
library(data.table)
library(ggforce)
library(wesanderson)
library(cowplot)
theme_set(theme_classic())
palette_groups <- c("#332211","#ffcc44", "#44aacc", "#bb2211")
palette_time <- c("#9A8822", "#F8AFA8", "#74A089")
palette_collapsed <- c("#FF0000", "#00A08A")

################################
#collapsed and only t24
################################
hm <- as.data.frame(fread("data/pathabundance_cpm_unstratified.tsv"))
hm <- hm %>%
  rename(pathway = `# Pathway`)
rownames(hm) <- hm$pathway
hm <- hm[,-1]
#read sample metadata
meta <- read.csv("cleaned_data/metadata.csv", nrows=36, row.names=1)
meta <- meta %>%
  mutate(strep_exposure = case_when(
    treatment == "FS" ~ "Y", 
    treatment == "S" ~ "Y",
    treatment == "F" ~ "N",
    treatment == "B" ~ "N"))
samplenames <- rownames(meta)
hm2 <- hm %>%
  rename_with(~ gsub("_Abundance", "", .)) %>%
  #select(contains("%")) %>% #filter to relative abundance entries
  dplyr::select(matches(paste(samplenames, collapse = "|"))) #select nisin samples
#dplyr::select(seq(2, ncol(.), by = 2)) %>% #select every other (forward/reverse reads)
#rename_with(~ gsub("_microbial_2.fastq %", "", .))
hm2 <- hm2[rowSums(hm2) >= 0.01,]
hm2 <- t(hm2)
rm(hm)
#need new ordination filtered to T24
meta_filt <- meta %>%
  filter(timepoint=="T_24")
hm2[rownames(hm2) %in% rownames(meta_filt),] -> hm_table_filt
#remove empty rows
#HUMAnN CPM - sum-normalizing RPKs to relative abundance and then multiplying by 1e6 to yield bigger numbers,
#hm_table_filt <- hm_table_filt[,colSums(hm_table_filt) >= 1000]
#mpa_table_filt_mat <- t(mpa_table_filt)
#robust aitchisons transformation
dist_mat_hm_filt = vegdist(hm_table_filt, method="robust.aitchison")
#plot
cmd_res = cmdscale(dist_mat_hm_filt, 
                   k = (nrow(hm_table_filt) - 1),
                   eig = TRUE)
# calculate the proportion of variance in the data which is explained by 
#the first two PCoA axes
PC1 <- round(100*(cmd_res$eig[1]/(sum(cmd_res$eig))), digits = 2)
PC2 <- round(100*(cmd_res$eig[2]/(sum(cmd_res$eig))), digits = 2)
#extract details of features
pcoa_df = tibble(PC1 = cmd_res$points[,1], 
                 PC2 = cmd_res$points[,2])
#combine pcoa with sample data
pcoa_meta = bind_cols(pcoa_df, meta_filt)
p4 = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=strep_exposure)) + 
  geom_mark_ellipse(aes(fill=strep_exposure),
                    alpha=0.05,
                    linetype = 2) +
  scale_color_manual(name="Inoculation",
                     labels=c(expression(paste("No ",italic("S. salivarius"), " exposure")), 
                              expression(paste("Any ",italic("S. salivarius"), " exposure"))), 
                     values = wes_palette("Darjeeling1", n = 2, type = "discrete")) +
  scale_fill_manual(name="Inoculation",
                    labels=c(expression(paste("No ",italic("S. salivarius"), " exposure")), 
                             expression(paste("Any ",italic("S. salivarius"), " exposure"))), 
                    values = wes_palette("Darjeeling1", n = 2, type = "discrete")) +
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
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=wes_palette("Darjeeling1", n = 2, type = "discrete"))),
         shape = guide_legend(title.position = "top",
                              nrow=2, byrow=TRUE,
                              override.aes=list(col=c("#FF0000", "#FF0000","#00A08A","#00A08A")))) +
  ggtitle(expression(paste("T_24 only ~ ",italic("S. salivarius"), " exposure")))
p4 <- p4 + expand_limits(x = c(-20,20),y = c(-12, 8)) + theme(legend.position = "none")
p4

###################
#'stratified'
#'
#################
rm(cmd_res, filt_df_all, hm_table_filt, hm2, md_res, meta, meta_filt, pcoa_df, pcoa_meta, dist_mat_hm_filt, dist_mat_ra)

###################
# Read and clean  #
###################
hm <- as.data.frame(fread("data/pathabundance_cpm_stratified.tsv"))
hm <- hm %>%
  rename(pathway = `# Pathway`)
rownames(hm) <- hm$pathway
hm <- hm[,-1]
#read sample metadata
meta <- read.csv("cleaned_data/metadata.csv", nrows=36, row.names=1)
meta <- meta %>%
  mutate(strep_exposure = case_when(
    treatment == "FS" ~ "Y", 
    treatment == "S" ~ "Y",
    treatment == "F" ~ "N",
    treatment == "B" ~ "N"))
samplenames <- rownames(meta)
hm2 <- hm %>%
  rename_with(~ gsub("_Abundance", "", .)) %>%
  #select(contains("%")) %>% #filter to relative abundance entries
  dplyr::select(matches(paste(samplenames, collapse = "|"))) #select nisin samples
#dplyr::select(seq(2, ncol(.), by = 2)) %>% #select every other (forward/reverse reads)
#rename_with(~ gsub("_microbial_2.fastq %", "", .))
hm2 <- hm2[rowSums(hm2) >= 0.01,]
hm2 <- t(hm2)
rm(hm)

meta_filt <- meta %>%
  filter(timepoint=="T_24")
hm2[rownames(hm2) %in% rownames(meta_filt),] -> hm_table_filt
#remove empty rows
#HUMAnN CPM - sum-normalizing RPKs to relative abundance and then multiplying by 1e6 to yield bigger numbers,
#hm_table_filt <- hm_table_filt[,colSums(hm_table_filt) >= 1000]
#mpa_table_filt_mat <- t(mpa_table_filt)
#robust aitchisons transformation
dist_mat_hm_filt = vegdist(hm_table_filt, method="robust.aitchison")

#plot
cmd_res = cmdscale(dist_mat_hm_filt, 
                   k = (nrow(hm_table_filt) - 1),
                   eig = TRUE)

# calculate the proportion of variance in the data which is explained by 
#the first two PCoA axes
PC1 <- round(100*(cmd_res$eig[1]/(sum(cmd_res$eig))), digits = 2)
PC2 <- round(100*(cmd_res$eig[2]/(sum(cmd_res$eig))), digits = 2)
#extract details of features
pcoa_df = tibble(PC1 = cmd_res$points[,1], 
                 PC2 = cmd_res$points[,2])

#combine pcoa with sample data
pcoa_meta = bind_cols(pcoa_df, meta_filt)
p5 = ggplot(pcoa_meta,aes(x = PC1, y = PC2, color=strep_exposure)) + 
  geom_mark_ellipse(aes(fill=strep_exposure),
                    alpha=0.05,
                    linetype = 2) +
  scale_color_manual(name="Inoculation",
                     labels=c(expression(paste("No ",italic("S. salivarius"), " exposure")), 
                              expression(paste("Any ",italic("S. salivarius"), " exposure"))), 
                     values = wes_palette("Darjeeling1", n = 2, type = "discrete")) +
  scale_fill_manual(name="Inoculation",
                    labels=c(expression(paste("No ",italic("S. salivarius"), " exposure")), 
                             expression(paste("Any ",italic("S. salivarius"), " exposure"))), 
                    values = wes_palette("Darjeeling1", n = 2, type = "discrete")) +
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
         col= guide_legend(title.position = "top",
                           nrow=2, byrow=TRUE,
                           override.aes=list(fill=wes_palette("Darjeeling1", n = 2, type = "discrete"))),
         shape = guide_legend(title.position = "top",
                              nrow=4, byrow=TRUE,
                              override.aes=list(col=c("#FF0000", "#FF0000","#00A08A","#00A08A")))) +
  ggtitle(expression(paste("T_24 only ~ ",italic("S. salivarius"), " exposure")))
leg1 <- get_legend(p5)
p5 <- p5 + expand_limits(x = c(-45, 45),y = c(-35, 30)) + theme(legend.position = "none")
p5

########################
#adjust titles
p4 <- p4 + ggtitle("HUMAnN4 (Unstratified)")
p5 <- p5 + ggtitle("HUMAnN4 (Stratified)")
p0 <- plot_grid(p4,p5, labels=c("D","E"))
p0
p1<- plot_grid(leg1)
ggsave("figures\\humann4_strat_unstrat_t24.png", p0, dpi=600)
ggsave("figures\\humann4_strat_unstrat_t24_leg.png", p1, dpi=600)
