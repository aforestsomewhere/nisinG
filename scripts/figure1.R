#beta - aitchison
library(ComplexHeatmap)
library(dplyr)
library(vegan)
library(tidyr)
library(ggplot2)
library(wesanderson)
#update.packages("wesanderson")
library(ggforce)
library(cowplot)
library(viridis)
#library(grid)

theme_set(theme_classic())
palette_groups <- c("#332211","#ffcc44", "#44aacc", "#bb2211")
palette_time <- c("#9A8822", "#F8AFA8", "#74A089")

meta <- read.csv("cleaned_data/metadata.csv", nrows=36, row.names=1) %>%
  janitor::clean_names() %>%
  mutate(timepoint=factor(timepoint, levels=c("T_0", "T_6", "T_24"))) %>%
  mutate(strep_exposure = case_when(
    treatment == "FS" ~ "Y", 
    treatment == "S" ~ "Y",
    treatment == "F" ~ "N",
    treatment == "B" ~ "N"))
mpa_table <- read.table("cleaned_data/mpa_table.tsv", comment.char = '#', sep = '\t', header = TRUE, row.names=NULL)
#filter to species level assignments
mpa_table <- mpa_table[grep('s__',mpa_table[,1]),]
mpa_table[,1] <- gsub(".+\\|s__", "", mpa_table[,1])
mpa_table <- mpa_table %>% filter(grepl('t__',row.names))
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]
#generate matrix of abundance data
mpa_mat = t(as.matrix(mpa_table))


#robust aitchisons transformation
dist_mat_ra = vegdist(mpa_mat, method="robust.aitchison")
cmd_res = cmdscale(dist_mat_ra, 
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
  ggtitle("All timepoints ~ timepoint")
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
  ggtitle("All timepoints ~ treatment")
p3
p3a <- p3 + expand_limits(x = c(-17,9), y = c(-15,10))
p3a
p3leg <- get_legend(p3a)
p3b <- p3a + theme(legend.position='none')
p3c <- ggdraw(plot_grid(p3b, p3leg,
                        nrow=2, rel_heights=c(10, 2),
                        rel_widths=c(17, 2),
                        axis = 'l', align = 'h'))


################################
#collapsed and only t24
################################
#need new ordination filtered to T24
meta_filt <- meta %>%
  filter(timepoint=="T_24")
mpa_table[,colnames(mpa_table) %in% rownames(meta_filt)] -> mpa_table_filt
#remove empty rows
mpa_table_filt <- mpa_table_filt[rowSums(mpa_table_filt) >= 0.000001,]
mpa_table_filt_mat <- t(mpa_table_filt)
#robust aitchisons transformation
dist_mat_ra_filt = vegdist(mpa_table_filt_mat, method="robust.aitchison")
cmd_res = cmdscale(dist_mat_ra_filt, 
                   k = (nrow(mpa_table_filt_mat) - 1),
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
p4
#theme_margin <- theme(legend.box.margin = margin(100, 10, 100, 10))
#legend  <- cowplot::get_legend(p4 + theme_margin)


p4a <- p4 + expand_limits(x = c(-15,15), y = c(-11,10))
p4leg <- get_legend(p4a)
p4b <- p4a + theme(legend.position='none')
p4c <- ggdraw(plot_grid(p4b, p4leg,
                nrow=2, rel_heights=c(10, 2),
                rel_widths=c(17, 2),
                axis = 'l', align = 'h'))
p4c
p0 <- ggdraw(plot_grid(p2c,p3c,p4c,
             align = "h", ncol = 1,
             labels = "AUTO"))
p0
ggsave("figures/figure1.png", width=11, height=29, dpi=300, units="cm")
ggsave("figures/figure1.pdf", width=11, height=29, dpi=300, units="cm")
##############################
# ComplexHeatmap::Heatmap abundance plots    #
##############################
#ComplexHeatmap::Heatmap at T24
#transform so species are vertical
mat <- t(mpa_table_filt_mat)
#subset to highly abundant species
mat <- mat[rowSums(mat) >= 5,]
#remove "t__" sgb assignment for cleaner rownames
#rename t__ level to species level assignments
sp1 <- as.data.frame(row.names(mat))
names(sp1) <- c("id")
sp1 <- sp1 %>%
  separate(id, c("species","sgb"), "__") %>%
  mutate(taxon = gsub("\\|t","", species))
#rejoin
rownames(mat) <- sp1$taxon
#base ComplexHeatmap::Heatmap
ht1 = ComplexHeatmap::Heatmap(mat,
        col=viridis(3),
        name = "T_24")
ht1=draw(ht1) #initialise ComplexHeatmap::Heatmap to extract row order to ensure consistent ordering
row_ord_ht1 <- row_order(ht1)
col_ord_ht1 <- column_order(ht1)
#make annotation
#subset to pertinent metadata
# F_NAME <- expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487"))
# F_NAME <- bquote(italic("F. nuc") ~ DSM15643 + italic("S. sal") ~ DPC6487
# F_NAME <- paste(expression(italic("F. nuc")), " DSM15643 + ", expression(italic("S. sal")), " DPC6487")
# class(F_NAME)
#expression(paste("No ",italic("S. salivarius"), " exposure"
my_meta <- meta_filt %>% dplyr::select(treatment) 
# %>%
#   mutate(treatment = case_when(
#     treatment == "FS" ~ "Y", 
#     treatment == "S" ~ "Y",
#     treatment == "F" ~ "N",
#     treatment == "B" ~ "N"))
htax <- ComplexHeatmap::HeatmapAnnotation(df=my_meta,
                          which="col",
                          simple_anno_size = unit(1, "cm"),
                          annotation_name_gp= gpar(fontsize = 20),
                          show_annotation_name = FALSE,
                          col=list(#strep_exposure=c("Y"="#00A08A", "N"="#FF0000"),
                                   treatment = c("B"= "#332211","F"="#ffcc44", "FS"="#44aacc", "S"="#bb2211")),
                          annotation_legend_param = list(
                            treatment = list(
                              title = "Treatment groups",
                              at = c("B","F","FS", "S"), #CUSTOMISING DISPLAYED TREATMENT GROUP NAMES
                              labels = c("Control", #b
                                         expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                                         expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                                         expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                              title_gp = gpar(fontsize = 20, fontface = "bold"), 
                              labels_gp = gpar(fontsize = 16), 
                              show_annotation_name = FALSE,
                              grid_width=unit(1.5,"cm"))))
                          
ht1 = ComplexHeatmap::Heatmap(mat, #need to recreate without the draw()
              col=viridis(3),
              name = "Relative abundance at T_24 (%)",
              #rect_gp = gpar(col = "white", lwd = 1),
              #column_names_side = "top",
              top_annotation = htax,
              column_order = col_ord_ht1,
              cluster_columns = TRUE,
              show_column_names = FALSE,
              show_column_dend = TRUE, 
              show_row_dend = FALSE,
              row_dend_reorder = FALSE,
              width = unit(15, "cm"), height = unit(30, "cm"),
              heatmap_legend_param = list(color_bar = "continuous",
                                          grid_height=unit(1,"cm"),
                                          #legend_height = unit(5,"cm"),
                                         legend_width = unit(10, "cm"), 
                                          legend_direction = "horizontal",
                                          title_position="topcenter",
                                          labels_gp = gpar(fontsize=20),
                                          legend_label_gp = gpar(fontsize = 22, fontface="bold")),
                                         # labels_gp = gpar(fontsize = 20)),
              row_names_gp = gpar(fontsize = 16, fontface="italic"))
ht1= draw(ht1,heatmap_legend_side="bottom", 
          #legend_title_gp = gpar(fontsize = 20, fontface = "bold"),
          annotation_legend_side = "top") # modify text size for ComplexHeatmap::Heatmap legend(rel abun)
draw(ht1)
#htg = p2 = grid.grabExpr(ht1)
#ggsave("figures/figure1.png", width=11, height=29, dpi=300, units="cm")
pdf("figures/T24_Heatmap.pdf", width=21, height=29)
png(file="figures/Heatmap.png", width=1000, height=2000)
#png(file="figures/ComplexHeatmap::Heatmap.png", width =1240, height=300)
draw(ht1)
dev.off()


#SEPARATE ComplexHeatmap::Heatmap FOR FUSO
fuso <- (t(mpa_table_filt_mat))
#subset to highly abundant species
# Get row names of the matrix
row_names <- rownames(fuso)
# Filter rows with partial string matches
partial_match_rows <- row_names[grepl("Fusobacterium", row_names, ignore.case = TRUE)]
# Subset the matrix based on the filtered row names
fusoy <- fuso[partial_match_rows, , drop = FALSE]
#rename t__ level to species level assignments
sp1 <- as.data.frame(row.names(fusoy))
names(sp1) <- c("id")
sp1 <- sp1 %>%
  separate(id, c("species","sgb"), "__") %>%
  mutate(taxon = gsub("\\|t","", species))
#rejoin
rownames(fusoy) <- sp1$taxon


# ht2 = ComplexHeatmap::Heatmap(fusoy,
#               col=viridis(3),
#               name = "T_24")

#need second annotation object to suppress display of legend
my_meta <- meta_filt %>% select(treatment) 
htax2 <- ComplexHeatmap::HeatmapAnnotation(df=my_meta,
                          which="col",
                          simple_anno_size = unit(1, "cm"),
                          annotation_name_gp= gpar(fontsize = 24),
                          col=list(#strep_exposure=c("Y"="#00A08A", "N"="#FF0000"),
                            treatment = c("B"= "#332211","F"="#ffcc44", "FS"="#44aacc", "S"="#bb2211")),
                          show_legend = FALSE,
                          show_annotation_name = FALSE,
                          annotation_legend_param = list(
                            treatment = list(
                              title = "Treatment groups",
                              at = c("B","F","FS", "S"), #CUSTOMISING DISPLAYED TREATMENT GROUP NAMES
                              labels = c("Control", #b
                                         expression(paste(italic("F. nucleatum"), " DSM15643")), #f
                                         expression(paste(italic("F. nuc"), " DSM15643 + ", italic("S. sal"), " DPC6487")), #fs
                                         expression(paste(italic("S. salivarius"), " DPC6487"))), #s
                              #title_gp = gpar(fontsize = 24, fontface = "bold"), 
                             # labels_gp = gpar(fontsize = 24), 
                              show_annotation_name = FALSE,
                              show_legend = FALSE,
                              grid_width=unit(1.5,"cm"))))
#ht2=draw(ht2) #initialise ComplexHeatmap::Heatmap to extract row order to ensure consistent ordering
#row_ord_ht2 <- row_order(ht2)
#col_ord_ht2 <- column_order(ht2)
ht2 = ComplexHeatmap::Heatmap(fusoy, #need to recreate without the draw()
              col=viridis(3),
              name = "Relative abundance at T_24 (%)",
              #rect_gp = gpar(col = "white", lwd = 1),
              #column_names_side = "top",
              bottom_annotation = htax2,
              #column_order = col_ord_ht1, #
              cluster_columns = TRUE,
              show_column_names = FALSE,
              show_column_dend = TRUE,
              column_dend_side = "bottom",
              show_row_dend = FALSE,
              row_dend_reorder = FALSE,
              width = unit(20, "cm"), height = unit(3, "cm"),
              heatmap_legend_param = list(color_bar = "continuous",
                                          grid_height=unit(1,"cm"),
                                          legend_width = unit(10, "cm"), 
                                          legend_direction = "horizontal",
                                          title_position="topcenter",
                                          #legend_title_gp = gpar(fontsize = 36, fontface = "bold")),
                                          labels_gp = gpar(fontsize=20),
                                         #legend_gp = gpar(fontsize = 34, fontface = "bold")),
                                          legend_label_gp = gpar(fontsize = 22, fontface = "bold")),
              row_names_gp = gpar(fontsize = 24, fontface="italic"))
ht2
#rm(ht2)
#ht2 = draw(ht2)
ht2= draw(ht2,heatmap_legend_side="bottom", legend_title_gp=gpar(fontsize=12))
draw(ht2)
          #nnotation_legend_side = NULL) # modify text size for ComplexHeatmap::Heatmap legend(rel abun)
# library(grid)
# htg = grid.grabExpr(draw(ht1))
# htg2 = grid.grabExpr(ht2,ComplexHeatmap::Heatmap_legend_side="bottom", 
#                      legend_title_gp = gpar(fontsize = 20, fontface = "bold"))
# ht2
#pdf("figures/T24_ComplexHeatmap::Heatmap_fuso.pdf", width=21, height=29)
png(file="figures/fuso.png", width = 4200, height =1500, res=300)
draw(ht2)
dev.off()

#ggsave("figures/T24_ComplexHeatmap::Heatmap.pdf", width=29, height=21, dpi=300, units="cm")


# g = plot_grid(htg, htg2, rel_heights = c(10, 1), 
#           nrow = 2, rel_widths = c(1, 1))
# ggdraw(g)
