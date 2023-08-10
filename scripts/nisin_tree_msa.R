library(dplyr)
library(ggtree)
library(ape)
library(phytools)
library(ggplot2)
library(ggnewscale)
library(viridis)
library(treeio)
library(ggrepel)
library(Biostrings)
library(ggmsa)
library(RColorBrewer)
library(pals)
library(ggnewscale)
library(msa)
library(seqinr)
#options(download.file.method = "wininet")
#if (!requireNamespace("devtools", quietly=TRUE))
#  install.packages("devtools")
#devtools::install_github("YuLab-SMU/ggmsa")


# tr <- treeio::read.tree("data/nisin_clustalw.nexus.treefile")
# ggtree(tr) + geom_tiplab() + geom_nodelab()
# tr
# tr <- read.newick("data/nisin_clustalw.nexus.contree", node.label='support')
# ggtree(tr, right = TRUE) + geom_tiplab() + 
#   geom_text2(aes(subset = !isTip, label=support))+ theme_tree2()
# 
# p <- ggtree(tr) + 
#   geom_label_repel(aes(subset = !isTip & !is.na(label), label=support, fill=support), force=0, color="white") + 
#   scale_fill_viridis_c() + geom_tiplab(aes(label=label, subset= !node %in% c(13)), align=TRUE) + 
#   geom_tiplab(aes(label=label, fontface="bold", subset=node %in% c(13)), align=TRUE) + theme_tree2()+
#   theme(legend.position = "bottom", legend.direction = "horizontal") +
#   guides(fill = guide_colourbar(title="Bootstrap Support", title.position = "top"))
# 
# ######
# msa <- read.nexus.data("data/nisin.fasta")
# msa <- readAAStringSet("data/nisin_clustalw.fasta")
# msa <- readAAStringSet("data/nisin.fasta")
# x <- ape::as.AAbin(x)
# names(msa)
# msa <- Biostrings::AAMultipleAlignment("data/nisin_clustalw.fasta")
# msa <- Biostrings::AAStringSet("data/nisin_clustalw.fasta")
# fasta <- system.file("data/nisin_clustalw.fasta", package="ggtree")
# aastring
# msaplot(ggtree(tr), Biostrings::BStringSet(msa), offset=3, width=2)
# fasta <- system.file("data/nisin_clustalw.fasta", package="ggtree")
# msaplot(p=ggtree(tr), fasta="data/nisin_clustalw.fasta")
# 
# data = tidy_msa(x)
# p + geom_facet(geom=geom_msa, data = data)


d <- read.csv("data/nisin_variants_names.csv")
#d <- d %>%
#  mutate(species = gsub(". ","_", species)) %>%
#  mutate(strain = gsub(" ","", strain)) %>%
#  mutate(strain = gsub("-","", strain))
d$size <- as.character(d$size)
d$strain <- as.character(d$strain)
# p <- ggtree(tr) + geom_tiplab(aes(label=label, subset= !node %in% c(13)), align=TRUE) + 
#   geom_tiplab(aes(label=label, fontface="bold", subset=node %in% c(13)), align=TRUE) + theme_tree2()
# p2 <-ggtree(tr) %<+% d +
#   geom_tiplab(aes(label=paste0('bold(',label,')~italic(', species, ')~', '"', strain, '"')), parse=TRUE)
# p2
#alt tree alignment

#BiocManager::install("msa")
mySequences <- readAAStringSet("data/nisin.fasta")
my_alignment <- msa(mySequences)
#my_alignment_sequence <- msaConvert(my_alignment, type="seqinr::alignment")
#distance_alignment <- dist.alignment(my_alignment_sequence)
## compute phylogenetic tree using neighbor joining
#tree <- bionj(distance_alignment)
f <- function(x) njs(dist.aa(x))
s <- as.AAbin(my_alignment)
tree2 <- f(s)
#ggtree(tree2)
#obtain 100 boostrapped trees
# bstrees <- boot.phylo(tree2, s, f, trees = TRUE)
# ## get proportions of each clade:
# clad <- prop.clades(tree2, bstrees, rooted = TRUE)
# ## get proportions of each bipartition:
# boot <- prop.clades(tree2, bstrees)
# tree2 <- ape::as.phylo(tree2, use.labels=TRUE)
# tree2$Nnode
# tree3 <- makeNodeLabel(tree2, "number", prefix="")


#a#a#a
myBoots <- boot.phylo(tree2, s, f, rooted = FALSE)
#myBoots[is.na(myBoots)] <- 0
myBoots <- c(rep(NA, 16), myBoots)
# tr <- ggtree(tree2) +
#   #theme_tree2() +
#   geom_tiplab() +
#   geom_nodelab(aes(label=myBoots), size = 3) +
#   xlim(0, 15)
# tr



# p2 <-ggtree(tree2) %<+% d +
#   geom_tiplab(aes(label=paste0('bold(',label,')~italic(', genus, ')~italic(', species,')~', '"', strain, '"')), parse=TRUE) +
#   geom_hilight(node=6, fill="#800020", alpha=0.2, extend=5) +
#   geom_nodelab(aes(label=myBoots), size = 3) +
#   geom_treescale(width=0.5, x=4.5, y=0)
# p2

##############
p2 <-ggtree(tree2) %<+% d +
  geom_tiplab(aes(label=paste0('bold(',label,')~italic(', genus, ')~italic(', species,')~', '"', strain, '"')), parse=TRUE, align=TRUE) +
  geom_hilight(node=6, fill="#800020", alpha=0.2, extend=32) +
  geom_label_repel(aes(subset = !is.na(label), label=myBoots, fill=myBoots), force=0, color="white", size=2) + scale_fill_viridis_c(name = "Bootstrap Support") +
  #geom_treescale(width=1, x=25, y=0) +
  theme(legend.position = "bottom", 
         legend.direction = "horizontal",
         legend.title = element_text(face="bold", hjust = 0.5)) #+ xlim_expand(c(0, 35), 'Dot')
  #hexpand(.4, direction = 1)
p2
     #+
  #guides(fill = guide_colourbar(title="Bootstrap Support", title.position = "top"))
# p2 + xlim_tree(1) + xlim_expand(c(0, 35), 'Dot')
# geom_label_repel(aes(subset = !isTip & !is.na(label), label=support, fill=support), force=0, color="white") + 
#   scale_fill_viridis_c()

protein_sequences <-"data/nisin_clustalw.fasta"
x <- readAAStringSet(protein_sequences)
data = tidy_msa(x)
letters <- unique(data$character)

#my_pal <- colorRampPalette(rev(brewer.pal(n = 9, name = "PuBuGn")))
my_pal <- colorRampPalette(c("#015746","#F6EDF5"))
my_custom <- data.frame(names = letters, 
                         color = my_pal(18), 
                         stringsAsFactors = FALSE)
my_custom[(which(my_custom$names=="C")),2] <- c("#FFC300")

p3 <- p2 + ggnewscale::new_scale_fill()
p4 <- p3 +#+ xlim_tree(3) +
  geom_facet(geom = geom_msa, 
             data = data,  
             panel = 'Multiple Sequence Alignment',
            # seq_name = TRUE,
             #consensus_views = TRUE, disagreement = FALSE, use_dot = TRUE, #color = 'Taylor_AA',
            color="Clustal",
            by_conversation = TRUE,
             #ref = "NisG",
            # posHighligthed = c(),
             #color = "Chemistry_AA",
             #custom_color = my_custom,
             font="DroidSansMono",
             char_width=0.5) + theme_tree2() + #+ xlim_tree(30) +
   theme(legend.position = "bottom", 
         legend.direction = "horizontal",
         legend.title = element_text(face="bold", hjust = 0.5))#+
           # guides(fill = guide_legend(title="Bootstrap Support", title.position = "top"))
p4
p4 <- p4  + xlim_tree(35) #+ xlim_expand(c(0, 35), 'Dot')
pp <- facet_labeller(p4, c(Tree = "Neighbor-joining tree", panel = "Multiple Sequence Alignment")) %>% facet_widths(widths=c(.5,1))
pp

# +
#   geom_facet(geom = ggmsa:::geom_seqlogo, data = data,  panel = 'msa', color = "Chemistry_AA", font = "helvetical", adaptive = F) + 
#   xlim_tree(1)
ggsave("figures/nisin_msa.pdf", width = 45, height = 20, units = "cm")
dev.off()
# q <-ggmsa(protein_sequences, color = "Chemistry_AA", seq_name = F, font = "DroidSansMono", char_width = 0.5) + 
#   geom_seqlogo(color = "Chemistry_AA", font = "helvetical", adaptive = F )+ theme_void()
# q
