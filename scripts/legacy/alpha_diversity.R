#' Final plots P4 (T_24 only, collapsing B with F and FS with S
#' BC, TSNE, UNIFRAC, ADONIS, ANOSIM, BETA DISP
#' Ascertain S and FS are more similiar, and likewise FS/B
#' Then (only then) collapse
#' 22/03/23

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
library(rbiom)
library(phyloseq)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)
library(wesanderson)
library(ggrepel)
library(car)
library(Rtsne)
library(tsnemicrobiota)
library(phylosmith)
library(microbiome)
library(pairwiseAdonis)
library(ggpubr)
library(rstatix)
library(ggsignif)

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
mpa_table <- read.table("data\\metaphlan.tsv", comment.char = '#', sep = '\t', header = TRUE)
mpa_table <- mpa_table[grep('s__',mpa_table[,1]),]
mpa_table[,1] <- gsub(".+\\|s__", "", mpa_table[,1])
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]
info <- read.csv("data\\enriq_metadata.csv", nrows=62, row.names=1)
names(info)[names(info) == 'experiement'] <- "experiment"
info$treatment_replicate <- info$condition #split condition columns for subsequent grouping
info %>% separate(condition, c("treatment", "replicate"), "_") %>% 
  mutate(exp_group = paste(treatment, timepoint, sep ='_')) %>%
  mutate(treatment = str_replace(treatment, "FS", "S")) %>%
  mutate(treatment = str_replace(treatment, "F", "B")) %>%
  filter(experiment=="nisin") %>%
  filter(timepoint=="T_24") -> info #filter only nisin experiment data, only T_24
rownames(info) -> nisin_sam
mpa_table[,colnames(mpa_table) %in% nisin_sam] -> mpa_table
#mpa_table <- mpa_table[rowSums(mpa_table) >= 0.000001,]

bug_mat = t(as.matrix(mpa_table))
bug_mat <- bug_mat[sort(rownames(bug_mat)),] #order samples (rows) to properly merge with infO
bug_mat_t <- t(bug_mat)

#ALPHA
richness_vec <-(specnumber(bug_mat)) #ok
shannon_vec <- vegan::diversity(bug_mat, index = "shannon")
simp_vec <- vegan::diversity(bug_mat, index = "simpson")
inv_simp_vec <- vegan::diversity(bug_mat, index = "invsimpson")
tidy_df <- cbind.data.frame(comb_id = rownames(info), richness = richness_vec, shannon = shannon_vec, simpson = simp_vec, inv_simpson = inv_simp_vec, meta = info)
tidy_df$meta.treatment = factor(tidy_df$meta.treatment)

# #HISTOGRAMS FOR NORMALITY
# par(mfrow=c(2,2))
# hist(tidy_df$shannon)
# hist(tidy_df$richness)
# hist(tidy_df$simpson)
# hist(tidy_df$inv_simpson)
# #TREATMENT=B
# tidy_b <- tidy_df %>%
#   filter(meta.treatment=="B")
# par(mfrow=c(2,2))
# hist(tidy_b$shannon)
# hist(tidy_b$richness)
# hist(tidy_b$simpson)
# hist(tidy_b$inv_simpson)
# #TREATMENT=S
# tidy_s <- tidy_df %>%
#   filter(meta.treatment=="S")
# par(mfrow=c(2,2))
# hist(tidy_s$shannon)
# hist(tidy_s$richness)
# hist(tidy_s$simpson)
# hist(tidy_s$inv_simpson)

#########Alpha boxplots
alpha_all <- data.frame(Data="MetaPhlAn4", 
                        Metric="Metric",
                        Comparing="Comparing",
                        Test="Wilcox",
                        Test_Statistic="F",
                        P_value="P",
                        p_adjust_method="BH",
                        Significance_codes="NA") #initialise output DF
dummy_legend <- ggplot(data = tidy_df, aes(x=meta.treatment, y=richness, fill=meta.treatment))+
  geom_boxplot() + geom_jitter() 
 # scale_fill_manual(name="Treatment",
                    #labels=c(expression(paste("No ", italic("S. salivarius")," exposure")),
                    #         expression(paste(italic("S. salivarius")," treated"))),
                    #values = wes_palette("Darjeeling1", n = 2, type = "discrete"))+
 # theme(legend.position = "bottom", legend.box="horizontal",
  #      legend.text = element_text(size=12),
  #      legend.title = element_text(size=14))
dummy_legend
legend <- cowplot::get_legend(dummy_legend)


#Richness boxplot
bxp <- ggboxplot(tidy_df, x = "meta.treatment", y = "richness", 
                 fill = "meta.treatment", xlab="", ylab="")
res.stat <- tidy_df %>% 
  rstatix::wilcox_test(data = ., richness ~ meta.treatment) %>% 
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "meta.treatment", fun = "max") 
#compile results to report
alpha_df <- c("MetaPhlAn4", res.stat$.y., paste(res.stat$group1," v.s ", res.stat$group2), "Wilcoxon rank sum",res.stat$statistic,res.stat$p.adj,"BH",res.stat$p.adj.signif)
alpha_all <-rbind(alpha_all, alpha_df)
p1<- bxp + stat_pvalue_manual(res.stat, label = "p.adj = {p.adj}", tip.length = 0.01)+
  scale_fill_manual(name="Treatment",
                    labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = wes_palette("Darjeeling1", n = 2, type = "discrete"))+
  geom_jitter() + ggtitle("Observed Species") + theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank())
p1
#Shannon boxplot
bxp <- ggboxplot(tidy_df, x = "meta.treatment", y = "shannon", 
                 fill = "meta.treatment", xlab="", ylab="")
res.stat <- tidy_df %>% 
  rstatix::wilcox_test(data = ., shannon ~ meta.treatment) %>% 
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "meta.treatment", fun = "max") 
#compile results to report
alpha_df <- c("MetaPhlAn4", res.stat$.y., paste(res.stat$group1," v.s ",res.stat$group2), "Wilcoxon rank sum",res.stat$statistic,res.stat$p.adj,"BH",res.stat$p.adj.signif)
alpha_all <-rbind(alpha_all, alpha_df)
p2<- bxp + stat_pvalue_manual(res.stat, label = "p.adj = {p.adj}", tip.length = 0.01)+
  scale_fill_manual(name="Treatment",
                    labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = wes_palette("Darjeeling1", n = 2, type = "discrete"))+
  geom_jitter() + ggtitle("Shannon") + theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank())
p2

#Simpson boxplot
bxp <- ggboxplot(tidy_df, x = "meta.treatment", y = "simpson", 
                 fill = "meta.treatment", xlab="", ylab="")
res.stat <- tidy_df %>% 
  rstatix::wilcox_test(data = ., simpson ~ meta.treatment) %>% 
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "meta.treatment", fun = "max") 
alpha_df <- c("MetaPhlAn4", res.stat$.y., paste(res.stat$group1," v.s ",res.stat$group2), "Wilcoxon rank sum",res.stat$statistic,res.stat$p.adj,"BH",res.stat$p.adj.signif)
alpha_all <-rbind(alpha_all, alpha_df)
p3<- bxp + stat_pvalue_manual(res.stat, label = "p.adj = {p.adj}", tip.length = 0.01)+
  scale_fill_manual(name="Treatment",
                    labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = wes_palette("Darjeeling1", n = 2, type = "discrete"))+
  geom_jitter() + ggtitle("Simpson") + theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank())
p3
#############################################
#Inverse Simpson
bxp <- ggboxplot(tidy_df, x = "meta.treatment", y = "inv_simpson", 
                 fill = "meta.treatment", xlab="", ylab="")
res.stat <- tidy_df %>% 
  rstatix::wilcox_test(data = ., inv_simpson ~ meta.treatment) %>% 
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "meta.treatment", fun = "max") 
alpha_df <- c("MetaPhlAn4", res.stat$.y., paste(res.stat$group1," v.s ",res.stat$group2), "Wilcoxon rank sum",res.stat$statistic,res.stat$p.adj,"BH",res.stat$p.adj.signif)
alpha_all <-rbind(alpha_all, alpha_df)
p4<- bxp + stat_pvalue_manual(res.stat, label = "p.adj = {p.adj}", tip.length = 0.01)+
  scale_fill_manual(name="Treatment",
                    labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = wes_palette("Darjeeling1", n = 2, type = "discrete"))+
  geom_jitter() + ggtitle("Inverse simpson") + theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank())

plots <- plot_grid(p1,p2,p3,p4, nrow=1, align='v')
#plots
ggdraw(plot_grid(plots, legend, 
                 nrow=2, align='v', rel_heights=c(1, 0.2),
                 rel_widths=c(1, 2))) + theme(plot.background = element_rect(fill="#FFFFFF", color = NA))
ggsave("figures\\alpha_boxplot_sig_T24.emf", width = 10, height = 10, dpi = 600, device = {function(filename, ...) devEMF::emf(file = filename, ...)})

#outputstats
alpha_all<- alpha_all[-1,]
write.csv(alpha_all, "reports\\alpha_tested_T24.csv", row.names = FALSE)


















######################################################################
#stats testing alpha

alpha_all <- data.frame(Data="MetaPhlAn4", 
                        Metric="Shannon",
                        Comparing=NA,
                        Test="Wilcox",
                        Test_Statistic=NA,
                        P_value=NA,
                        p_adjust_method=NA,
                        Significance_codes="NA") #initialise output DF
#shannon
wilcox <- wilcox.test(tidy_df$shannon, tidy_df$treatment, p.adjust.method = "BH")
#anosim_treat <- anosim(dist_mat,info$treatment, permutations = 9999) #TIMEPOINT
P<- wilcox$p.value
F<-as.numeric(wilcox$statistic)
S<-sigcode(P)
T<-wilcox$method
M<- gsub(".*[$]","",wilcox$data.name) #name of diversity metric
G<- paste(as.character(levels(tidy_df$meta.treatment))[1],"v.s",as.character(levels(tidy_df$meta.treatment))[2])
alpha_df <- c("MetaPhlAn4", M, G, T,F,P,"BH",S)
alpha_all <-rbind(alpha_all, alpha_df)
alpha_all <- alpha_all[-1,]  #remove first helper shaper row

#richness
wilcox <- wilcox.test(tidy_df$richness, tidy_df$treatment, p.adjust.method = "BH")
P<- wilcox$p.value
F<-as.numeric(wilcox$statistic)
S<-sigcode(P)
T<-wilcox$method
M<- gsub(".*[$]","",wilcox$data.name) #name of diversity metric
G<- paste(as.character(levels(tidy_df$meta.treatment))[1],"v.s",as.character(levels(tidy_df$meta.treatment))[2])
alpha_df <- c("MetaPhlAn4", M, G, T,F,P,"BH",S)
alpha_all <-rbind(alpha_all, alpha_df)

#simpson
wilcox <- wilcox.test(tidy_df$simpson, tidy_df$treatment, p.adjust.method = "BH")
P<- wilcox$p.value
F<-as.numeric(wilcox$statistic)
S<-sigcode(P)
T<-wilcox$method
M<- gsub(".*[$]","",wilcox$data.name) #name of diversity metric
G<- paste(as.character(levels(tidy_df$meta.treatment))[1],"v.s",as.character(levels(tidy_df$meta.treatment))[2])
alpha_df <- c("MetaPhlAn4", M, G, T,F,P,"BH",S)
alpha_all <-rbind(alpha_all, alpha_df)



















dev.off()

#plot(legend)
p1<- ggplot(data = tidy_df, aes(x=meta.treatment, y=richness, fill=meta.treatment))+
  geom_boxplot() + geom_jitter() + ggtitle("Observed \n Species")+
  scale_fill_manual(name="Treatment",
                    labels=c('Control','Strep'),
                    values = wes_palette("Darjeeling1", n = 2, type = "discrete"))+
  theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("") +ylab("")
p2 <- ggplot(data = tidy_df, aes(x=meta.treatment, y=shannon, fill=meta.treatment))+
  geom_boxplot() + geom_jitter() + ggtitle("Shannon")+
  scale_fill_manual(name="Treatment",
                    labels=c('Control','Strep'),
                    values = wes_palette("Darjeeling1", n = 2, type = "discrete"))+
  theme(legend.position='none',axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("") +ylab("")
p3 <-ggplot(data = tidy_df, aes(x=meta.treatment, y=simp_vec, fill=meta.treatment))+
  geom_boxplot() + geom_jitter() + ggtitle("Simpson")+
  scale_fill_manual(name="Treatment",
                    labels=c('Control','Strep'),
                    values = wes_palette("Darjeeling1", n = 2, type = "discrete"))+
  theme(legend.position='none',axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("") +ylab("")
p4 <-ggplot(data = tidy_df, aes(x=meta.treatment, y=inv_simpson, fill=meta.treatment))+
  geom_boxplot() + ggtitle("Inverse \nSimpson") + geom_jitter()+
  #stat_pvalue_manual(rstatix::wilcox_test(inv_simpson ~ meta.treatment),
  #                  rstatix::add_xy_position()) +
  scale_fill_manual(name="Treatment",
                    labels=c(expression(paste("No ", italic("S. salivarius"))),expression(italic("S. salivarius"))),
                    values = wes_palette("Darjeeling1", n = 2, type = "discrete"))+
  theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  xlab("") +ylab("")
p4









######################################################
wilcox.test(richness ~ meta.treatment, data = tidy_df)
wilcox.test(shannon ~ meta.treatment, data = tidy_df)
wilcox.test(simpson ~ meta.treatment, data = tidy_df)
wilcox.test(inv_simpson ~ meta.treatment, data = tidy_df)
t.test(richness ~ meta.treatment, data =tidy_df)
t.test(shannon ~ meta.treatment, data =tidy_df)
t.test(simpson ~ meta.treatment, data =tidy_df)
t.test(inv_simpson ~ meta.treatment, data =tidy_df)

ggsave("figures\\ALPHA_boxplots_n3_T24.emf", width = 10, height = 10, dpi = 600, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
