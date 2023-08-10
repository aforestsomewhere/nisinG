#' Script to generate plots of dominance on NisinG samples
#' Developed by A. Kate Falà (https://github.com/aforestsomewhere)
#' Uncollapsed, showing all groups 
#' All timepoints
#' l

library(cowplot)
library(ggplot2)
library(dplyr)
palette_time <- c("#9A8822", "#F8AFA8", "#74A089")
palette_groups <- c("#332211","#ffcc44", "#44aacc", "#bb2211")
palette_collapsed <- c("#FF0000", "#00A08A")

df1 <- read.csv("reports/dominance_all_values.csv")
df1 <- df1[,-1]
# set_theme(theme_bw())
# ggplot2::theme_set(theme_bw())
#first plto to extract common legend
#palette_groups <- c("#FF0000", "#00A08A")
df1 <- df1 %>% janitor::clean_names() %>%
  mutate(timepoint=factor(timepoint, levels=c("T_0", "T_6", "T_24")))

# old <- theme_set(theme_classic())
# theme_set(old)
legend <-ggplot(data = df1, aes(x=treatment, y=dbp, fill=treatment))+
  geom_violin() +geom_point(color="black", alpha=0.8)  +
  scale_fill_manual(name="Treatment group",
                    labels=c('Control', 
                             expression(paste(italic("F. nucleatum"))), 
                             expression(paste(italic("F. nucleatum/S. salivarius"))), 
                             expression(paste(italic('S. salivarius')))),
                    values = palette_groups) + theme_classic() +
  theme(legend.box = "vertical", 
        legend.position = "bottom",
        legend.title = element_text(face="bold", ),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + facet_wrap(~timepoint) +
  guides(fill = guide_legend(title.position = "top",nrow=4, byrow=TRUE))
legend <- get_legend(legend)

#function to plot dominance
plot_dom <- function(df, metric, color_by, wrap_by){
  df %>% ggplot(aes_string(x=color_by, y=metric, fill=color_by))+
    geom_violin(trim=F) +geom_point(color="black", alpha=0.8, size=.5)  +
    scale_fill_manual(name="Inoculation",
                      labels=c(expression(paste("No ",italic("S. salivarius"))),
                               expression(paste("Any ",italic("S. salivarius")))),                       
                      values = palette_groups) + 
    xlab("") + theme_classic() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(face="bold")) + facet_wrap(wrap_by) +
    guides(fill="none", alpha="none", color="none")
}
#fetch list of dominance indices
subset <- df1[,1:7]
dominance_indices <- names(subset)
rm(subset)
#make plot for each dominance index
myplots <- vector("list", 7)

for (i in 1:length(dominance_indices)){
  myplots[[i]] <- local({
    i <- i
    p1 <- plot_dom(df=df1, metric = dominance_indices[i], color_by = "treatment", wrap_by="timepoint")
    print(p1)
  })
}
# old <- theme_set(theme_classic())
# theme_set(old)
final_plot <- plot_grid(nrow=4, 
                        plotlist = myplots, legend, 
                        labels = c('','A','B','C','D','E','F','G'))
final_plot
ggsave("figures/dominance_uncollapsed.pdf", width=21, height=29, dpi=300, units="cm")

#save individual plot of dominance across timepoints for main paper
#df_dbp <- df1 %>% select(db
dbp_plot <- plot_dom(df=df1, metric = "dbp", color_by="treatment", wrap_by="timepoint")
dbp_plot + geom_point(size=2) + 
  ylab("DBP Berger–Parker index") +
  theme(axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=14),
        strip.text.x = element_text(size = 24, face="bold"))
ggsave("figures/dominance_alltimepoints.pdf", width=29, height=12, dpi=300, units="cm")
ggsave("figures/dominance_alltimepoints.png", width=29, height=12, dpi=600, units="cm")

