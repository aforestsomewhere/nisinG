#' Script to generate plots of dominance on NisinG samples
#' Developed by A. Kate Falà (https://github.com/aforestsomewhere)
#' Collapsed
#' All timepoints

library(wesanderson)
library(FSA)
####################
### Functions    ###
####################
#function to translate p values to significance codes
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
palette_time <- c("#9A8822", "#F8AFA8", "#74A089")
palette_groups <- c("#332211","#ffcc44", "#44aacc", "#bb2211")
palette_collapsed <- c("#FF0000", "#00A08A")

#palette_time <- wes_palette("Darjeeling2", n=3, type="discrete")
#strep_pal <- wes_palette("Darjeeling1", n=2, type="discrete")
#function to plot dominance
plot_dom <- function(df, metric, color_by, wrap_by, pal){
  df %>% ggplot(aes_string(x=color_by, y=metric, fill=color_by))+
    geom_violin(trim=F) + geom_point(color="black", alpha=0.8)  +
    scale_fill_manual(name="Inoculation",
                      labels=c(expression(paste("No ",italic("S. salivarius"))),
                               expression(paste("Any ",italic("S. salivarius")))), 
                      values = pal) +    
    xlab("") + theme_classic() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(face="bold")) + facet_wrap(wrap_by) +
    guides(fill="none", alpha="none", color="none")
}
#similiar but without faceting
plot_dom_no_facet <- function(df, metric, color_by, pal){
  df %>% ggplot(aes_string(x=color_by, y=metric, fill=color_by))+
    geom_violin(trim=F) +geom_point(color="black", alpha=0.8)  +
    scale_fill_manual(name="Inoculation",
                      labels=c(expression(paste("No ",italic("S. salivarius"))),
                               expression(paste("Any ",italic("S. salivarius")))), 
                      values = pal) +    
    xlab("")+ theme_classic() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(face="bold")) +
    guides(fill="none", alpha="none", color="none")
}

####################
### Main    ###
####################

df1 <- read.csv("reports/dominance_all_values.csv")
df1 <- df1[,-1]
theme_set(theme_bw())
#first plot to extract common legend
#pal1 <- c("#FF0000", "#00A08A")
df1 <- df1 %>% janitor::clean_names() %>%
  mutate(timepoint=factor(timepoint, levels=c("T_0", "T_6", "T_24"))) %>%
  mutate(strep_exposure = case_when(
    treatment == "FS" ~ "Y", 
    treatment == "S" ~ "Y",
    treatment == "F" ~ "N",
    treatment == "B" ~ "N"))

#fetch list of dominance indices
subset <- df1[,1:7]
dominance_indices <- names(subset)
rm(subset)

########################
# Testing significance #
########################

#compare treatment, strep_exposure, timepoint groups
grouping_var <- c("treatment","strep_exposure", "timepoint")
results_collapsed <- data.frame()
for (i in (1:length(dominance_indices))){
  for (j in (1:length(grouping_var))){
    result <- df1 %>%
      do(tidy(kruskal.test(x = .[[dominance_indices[i]]], g = .[[grouping_var[j]]])))
    result <- result %>% 
      mutate(metric=dominance_indices[i]) %>%
      mutate(group=grouping_var[j]) %>%
      mutate(sig_code=sigcode(p.value))
    results_collapsed <- bind_rows(results_collapsed, result)
  }
}

#follow up test to quantify effect of timepoint
c3 <- dunnTest(dbp ~ timepoint, data = df1, method="bh")
dunn_result2 <- c3$res %>%
  mutate(method="Dunn's Test") %>%
  mutate(correction="BH") %>%
  rowwise() %>%
  mutate(sig_code=sigcode(as.numeric(P.adj)))

#plot to show timepoint
legend <-ggplot(data = df1, aes(x=timepoint, y=dbp, fill=timepoint))+
  geom_violin() +geom_point(color="black", alpha=0.8)  +
  scale_fill_manual(name="Timepoint",
                    labels=c("T_0", "T_6", "T_24"),
                    values = palette_time) +
  theme(legend.box = "vertical", 
        legend.position = "bottom",
        legend.title = element_text(face="bold", ),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  guides(fill = guide_legend(title.position = "top",nrow=4, byrow=TRUE))
legend <- get_legend(legend)


#make plot for each dominance index
myplots <- vector("list", 7)

for (i in 1:length(dominance_indices)){
  myplots[[i]] <- local({
    i <- i
    p1 <- plot_dom_no_facet(df=df1, metric = dominance_indices[i], color_by = "timepoint", pal=palette_time)
    print(p1)
  })
}

final_plot <- plot_grid(nrow=2, 
                        plotlist = myplots, legend, 
                        labels = c('','A','B','C','D','E','F','G'))
final_plot
ggsave("figures/dominance_collapsed_bytimepoint.pdf", width=29, height=21, dpi=300, units="cm")



#testing within timepoint groups
grouping_var <- c("strep_exposure", "dominant_species")
results_collapsed <- data.frame()
for (i in (1:length(indices))){
  for (j in (1:length(grouping_var))){
    result <- df1 %>%
      group_by(timepoint) %>%
      do(tidy(kruskal.test(x = .[[indices[i]]], g = .[[grouping_var[j]]])))
    result <- result %>% 
      mutate(metric=indices[i]) %>%
      mutate(group=grouping_var[j]) %>%
      mutate(sig_code=sigcode(p.value))
    results_collapsed <- bind_rows(results_collapsed, result)
  }
}
#dominant species approaches significance level - follow up
c2 <- dunnTest(dbp ~ dominant_species, data = dom_big, method="bh")
dunn_result2 <- c2$res %>%
  mutate(method="Dunn's Test") %>%
  mutate(correction="BH") %>%
  rowwise() %>%
  mutate(sig_code=sigcode(as.numeric(P.adj)))


legend <-ggplot(data = df1, aes(x=timepoint, y=dbp, fill=timepoint))+
  geom_violin() +geom_point(color="black", alpha=0.8)  +
  scale_fill_manual(name="Inoculation",
                    labels=c(expression(paste("No ",italic("S. salivarius"))),
                             expression(paste("Any ",italic("S. salivarius")))), 
                    values = wes_palette("Darjeeling1", n = 4, type = "discrete")) +
  theme(legend.box = "vertical", 
        legend.position = "bottom",
        legend.title = element_text(face="bold", ),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + facet_wrap(~timepoint) +
  guides(fill = guide_legend(title.position = "top",nrow=4, byrow=TRUE))
legend <- get_legend(legend)


#make plot for each dominance index
myplots <- vector("list", 7)

for (i in 1:length(dominance_indices)){
  myplots[[i]] <- local({
    i <- i
    p1 <- plot_dom(df=df1, metric = dominance_indices[i], color_by = "strep_exposure", wrap_by="timepoint")
    print(p1)
  })
}

final_plot <- plot_grid(nrow=4, 
                        plotlist = myplots, legend, 
                        labels = c('','A','B','C','D','E','F','G'))
final_plot
ggsave("figures/dominance_collapsed.pdf", width=21, height=29, dpi=300, units="cm")


#######################
#generate output file #
#######################
list_of_datasets <- list("Dominance Indices" = dom_big,
                         "Kruskal Wallis, all groups" = results_uncollapsed,
                         "Dunn Testing, Timepoint" = dunn_result,
                         "Kruskal Wallis, collapsed" = results_collapsed,
                         "Dunn Testing, dominant species" = dunn_result2)
openxlsx::write.xlsx(list_of_datasets, file = "reports/dominance.xlsx")




########################
# Generating plots     #
########################
  
legend <-ggplot(data = df1, aes(x=strep_exposure, y=dbp, fill=strep_exposure))+
  geom_violin(trim=F) +geom_point(color="black", alpha=0.8)  +
  scale_fill_manual(name="Inoculation",
                    labels=c(expression(paste("No ",italic("S. salivarius"))),
                             expression(paste("Any ",italic("S. salivarius")))), 
                    values = wes_palette("Darjeeling1", n = 4, type = "discrete")) + theme_classic() +
  theme(legend.box = "vertical", 
        legend.position = "bottom",
        legend.title = element_text(face="bold", ),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + facet_wrap(~timepoint) +
  guides(fill = guide_legend(title.position = "top",nrow=4, byrow=TRUE))
legend <- get_legend(legend)


#fetch list of dominance indices
subset <- df1[,1:7]
dominance_indices <- names(subset)
rm(subset)
#make plot for each dominance index
myplots <- vector("list", 7)

for (i in 1:length(dominance_indices)){
  myplots[[i]] <- local({
    i <- i
    p1 <- plot_dom(df=df1, metric = dominance_indices[i], color_by = "strep_exposure", wrap_by="timepoint", pal=strep_pal)
    print(p1)
  })
}

final_plot <- plot_grid(nrow=4, 
                        plotlist = myplots, legend, 
                        labels = c('','A','B','C','D','E','F','G'))
final_plot
ggsave("figures/dominance_collapsed_alltimes.pdf", width=21, height=29, dpi=300, units="cm")
