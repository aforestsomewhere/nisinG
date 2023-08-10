#alpha all like dom


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
time_pal <- wes_palette("Darjeeling2", n=3, type="discrete")
strep_pal <- wes_palette("Darjeeling1", n=2, type="discrete")
#function to plot dominance
plot_dom <- function(df, metric, color_by, wrap_by, pal){
  df %>% ggplot(aes_string(x=color_by, y=metric, fill=color_by))+
    geom_violin() +geom_point(color="black", alpha=0.8)  +
    scale_fill_manual(name="Inoculation",
                      labels=c(expression(paste("No ",italic("S. salivarius"))),
                               expression(paste("Any ",italic("S. salivarius")))), 
                      values = pal) +    
    xlab("")+ 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(face="bold")) + facet_wrap(wrap_by) +
    guides(fill="none", alpha="none", color="none")
}
#similiar but without faceting
plot_dom_no_facet <- function(df, metric, color_by, pal){
  df %>% ggplot(aes_string(x=color_by, y=metric, fill=color_by))+
    geom_violin() +geom_point(color="black", alpha=0.8)  +
    scale_fill_manual(name="Inoculation",
                      labels=c(expression(paste("No ",italic("S. salivarius"))),
                               expression(paste("Any ",italic("S. salivarius")))), 
                      values = pal) +    
    xlab("")+ 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(face="bold")) +
    guides(fill="none", alpha="none", color="none")
}

####################
### Main    ###
####################

df1 <- read.csv("reports/alpha_all_metrics.csv")
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

#generate list of alpha diversity metrics
indices <- c("richness","shannon","simpson", "inv_simpson")

####################################
# Testing significance uncollapsed #
####################################

#compare treatment, strep_exposure, timepoint groups
grouping_var <- c("treatment","strep_exposure", "timepoint")
results_collapsed <- data.frame()
for (i in (1:length(indices))){
  for (j in (1:length(grouping_var))){
    result <- df1 %>%
      do(tidy(kruskal.test(x = .[[indices[i]]], g = .[[grouping_var[j]]])))
    result <- result %>% 
      mutate(metric=indices[i]) %>%
      mutate(group=grouping_var[j]) %>%
      mutate(sig_code=sigcode(p.value))
    results_collapsed <- bind_rows(results_collapsed, result)
  }
}

#follow up test to quantify effect of timepoint
c3 <- dunnTest(shannon ~ timepoint, data = df1, method="bh")
dunn_result2 <- c3$res %>%
  mutate(method="Dunn's Test") %>%
  mutate(correction="BH") %>%
  rowwise() %>%
  mutate(sig_code=sigcode(as.numeric(P.adj)))

####################################
# Testing significance collapsed   #
####################################
#can group and check within timepoints with n=6 for collapsed (strep_exposure)
#testing within timepoint groups
grouping_var <- c("strep_exposure")
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

