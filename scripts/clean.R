#Script to clean and subset metadata and MPA4 output to NisinG samples
#Developed by A. Kate Falà (https://github.com/aforestsomewhere)

#metadata file
meta <- read.csv("data\\enriq_metadata.csv", nrows=62, row.names=1)
#correct typo
names(meta)[names(meta) == 'experiement'] <- "experiment"
#split condition columns for subsequent grouping
meta$treatment_replicate <- meta$condition
#clean table and subset to NisinG samples
meta %>% separate(condition, c("treatment", "replicate"), "_") %>% 
  mutate(exp_group = paste(treatment, timepoint, sep ='_')) %>%
  filter(experiment=="nisin") %>%
  select(-experiment) -> meta
#save NisinG sample names
rownames(meta) -> nisin_sam
#Read MetaPhlAn4 table
mpa_table <- read.table("data\\metaphlan.tsv", comment.char = '#', sep = '\t', header = TRUE)
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]
#Subset to NisinG samples
mpa_table[,colnames(mpa_table) %in% nisin_sam] -> mpa_table
#write outputs
write.csv(meta, "cleaned_data/metadata.csv")
write.table(mpa_table, "cleaned_data/mpa_table.tsv", sep = '\t')
