#basic rel abun

mpa_table <- read.table("cleaned_data/mpa_table.tsv", comment.char = '#', sep = '\t', header = TRUE, row.names=NULL)
#filter to species level assignments
mpa_table <- mpa_table[grep('s__',mpa_table[,1]),]
mpa_table[,1] <- gsub(".+\\|s__", "", mpa_table[,1])
mpa_table <- mpa_table %>% filter(grepl('t__',row.names))
rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[,-1]
#generate matrix of abundance data
#remove empty rows
mpa_table_filt <- mpa_table[rowSums(mpa_table) >= 0.000001,]

# Get row names of the matrix
row_names <- rownames(mpa_table_filt)
# Filter rows with partial string matches
partial_match_rows <- row_names[grepl("Fusobacterium", row_names, ignore.case = TRUE)]
partial_match_rows_S <- row_names[grepl("Streptococcus", row_names, ignore.case = TRUE)]
partial_match_rows_B <- row_names[grepl("Bifidobacterium_adolescentis", row_names, ignore.case = TRUE)]

# Subset the matrix based on the filtered row names
fuso <- mpa_table_filt[partial_match_rows, , drop = FALSE]
strep <- mpa_table_filt[partial_match_rows_S, , drop = FALSE]
bif <- mpa_table_filt[partial_match_rows_B, , drop = FALSE]
#strep <- t(strep)

#FUSO
#rename t__ level to species level assignments
sp1 <- as.data.frame(row.names(fuso))
names(sp1) <- c("id")
sp1 <- sp1 %>%
  separate(id, c("species","sgb"), "__") %>%
  mutate(taxon = gsub("\\|t","", species))

rownames(fuso) <- sp1$taxon
fuso <- t(fuso)
#STREP
#rename t__ level to species level assignments
sp1 <- as.data.frame(row.names(strep))
names(sp1) <- c("id")
sp1 <- sp1 %>%
  separate(id, c("species","sgb"), "__") %>%
  mutate(taxon = gsub("\\|t","", species))
rownames(strep) <- sp1$taxon
strep <- t(strep)
#BIF
sp1 <- as.data.frame(row.names(bif))
names(sp1) <- c("id")
sp1 <- sp1 %>%
  separate(id, c("species","sgb"), "__") %>%
  mutate(taxon = gsub("\\|t","", species))
rownames(bif) <- sp1$taxon
bif <- t(bif)

meta <- read.csv("cleaned_data/metadata.csv", nrows=36, row.names=1)
meta <- meta  %>%
  mutate(strep_exposure = case_when(
    treatment == "FS" ~ "Y", 
    treatment == "S" ~ "Y",
    treatment == "F" ~ "N",
    treatment == "B" ~ "N"))

df <- merge(meta, fuso, by = 'row.names')
df <- df %>% select(exp_group, Fusobacterium_nucleatum) %>% 
  group_by(exp_group) %>%
  summarise('Mean Relative Abundance(n=3)' = round(mean(Fusobacterium_nucleatum, na.rm = TRUE), digits = 2),
            S.d. = round(sd(Fusobacterium_nucleatum, na.rm = TRUE), digits=3)) %>%
  mutate(Species="Fusobacterium_nucleatum") 
df2 <- df %>% 
  separate(exp_group, c("treatment","t","timepoint"), "_") %>%
  group_by(timepoint) %>%
  mutate(timepoint = paste0("T_", timepoint)) %>%
  ungroup() %>% select(-t) %>% 
  mutate(treatment=case_when(
    treatment =="B" ~"Control",
    treatment == "F" ~ "F. nucleatum DSM15643",
    treatment == "FS" ~"F. nuc DSM15643 + S. sal DPC6487",
    treatment == "S" ~ "S. salivarius DPC6487")) %>%
  select(5,1,2,3,4)

write.csv(df2,"reports/abundance_f_nucleatum.csv")

######################################

ggplot(df, aes(x = treatment_replicate, y= Streptococcus_salivarius)) + geom_bar(stat="identity") +
  facet_wrap(~timepoint)

df %>% group_by(exp_group) %>% mutate(mean=mean(Streptococcus_salivarius)) -> test
df %>% filter(timepoint=="T_24") %>% group_by(strep_exposure) %>% mutate(mean=mean(Streptococcus_salivarius)) %>%
  mutate(sd=sd(Streptococcus_salivarius)) -> test2

df_bif <- merge(meta, bif, by='row.names')
df_bif %>% filter(timepoint=="T_24") %>% group_by(strep_exposure) %>% mutate(mean=mean(Bifidobacterium_adolescentis )) %>%
  mutate(sd=sd(Bifidobacterium_adolescentis)) -> test2
