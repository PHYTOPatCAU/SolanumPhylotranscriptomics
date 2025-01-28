## QDR RNAseq Solanum species
## Generate a GO-annotation for OGs 
## Approach: 
## 1. Infer gene ID and OG ID
## 2. Read the GO files of the five species
## 3. Filter the GO terms based on PPV > .4
## 4. Merge the GO terms of the five species and keep the shared GO terms (? might get complicated)
## 5. Write the merged GO terms to a file
# Severin Einspanier

rm(list=ls())
library(tidyverse)

##### 1. Infer gene ID and OG ID #####

sng_cp_ogs <- read.csv("orthogroups_singlecopy_goid.csv")

# Function to process OG data
process_og_data <- function(filename) {
  read.csv(filename) %>% 
    filter(OG %in% sng_cp_ogs$OG) %>% 
    filter(GeneID != "") %>% 
    mutate(gene = gsub("GeneExt~", "", GeneID)) %>%  
    mutate(gene = gsub("mRNA_", "", gene)) %>%  
    mutate(gene = gsub("t\\.peak", "g.peak", gene)) %>% 
    mutate(gene = gsub("t\\.minus", "g.minus", gene)) %>% 
    mutate(gene = gsub("t\\.plus", "g.plus", gene)) %>% 
    mutate(gene = gsub("\\.[1-9].*|\\.p[1-9].*", "", gene))
}

OG_data_slyco <- process_og_data("genid2goid_slyco.csv")
OG_data_spen <- process_og_data("genid2goid_spen.csv")
OG_data_shabro <- process_og_data("genid2goid_shabro.csv")
OG_data_spimp <- process_og_data("genid2goid_spimp.csv")
OG_data_schil <- process_og_data("genid2goid_schil.csv")

#### 2, 3. Read the GO files of the five species and filter the GO terms based on PPV > .4 ####

# Function to process GO terms
process_go_terms <- function(filename, og_data) {
  read.table(filename, header = TRUE, sep = "\t", 
             colClasses = c("character", "character", "character", "character", "numeric", "numeric", "numeric", "numeric")) %>% 
    filter(ARGOT_PPV > 0.4) %>% 
    select(qpid, goid) %>% 
    mutate(gene = gsub("GeneExt~", "", qpid)) %>%  
    mutate(gene = gsub("mRNA_", "", gene)) %>%  
    mutate(gene = gsub("t\\.peak", "g.peak", gene)) %>% 
    mutate(gene = gsub("t\\.minus", "g.minus", gene)) %>% 
    mutate(gene = gsub("t\\.plus", "g.plus", gene)) %>% 
    mutate(gene = gsub("\\.[1-9].*|\\.p[1-9].*", "", gene)) %>% 
    right_join(og_data, by = c("gene" = "gene"))
}

GOterms_slyco <- process_go_terms("slyco_GO.out", OG_data_slyco)
GOterms_spen <- process_go_terms("spen_GO.out", OG_data_spen)
GOterms_shabro <- process_go_terms("shabr_GO.out", OG_data_shabro)
GOterms_spimp <- process_go_terms("spimp_GO.out", OG_data_spimp)
GOterms_schil <- process_go_terms("schil_GO.out", OG_data_schil)

#### 4. Merge the GO terms of the five species and keep the shared GO terms ####

GOterms_all <- GOterms_slyco %>% 
  bind_rows(GOterms_spen) %>% 
  bind_rows(GOterms_shabro) %>% 
  bind_rows(GOterms_spimp) %>% 
  bind_rows(GOterms_schil) %>% 
  select(OG, goid) %>% 
  group_by(OG, goid) %>% 
  drop_na(goid) %>% 
  summarise(count = n(), goid = goid) %>% 
  unique() %>% 
  filter(count > 2)

# How many GO terms has each OG? 
(count_GO_per_OG <- GOterms_all %>% 
  group_by(OG) %>% 
  summarise(n_GO = n()) %>% 
  arrange(desc(n_GO)) %>% 
  ggplot(aes(x = n_GO)) +
  geom_histogram(binwidth = 1)
)

### 5. Write the merged GO terms to a file ####

# Select the required columns and modify the goid column
out_file <- GOterms_all %>% 
  select(OG, goid) %>% 
  mutate(goid = paste0("= ", goid))

# Define the output file path
output_file_path <- "OG_GO_terms.txt"

# Write the header lines to the file
writeLines(c("(species=SolanumMix) (type=Consortium) (curator=GO)"), con = output_file_path)

# Append the selected data to the file
write.table(out_file, file = output_file_path, append = TRUE, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)