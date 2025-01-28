# QDR RNAseq Solanum species
# Filtering GO terms of all species, keep only PPV > 0.4 
# Severin Einspanier

rm(list=ls())
library(tidyverse)

setwd("")

# Function to filter GO terms and write output files
filter_go_terms <- function(species, input_path, output_filtered, output_ids, output_bingo) {
  GO_data <- read.csv(input_path, header = T, sep = "\t")
  GO_filtered <- GO_data %>% 
    mutate(PPV = as.numeric(PPV)) %>% 
    filter(PPV > .4)
  GO_filtered_ids <- GO_filtered %>% 
    select(qpid) %>% 
    unique()
  
  write_tsv(GO_filtered, output_filtered)
  write_tsv(GO_filtered_ids, output_ids)
  
  out_file <- read.delim(gsub("anno.out", "GO.out", input_path)) %>%
    filter(ARGOT_PPV > 0.4) %>% 
    select(qpid, goid) %>% 
    mutate(gene = gsub("GeneExt~", "", qpid)) %>%  
    mutate(gene = gsub("mRNA_", "", gene)) %>%  
    mutate(gene = gsub("t\\.peak", "g.peak", gene)) %>% 
    mutate(gene = gsub("t\\.minus", "g.minus", gene)) %>% 
    mutate(gene = gsub("t\\.plus", "g.plus", gene)) %>% 
    mutate(gene = gsub("\\.[1-9].*|\\.p[1-9].*", "", gene)) %>% 
    select(gene, goid) %>% 
    mutate(goid = paste0("= ", goid))
  
  writeLines(c(paste0("(species=", species, ") (type=Consortium) (curator=GO)")), con = output_bingo)
  write.table(out_file, file = output_bingo, append = TRUE, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

# Filter GO terms for each species
filter_go_terms("Solanum lycopersicoides", "slyco/anno.out", "slyco_GO_filtered.txt", "slyco_GO_filtered_IDs.txt", "slyco_goterms.txt")
filter_go_terms("Solanum pennellii", "spen/anno.out", "spen_GO_filtered.txt", "spen_GO_filtered_IDs.txt", "spen_goterms.txt")
filter_go_terms("Solanum pimpinellifolium", "spimp/anno.out", "spim_GO_filtered.txt", "spim_GO_filtered_IDs.txt", "spimp_goterms.txt")
filter_go_terms("Solanum habrochaites", "shabr/anno.out", "shabro_GO_filtered.txt", "shabro_GO_filtered_IDs.txt", "shabro_goterms.txt")
filter_go_terms("Solanum chilense", "schil/anno.out", "schil_GO_filtered.txt", "schil_GO_filtered_IDs.txt", "schil_goterms.txt")