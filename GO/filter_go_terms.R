# Filtering GO terms of all species 

rm(list=ls())
library(tidyverse)

setwd("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/pannzer2_out/manual_pannzer2/")
#slyco
GO_slyco <- read.csv("slyco/anno.out", header = T, sep = "\t")
GO_slyco_filtered <- GO_slyco %>% 
  mutate(PPV = as.numeric(PPV)) %>% 
  filter(PPV > .4)
GO_slyco_filtered_ids <- GO_slyco_filtered %>% 
  select(qpid) %>% 
  unique() 


write_tsv(GO_slyco_filtered, "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/GO/data/slyco_GO_filtered.txt")
write_tsv(GO_slyco_filtered_ids, "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/GO/data/slyco_GO_filtered_IDs.txt")


# for bingo 

out_file <- read.delim("slyco/GO.out") %>%
  filter(ARGOT_PPV >0.4) %>% 
  select(qpid, goid) %>%
  mutate(gene=gsub("GeneExt~", "", qpid))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>% 
  select(gene, goid) %>% 
  mutate(goid = paste0("= ", goid))

# Define the output file path
output_file_path <- "C:/Users/suaph281/Desktop/nesh_local/LDT_RNAseq/pannzer2_out/FILTERED_GO_TERMS/slyco_goterms.txt"

# Write the header lines to the file
writeLines(c("(species=Solanum lycopersicoides) (type=Consortium) (curator=GO)"), con = output_file_path)

# Append the selected data to the file
write.table(out_file, file = output_file_path, append = TRUE, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


#spen 
GO_spen <- read.csv("spen/anno.out", header = T, sep = "\t")
GO_spen_filtered <- GO_spen %>% 
  mutate(PPV = as.numeric(PPV)) %>% 
  filter(PPV > .4)
GO_spen_filtered_ids <- GO_spen_filtered %>% 
  select(qpid) %>% 
  unique() 
# should I remove isoforms? 

write_tsv(GO_spen_filtered, "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/GO/data/spen_GO_filtered.txt")
write_tsv(GO_spen_filtered_ids, "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/GO/data/spen_GO_filtered_IDs.txt")

# for bingo 

out_file <- read.delim("spen/GO.out") %>%
  filter(ARGOT_PPV >0.4) %>% 
  select(qpid, goid) %>%
  mutate(gene=gsub("GeneExt~", "", qpid))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>% 
  select(gene, goid) %>% 
  mutate(goid = paste0("= ", goid))

# Define the output file path
output_file_path <- "C:/Users/suaph281/Desktop/nesh_local/LDT_RNAseq/pannzer2_out/FILTERED_GO_TERMS/spen_goterms.txt"

# Write the header lines to the file
writeLines(c("(species=Solanum pennellii) (type=Consortium) (curator=GO)"), con = output_file_path)

# Append the selected data to the file
write.table(out_file, file = output_file_path, append = TRUE, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)




#spimp

GO_spimp <- read.csv("spimp/anno.out", header = T, sep = "\t")
GO_spimp_filtered <- GO_spimp %>% 
  mutate(PPV = as.numeric(PPV)) %>% 
  filter(PPV > .4)
GO_spimp_filtered_ids <- GO_spimp_filtered %>% 
  select(qpid) %>% 
  unique()

write_tsv(GO_spimp_filtered, "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/GO/data/spim_GO_filtered.txt")
write_tsv(GO_spimp_filtered_ids, "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/GO/data/spim_GO_filtered_IDs.txt")

# for bingo 

out_file <- read.delim("spimp/GO.out") %>%
  filter(ARGOT_PPV >0.4) %>% 
  select(qpid, goid) %>%
  mutate(gene=gsub("GeneExt~", "", qpid))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>% 
  select(gene, goid) %>% 
  mutate(goid = paste0("= ", goid))

# Define the output file path
output_file_path <- "C:/Users/suaph281/Desktop/nesh_local/LDT_RNAseq/pannzer2_out/FILTERED_GO_TERMS/spimp_goterms.txt"

# Write the header lines to the file
writeLines(c("(species=Solanum pimpinellifolium) (type=Consortium) (curator=GO)"), con = output_file_path)

# Append the selected data to the file
write.table(out_file, file = output_file_path, append = TRUE, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


#shabro
GO_shabro <- read.csv("shabr/anno.out", header = T, sep = "\t")
GO_shabro_filtered <- GO_shabro %>% 
  mutate(PPV = as.numeric(PPV)) %>% 
  filter(PPV > .4)
GO_shabro_filtered_ids <- GO_shabro_filtered %>% 
  select(qpid) %>% 
  unique()

write_tsv(GO_shabro_filtered, "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/GO/data/shabro_GO_filtered.txt")
write_tsv(GO_shabro_filtered_ids, "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/GO/data/shabro_GO_filtered_IDs.txt")
# for bingo 

out_file <- read.delim("shabr/GO.out") %>%
  filter(ARGOT_PPV >0.4) %>% 
  select(qpid, goid) %>%
  mutate(gene=gsub("GeneExt~", "", qpid))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>% 
  select(gene, goid) %>% 
  mutate(goid = paste0("= ", goid))

# Define the output file path
output_file_path <- "C:/Users/suaph281/Desktop/nesh_local/LDT_RNAseq/pannzer2_out/FILTERED_GO_TERMS/shabro_goterms.txt"

# Write the header lines to the file
writeLines(c("(species=Solanum habrochaites) (type=Consortium) (curator=GO)"), con = output_file_path)

# Append the selected data to the file
write.table(out_file, file = output_file_path, append = TRUE, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#schil
GO_schil <- read.csv("schil/anno.out", header = T, sep = "\t")
GO_schil_filtered <- GO_schil %>% 
  mutate(PPV = as.numeric(PPV)) %>% 
  filter(PPV > .4)
GO_schil_filtered_ids <- GO_schil_filtered %>%
  select(qpid) %>% 
  unique()
  
write_tsv(GO_schil_filtered, "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/GO/data/schil_GO_filtered.txt")
write_tsv(GO_schil_filtered_ids, "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/GO/data/schil_GO_filtered_IDs.txt")


# for bingo 

out_file <- read.delim("schil/GO.out") %>%
  filter(ARGOT_PPV >0.4) %>% 
  select(qpid, goid) %>%
  mutate(gene=gsub("GeneExt~", "", qpid))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>% 
  select(gene, goid) %>% 
  mutate(goid = paste0("= ", goid))

# Define the output file path
output_file_path <- "C:/Users/suaph281/Desktop/nesh_local/LDT_RNAseq/pannzer2_out/FILTERED_GO_TERMS/schil_goterms.txt"

# Write the header lines to the file
writeLines(c("(species=Solanum chilense) (type=Consortium) (curator=GO)"), con = output_file_path)

# Append the selected data to the file
write.table(out_file, file = output_file_path, append = TRUE, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

