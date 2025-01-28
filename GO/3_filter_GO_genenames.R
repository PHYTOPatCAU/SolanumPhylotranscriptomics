# QDR RNAseq Solanum species
# Filter anno.out for node-annotation (GO-IDs must match GENE names)
# Severin Einspanier


rm(list=ls())
library(tidyverse)

GO_TERM <- read.delim("spen/anno.out")%>% 
    mutate(gene=gsub(">", "", qpid)) %>%
    mutate(gene=gsub("GeneExt~", "", gene))%>% 
    mutate(gene=gsub("mRNA_", "", gene))%>%  
    mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
    mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
    mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
    mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>%
    filter(PPV < .4 & type =="DE") %>%
    dplyr::select(gene, desc)
dim(GO_TERM)

write.csv(GO_TERM, "spen/filtered_protein_names_spen.csv")
head(GO_TERM)

GO_TERM <- read.delim("slyco/anno.out")%>% 
    mutate(gene=gsub(">", "", qpid)) %>%
    mutate(gene=gsub("GeneExt~", "", gene))%>% 
    mutate(gene=gsub("mRNA_", "", gene))%>%  
    mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
    mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
    mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
    mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>%
    filter(PPV < .4 & type =="DE") %>%
    dplyr::select(gene, desc)
dim(GO_TERM)

write.csv(GO_TERM, "slyco/filtered_protein_names_slyco.csv")
head(GO_TERM)

GO_TERM <- read.delim("spimp/anno.out")%>% 
    mutate(gene=gsub(">", "", qpid)) %>%
    mutate(gene=gsub("GeneExt~", "", gene))%>% 
    mutate(gene=gsub("mRNA_", "", gene))%>%  
    mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
    mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
    mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
    mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>%
    filter(PPV < .4 & type =="DE") %>%
    dplyr::select(gene, desc)

dim(GO_TERM)

write.csv(GO_TERM, "spimp/filtered_protein_names_spimp.csv")
head(GO_TERM)

GO_TERM <- read.delim("shabr/anno.out")%>% 
    mutate(gene=gsub(">", "", qpid)) %>%
    mutate(gene=gsub("GeneExt~", "", gene))%>% 
    mutate(gene=gsub("mRNA_", "", gene))%>%  
    mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
    mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
    mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
    mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>%
    filter(PPV < .4 & type =="DE") %>%
    dplyr::select(gene, desc)
dim(GO_TERM)

write.csv(GO_TERM, "shabr/filtered_protein_names_shabro.csv")

head(GO_TERM)

GO_TERM <- read.delim("schil/anno.out")%>% 
    mutate(gene=gsub(">", "", qpid)) %>%
    mutate(gene=gsub("GeneExt~", "", gene))%>% 
    mutate(gene=gsub("mRNA_", "", gene))%>%  
    mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
    mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
    mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
    mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>%
    filter(PPV < .4 & type =="DE") %>%
    dplyr::select(gene, desc)
dim(GO_TERM)

write.csv(GO_TERM, "schil/filtered_protein_names_schil.csv")
head(GO_TERM)
