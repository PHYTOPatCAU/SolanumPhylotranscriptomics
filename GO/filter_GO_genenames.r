# Quickly filter anno.out for node-annotation 
rm(list=ls())
library(tidyverse)

GO_TERM <- read.delim("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/GO_terms/spen/anno.out")%>% 
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

write.csv(GO_TERM, "/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/GO_terms/spen/filtered_protein_names_spen.csv")
head(GO_TERM)

GO_TERM <- read.delim("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/GO_terms/slyco/anno.out")%>% 
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

write.csv(GO_TERM, "/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/GO_terms/slyco/filtered_protein_names_slyco.csv")
head(GO_TERM)

GO_TERM <- read.delim("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/GO_terms/spimp/anno.out")%>% 
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

write.csv(GO_TERM, "/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/GO_terms/spimp/filtered_protein_names_spimp.csv")
head(GO_TERM)

GO_TERM <- read.delim("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/GO_terms/shabr/anno.out")%>% 
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

write.csv(GO_TERM, "/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/GO_terms/shabr/filtered_protein_names_shabro.csv")

head(GO_TERM)

GO_TERM <- read.delim("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/GO_terms/schil/anno.out")%>% 
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

write.csv(GO_TERM, "/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/GO_terms/schil/filtered_protein_names_schil.csv")
head(GO_TERM)
