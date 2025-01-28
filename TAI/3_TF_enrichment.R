# QDR RNAseq Solanum species
# plot signature for the of TF-GRNs
# Severin Einspanier
rm(list=ls())
pacman::p_load(tidyverse, myTAI, ggpubr)
setwd("")

set_one <- read.table("Sopen09g001470_grn_network_nodes_filtered.txt",
                      header=T, sep="\t")

# load meta data for annotation 
datTraits <- read.csv2("DeSeq/data/sample_infos.csv") %>%
  filter(species == "spen") %>%
  mutate(pseudofactor=paste0(genotype, "_", treatment, "_", rep)) %>% 
  select(X, pseudofactor)

spen_tai <- read.csv("TAI/data/2024_11_25_spen_TA_reads.csv") %>% 
  pivot_longer(cols = 3:18, names_to = "X", values_to = "exp") %>% 
  left_join(datTraits, by = c("X" = "X")) %>% 
  select(!X) %>% 
  pivot_wider(names_from = pseudofactor, values_from = exp) 

myTAI::PlotEnrichment(
  ExpressionSet=spen_tai, 
  test.set=set_one$genes,
  complete.bg = T
  
)

# Set Two

set_two <- read.table("Sopen10g006210_grn_network_nodes_filtered.txt",
                      header=T, sep="\t")

myTAI::PlotEnrichment(
  ExpressionSet=spen_tai, 
  test.set=set_two$genes,
  complete.bg = T
  
)

PlotSignature(ExpressionSet=(spen_tai %>% filter(GeneID%in%tolower(set_two$genes)))
              )

# Set Three

set_three <- read.table("Sopen05g003630_grn_network_nodes_filtered.txt",
                      header=T, sep="\t")

myTAI::PlotEnrichment(
  ExpressionSet=spen_tai, 
  test.set=set_three$genes,
  complete.bg = T
  
)

PlotSignature(ExpressionSet=(spen_tai %>% filter(GeneID%in%tolower(set_three$genes)))
)


#########
## Repeat for Slyco
#########

datTraits <- read.csv2("DeSeq/data/sample_infos.csv") %>%
  filter(species == "slyco") %>%
  mutate(pseudofactor=paste0(genotype, "_", treatment, "_", rep)) %>% 
  select(X, pseudofactor)

spen_tai <- read.csv("TAI/data/2024_11_25_slyco_TA_reads.csv") %>% 
  pivot_longer(cols = 3:18, names_to = "X", values_to = "exp") %>% 
  left_join(datTraits, by = c("X" = "X")) %>% 
  select(!X) %>% 
  pivot_wider(names_from = pseudofactor, values_from = exp) 


set_one <- read.table("Solyd02g064330_grn_network_nodes_filtered.txt",
                      header=T, sep="\t")

myTAI::PlotEnrichment(
  ExpressionSet=spen_tai, 
  test.set=set_one$genes,
  complete.bg = T
  
)


set_two <- read.table("Solyd09g071070_grn_network_nodes_filtered.txt",
                      header=T, sep="\t")

myTAI::PlotEnrichment(
  ExpressionSet=spen_tai, 
  test.set=set_two$genes,
  complete.bg = T
  
)

set_three <- read.table("Solyd03g050610_grn_network_nodes_filtered.txt",
                      header=T, sep="\t")

myTAI::PlotEnrichment(
  ExpressionSet=spen_tai, 
  test.set=set_three$genes,
  complete.bg = T
  
)
