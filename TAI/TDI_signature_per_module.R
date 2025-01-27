# Plot the per-module TDI signature of S. pennellii and S. lycopersicoides
# Severin Einspanier
# 2025_01_05
rm(list=ls())
pacman::p_load(tidyverse, myTAI, ggpubr)
setwd("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/")

# spen

modules <- c("brown", "yellow", "blue", "turquoise", "black", "green", "red", "pink")

# load meta data for annotation 
datTraits <- read.csv2("DeSeq/data/sample_infos.csv") %>%
  filter(species == "spen") %>%
  mutate(pseudofactor=paste0(genotype, "_", treatment, "_", rep)) %>% 
  select(X, pseudofactor)

spen_TDI <- read.csv("TDI/data/2024_11_25_spen_DS_reads.csv") %>% 
  pivot_longer(cols = 3:18, names_to = "X", values_to = "exp") %>% 
  left_join(datTraits, by = c("X" = "X")) %>% 
  select(!X) %>% 
  pivot_wider(names_from = pseudofactor, values_from = exp)

spen_resistance_module <- read.delim("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/WGCNA/spen/TOM_INFL_filtered_node_2024-12-17_filtered.tsv") %>% 
  mutate(nodeName = tolower(nodeName)) %>%
  rename(module_name=nodeAttr.nodesPresent...) %>% 
  unique()

plots <- list()

for (module in modules) {
  module_resistance <- spen_resistance_module %>% 
    filter(module_name == module)
  
  spen_TDI_filtered <- spen_TDI %>% 
    filter(GeneID %in% module_resistance$nodeName)
  
  plot <- PlotSignature(ExpressionSet = spen_TDI_filtered,
                        permutations = 20000,
                        measure = "TDI", xlab="Sample",
                        y.ticks = 6)
  
  plots[[module]] <- plot
}
svg("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/TEST_IMG_TDI_SPEN_all.svg", 
    width = 10, height = 20)
ggarrange(plotlist = plots, ncol = 2, nrow = 4, 
          labels=modules)
dev.off()



# Repeat with slyco

modules <- c("magenta", "brown", "yellow", "blue", "turquoise", "black", "green", "red", "pink")


datTraits <- read.csv2("DeSeq/data/sample_infos.csv") %>%
  filter(species == "slyco") %>%
  mutate(pseudofactor=paste0(genotype, "_", treatment, "_", rep)) %>% 
  select(X, pseudofactor)

slyco_TDI <- read.csv("TDI/data/2024_11_25_slyco_DS_reads.csv") %>% 
  pivot_longer(cols = 3:18, names_to = "X", values_to = "exp") %>% 
  left_join(datTraits, by = c("X" = "X")) %>% 
  select(!X) %>% 
  pivot_wider(names_from = pseudofactor, values_from = exp)


slyco_resistance_module <- read.delim("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/WGCNA/slyco/TOM_INFL_filtered_node_2024-12-13.tsv") %>% 
  mutate(nodeName = tolower(nodeName)) %>%
  rename(module_name=nodeAttr.nodesPresent...) %>% 
  unique()

plots <- list()

for (module in modules) {
  module_resistance <- slyco_resistance_module %>% 
    filter(module_name == module)
  
  slyco_TDI_filtered <- slyco_TDI %>% 
    filter(GeneID %in% module_resistance$nodeName)
  
  plot <- PlotSignature(ExpressionSet = slyco_TDI_filtered,
                        permutations = 20000,
                        measure = "TDI", 
                        xlab="Sample",
                        y.ticks=6)
  
  plots[[module]] <- plot
}
svg("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/TEST_IMG_TDI_SLYCO_all.svg", 
    width = 10, height = 20)
ggarrange(plotlist = plots, ncol = 2, nrow = 5, 
          labels=modules)
dev.off()
