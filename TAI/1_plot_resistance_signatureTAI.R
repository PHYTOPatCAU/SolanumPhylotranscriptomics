# QDR RNAseq Solanum species
# plot signature of the resistance-modules 
# Severin Einspanier

rm(list=ls())
pacman::p_load(myTAI, tidyverse)
setwd("")

# only focus on main module
resistance_module_spen <- "red"
spen_resistance_module <-  read.delim("TOM_INFL_filtered_node_2024-12-17_filtered.tsv") %>% 
  filter(nodeAttr.nodesPresent... == resistance_module_spen) %>%
  mutate(nodeName=tolower(nodeName)) %>% 
  select(nodeName)%>%
  unique()
spen_tai <- read.csv("TAI/data/2024_11_25_spen_TA_reads_collapsed.csv")

PlotCategoryExpr(ExpressionSet = spen_tai,
                 legendName    = "PS",
                 test.stat     = TRUE,
                 type          = "category-centered",
                 distr.type    = "boxplot",
                 gene.set      = spen_resistance_module$nodeName)
# make signature plot
spen_tai_red <- spen_tai %>% 
  filter(GeneID %in% spen_resistance_module$nodeName)

PlotSignature(ExpressionSet = spen_tai_red,
              #gene.set      = spen_resistance_module$nodeName,
              permutations=20000 ,
              measure    = "TAI")

# repeat for S lyco as well

resistance_module_slyco <- "turquoise"
slyco_resistance_module <-  read.delim("TOM_INFL_filtered_node_2024-12-13.tsv") %>% 
  filter(nodeAttr.nodesPresent... == resistance_module_slyco) %>%
  mutate(nodeName=tolower(nodeName)) %>% 
  select(nodeName)%>%
  unique()

slyco_tai <- read.csv("TAI/data/2024_11_25_slyco_TA_reads_collapsed.csv")

PlotCategoryExpr(ExpressionSet = slyco_tai,
                 legendName    = "PS",
                 test.stat     = TRUE,
                 type          = "category-centered",
                 distr.type    = "boxplot",
                 gene.set      = slyco_resistance_module$nodeName)

slyco_tai_turquoise <- slyco_tai %>% 
  filter(GeneID %in% slyco_resistance_module$nodeName)

PlotSignature(ExpressionSet = slyco_tai_turquoise,
              #gene.set      = slyco_resistance_module$nodeName,
              permutations=20000,
              measure    = "TAI")
