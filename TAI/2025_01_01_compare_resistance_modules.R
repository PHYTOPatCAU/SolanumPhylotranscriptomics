# Screen TAI for resistance hubs 

# Severin Einspanier
# 2025_01_01

#Calculate TAI / TDI for each of these gene sets:
  
#   The entire “resistance module” from each species.
#   Overlapping genes (if any) between the two modules.
#   Non-overlapping genes (that might show up in a different module in the second species).

###############################################################
# 0) Setup
###############################################################
rm(list=ls())
pacman::p_load(tidyverse, myTAI)

setwd("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/")

###############################################################
# 1) Read and clean input data
###############################################################

## 1.1) S. pennellii
spen_genid2goid <- read.csv("OGs/data/genid2goid_spen.csv", row.names = 1) %>% 
  filter(GeneID != "") %>% 
  mutate(gene = GeneID %>%
           gsub("GeneExt~", "", .) %>%
           gsub("mRNA_", "", .) %>%
           gsub("t\\.peak", "g.peak", .) %>%
           gsub("t\\.minus", "g.minus", .) %>%
           gsub("t\\.plus", "g.plus", .) %>%
           gsub("\\.[1-9].*|\\.p[1-9].*", "", .)) %>% 
  select(OG, gene) %>% 
  distinct()

spen_nodes <- read.delim("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/WGCNA/spen/TOM_INFL_filtered_node_2024-12-17_filtered.tsv") %>% 
  left_join(spen_genid2goid, by = c("nodeName" = "gene")) %>% 
  distinct()

## 1.2) S. lycopersicoides
slyco_genid2goid <- read.csv("OGs/data/genid2goid_slyco.csv", row.names = 1) %>% 
  filter(GeneID != "") %>% 
  mutate(gene = GeneID %>%
           gsub("GeneExt~", "", .) %>%
           gsub("mRNA_", "", .) %>%
           gsub("t\\.peak", "g.peak", .) %>%
           gsub("t\\.minus", "g.minus", .) %>%
           gsub("t\\.plus", "g.plus", .) %>%
           gsub("\\.[1-9].*|\\.p[1-9].*", "", .)) %>% 
  select(OG, gene) %>% 
  distinct()

slyco_nodes <- read.delim("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/WGCNA/slyco/TOM_INFL_filtered_node_2024-12-19_filtered.tsv") %>% 
  left_join(slyco_genid2goid, by = c("nodeName" = "gene")) %>% 
  distinct()


## Get TAI data 
spen_tai <- read.csv("TAI/data/2024_11_25_spen_TA_reads_collapsed.csv")
slyco_tai <- read.csv("TAI/data/2024_11_25_slyco_TA_reads_collapsed.csv")

###############################################################
# 2) Get TAI of the whole resistance-modules
###############################################################

# 2.1) Identify spen "red" OGs
resistance_module_spen <- c("red", "pink", "blue")

spen_resistance_module <- spen_nodes %>%
  filter(nodeAttr.nodesPresent... %in% resistance_module_spen) %>% 
  select(OG, nodeName)%>%
  unique()%>% 
  drop_na(OG)

# 2.2) Identify slyco "turquoise" OGs

resistance_module_slyco <- c("turquoise", "pink","green")
slyco_resistance_module <- slyco_nodes %>%
  filter(nodeAttr.nodesPresent... %in% resistance_module_slyco) %>% 
  select(OG, nodeName)%>%
  unique() %>% 
  drop_na(OG)

red_expression <- spen_tai %>%
  filter(GeneID %in% tolower(spen_resistance_module$nodeName))

TAI(red_expression)

#    LA1303mock LA1303sclero   LA2963mock LA2963sclero 
# 1.765792     1.775637     1.770824     1.780546 

turquoise_expression <- slyco_tai %>%
  filter(GeneID %in% tolower(slyco_resistance_module$nodeName))

TAI(turquoise_expression)
#   LA2777mock LA2777sclero   LA2951mock LA2951sclero 
# 1.978264     2.005911     1.958683     1.982362 

###############################################################
# 3) Get TAI of the overlapping genes
###############################################################

# 3.1) Identify overlapping genes
overlapping_genes <- spen_resistance_module %>% 
  rename("spen_GENEID"=nodeName) %>% 
  left_join(slyco_resistance_module, by = "OG") %>% 
  rename("slyco_GENEID"=nodeName) %>% 
  drop_na(spen_GENEID, slyco_GENEID)

spen_overlap_TAI <- spen_tai %>%
  filter(GeneID %in% tolower(overlapping_genes$spen_GENEID))

TAI(spen_overlap_TAI)

#     LA1303mock LA1303sclero   LA2963mock LA2963sclero 
# 1.681193     1.695101     1.687719     1.699765 

slyco_overlap_TAI <- slyco_tai %>%
  filter(GeneID %in% tolower(overlapping_genes$slyco_GENEID))

TAI(slyco_overlap_TAI)

#    LA2777mock LA2777sclero   LA2951mock LA2951sclero 
# 1.733932     1.747360     1.724876     1.741134 

###############################################################
# 4) Get TAI of the non-overlapping genes
###############################################################

# 4.1) Identify non-overlapping genes
spen_non_overlapping_genes <- spen_resistance_module %>% 
  rename("spen_GENEID"=nodeName) %>% 
  anti_join(overlapping_genes, by = "spen_GENEID") %>% 
  select(OG, spen_GENEID) 

spen_non_overlap_TAI <- spen_tai %>%
  filter(GeneID %in% tolower(spen_non_overlapping_genes$spen_GENEID))  

TAI(spen_non_overlap_TAI)
#     LA1303mock LA1303sclero   LA2963mock LA2963sclero 
# 1.915807     1.924948     1.917388     1.931397 

slyco_non_overlapping_genes <- slyco_resistance_module %>% 
  rename("slyco_GENEID"=nodeName) %>% 
  anti_join(overlapping_genes, by = "slyco_GENEID") %>% 
  select(OG, slyco_GENEID)

slyco_non_overlap_TAI <- slyco_tai %>%
  filter(GeneID %in% tolower(slyco_non_overlapping_genes$slyco_GENEID))  
TAI(slyco_non_overlap_TAI)
#     LA2777mock LA2777sclero   LA2951mock LA2951sclero 
# 2.401127     2.453106     2.367097     2.404728

# plot contribution 

PlotEnrichment( ExpressionSet = spen_tai,
                test.set = spen_non_overlap_TAI$GeneID,
                  legendName    = "PS",
                measure = "log-foldchange")

PlotCategoryExpr(ExpressionSet = spen_tai,
                 legendName    = "PS",
                 test.stat     = TRUE,
                 type          = "stage-centered",
                 distr.type    = "boxplot",
                 log.expr      = TRUE,
                 gene.set      = spen_non_overlap_TAI$GeneID)
