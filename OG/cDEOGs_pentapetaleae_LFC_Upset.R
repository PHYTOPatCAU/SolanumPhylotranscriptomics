# Are cDEOGs also DEGs in remaining pentapetaleae?
# 2024_12_14 Severin Einspanier
# Get cDEOG ID, are the gene-IDs also DEGs for the remaining species
# in inf-mock? 
# if yes = 1, else = 0. Then upsetplot
rm(list=ls())
pacman::p_load(tidyverse)
setwd("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/")
cDEOGs_id <- read.csv("OGs/data/coreDEOGs_joined_pentapetals.csv", row.names = 1) %>% 
  select(Orthogroup) %>% 
  unique()
# info: HOGs are npw introduced by orthofinder (hirachial orthogroups), but in the gene count, genes from a different
# Gene.Tree.Parent.Clade are also included. so i accept, that OG_sequence_ids has more rows thant cDEOGs_id
# load degs 

og_sequence_ids <- read.csv("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/orthofinder/NOV_2024_ITAG_PANNZER/Phylogenetic_Hierarchical_Orthogroups/N0.tsv",  
                              header=T, sep="\t", row.names=1, stringsAsFactors = F) %>% 
  #select(!Gene.Tree.Parent.Clade) %>% 
  filter(OG %in% cDEOGs_id$Orthogroup) %>% 
  unique() #%>% 
  

# repeat for all species


# now go through og_sequence_ids and check if the gene-IDs are DEGs in the remaining species
# I will not select a pvalue. we are only interested in: are they expressed?
# using a filter for baseMean > 10

ath_degs <- read.csv2("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/Florent_data/RNAseq_LFC_Infected_vs_Mock_Species_Sucher/Arabidopsis thaliana LFC FRT vs LNI.csv", 
                      header=T, sep=";") %>% 
  filter(baseMean > 10 )

bv_degs <- read.csv("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/Florent_data/RNAseq_LFC_Infected_vs_Mock_Species_Sucher/Beta vulgaris LFC FRT vs LNI.csv", 
                    header=T, sep=";")  %>% 
  filter(baseMean > 10 )

ha_degs <- read.csv("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/Florent_data/RNAseq_LFC_Infected_vs_Mock_Species_Sucher/Helianthus annus LFC FRT vs LNI.csv", 
                    header=T, sep=";")  %>% 
  filter(baseMean > 10 )

pv_degs <- read.csv("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/Florent_data/RNAseq_LFC_Infected_vs_Mock_Species_Sucher/Phaseolus vulgaris LFC FRT vs LNI.csv", 
                    header=T, sep=";")  %>% 
  filter(baseMean > 10 )

rc_degs <- read.csv("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/Florent_data/RNAseq_LFC_Infected_vs_Mock_Species_Sucher/Ricinus communis LFC FRT vs LNI.csv", 
                    header=T, sep=";") %>% 
  mutate(Gene.id = gsub("\\..*", "", Gene.id))  %>% 
  filter(baseMean > 10 )

sl_degs <- read.csv2("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/Florent_data/RNAseq_LFC_Infected_vs_Mock_Species_Sucher/Solanum lycopersicum LFC FRT vs LNI.csv", 
                     header=T, sep=";")  %>% 
  filter(baseMean > 10 )





ath_og_ids <- og_sequence_ids %>% 
  dplyr::select(Arabidopsis_thaliana_TAIR10_prot, OG) %>% 
  unique() %>% 
  separate_rows(Arabidopsis_thaliana_TAIR10_prot, sep = ", ") %>% 
  #filter(Arabidopsis_thaliana_TAIR10_prot != "") %>% 
  mutate(Arabidopsis_thaliana_TAIR10_prot = gsub("\\..*", "", Arabidopsis_thaliana_TAIR10_prot)) %>% 
  mutate(is_deg = if_else(Arabidopsis_thaliana_TAIR10_prot %in% ath_degs$Gene.id, 1, 0)) %>% 
  group_by(OG) %>% 
  summarize(is_deg = sum(is_deg)) %>% 
  mutate(ath = if_else(is_deg > 0, 1, 0)) %>% 
  select(OG, ath)


bv_og_ids <- og_sequence_ids %>% 
  dplyr::select(Beta_vulgaris_BETV122_prot, OG) %>% 
  unique() %>% 
  separate_rows(Beta_vulgaris_BETV122_prot, sep = ", ") %>% 
  mutate(Beta_vulgaris_BETV122_prot = gsub("\\..*", "", Beta_vulgaris_BETV122_prot)) %>% 
  mutate(is_deg = if_else(Beta_vulgaris_BETV122_prot %in% bv_degs$Gene.id, 1, 0)) %>%
  group_by(OG) %>%
  summarize(is_deg = sum(is_deg)) %>% 
  mutate(bv = if_else(is_deg > 0, 1, 0)) %>%
  select(OG, bv)


ha_og_ids <- og_sequence_ids %>%
  dplyr::select(Helianthus_annuus_HanXRQr1_prot, OG) %>% 
  unique() %>% 
  separate_rows(Helianthus_annuus_HanXRQr1_prot, sep = ", ") %>% 
  mutate(Helianthus_annuus_HanXRQr1_prot = gsub("\\..*", "", Helianthus_annuus_HanXRQr1_prot)) %>% 
  mutate(is_deg = if_else(Helianthus_annuus_HanXRQr1_prot %in% ha_degs$Gene.id, 1, 0)) %>%
  group_by(OG) %>%
  summarize(is_deg = sum(is_deg)) %>% 
  mutate(ha = if_else(is_deg > 0, 1, 0)) %>%
  select(OG, ha)



pv_og_ids <- og_sequence_ids %>%
  dplyr::select(Phaseolus_vulgaris_PVULGBAT93_prot, OG) %>% 
  unique() %>% 
  separate_rows(Phaseolus_vulgaris_PVULGBAT93_prot, sep = ", ") %>% 
  mutate(Phaseolus_vulgaris_PVULGBAT93_prot = gsub("\\..*", "", Phaseolus_vulgaris_PVULGBAT93_prot)) %>% 
  mutate(Phaseolus_vulgaris_PVULGBAT93_prot = gsub("T.*", "", Phaseolus_vulgaris_PVULGBAT93_prot)) %>%
  mutate(is_deg = if_else(Phaseolus_vulgaris_PVULGBAT93_prot %in% pv_degs$Gene.id, 1, 0)) %>%
  group_by(OG) %>%
  summarize(is_deg = sum(is_deg)) %>% 
  mutate(pv = if_else(is_deg > 0, 1, 0)) %>%
  select(OG, pv)



rc_og_ids <- og_sequence_ids %>%
  dplyr::select(Ricinus_communis_GCF000151685_prot, OG) %>% 
  unique() %>% 
  separate_rows(Ricinus_communis_GCF000151685_prot, sep = ", ") %>% 
  mutate(Ricinus_communis_GCF000151685_prot = gsub("\\..*", "", Ricinus_communis_GCF000151685_prot)) %>% 
  mutate(Ricinus_communis_GCF000151685_prot = gsub("XP_*", "XM_", Ricinus_communis_GCF000151685_prot)) %>% 
  mutate(is_deg = if_else(Ricinus_communis_GCF000151685_prot %in% rc_degs$Gene.id, 1, 0)) %>%
  group_by(OG) %>%
  summarize(is_deg = sum(is_deg)) %>% 
  mutate(rc = if_else(is_deg > 0, 1, 0)) %>%
  select(OG, rc)



sl_og_ids <- og_sequence_ids %>%
  dplyr::select(Solanum_lycopersicum_ITAG32_prot, OG) %>% 
  unique() %>% 
  separate_rows(Solanum_lycopersicum_ITAG32_prot, sep = ", ") %>% 
  mutate(Solanum_lycopersicum_ITAG32_prot = gsub("\\..*", "", Solanum_lycopersicum_ITAG32_prot)) %>% 
  mutate(is_deg = if_else(Solanum_lycopersicum_ITAG32_prot %in% sl_degs$Gene.id, 1, 0)) %>%
  group_by(OG) %>%
  summarize(is_deg = sum(is_deg)) %>% 
  mutate(sl = if_else(is_deg > 0, 1, 0)) %>%
  select(OG, sl)


merged <- ath_og_ids %>% 
  left_join(bv_og_ids, by = "OG") %>%
  left_join(ha_og_ids, by = "OG") %>%
  left_join(pv_og_ids, by = "OG") %>%
  left_join(rc_og_ids, by = "OG") %>%
  left_join(sl_og_ids, by = "OG") %>% 
  column_to_rownames("OG")
  
# Upset
selected_species=colnames(merged)
library(UpSetR)
svg(paste0("C:/Users/suaph281/Nextcloud/ResiDEvo/sequencing/figures/fig_3/", Sys.Date(),"_UpSet_coreDEOGs_broad.svg"), 
    width=6, height=4, bg="transparent")
upset(merged, 
              nsets = ncol(merged),
              #queries=list(upset_query(intersect=c('schil_comp', 'spen_comp'), color='orange')),
              queries = list(list(query=intersects, params=selected_species, color="green4", active=T)),
              keep.order = F,
              order.by = "freq"
)
dev.off()
