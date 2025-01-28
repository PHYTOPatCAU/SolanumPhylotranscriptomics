# QDR RNAseq Solanum species
# Draw heatmap of cicDEOGs expression pattern 
# Severin Einspanier

#1. get cicDEOGs
rm(list=ls())
library(tidyverse)
library(reshape2)
library(pheatmap)
setwd("")

DEGs_res <- read.csv("DeSeq/data/DeSeq_OUT/combined_RES_inf_mock_DEGs.csv")
DEGs_sus <- read.csv("DeSeq/data/DeSeq_OUT/combined_SUS_inf_mock_DEGs.csv")

core_DEOGs <- read.csv("OGs/data/OG_cDEOGs_IDs.csv", check.names = F)

induced_cDEOGs <- core_DEOGs %>% 
  dplyr::group_by(species) %>% 
  mutate(induced_sus = ifelse(GeneID %in% na.omit(DEGs_sus[['GeneID']]), 1, 0),
         induced_res = ifelse(GeneID %in% na.omit(DEGs_res[['GeneID']]), 1, 0)) %>% 
  select(OG, species, induced_sus, induced_res) %>% 
  mutate(icDEOGs=ifelse(induced_sus + induced_res > 0, 1, 0))  %>% 
  select(OG, species, icDEOGs) %>% 
  pivot_wider(names_from = species, values_from = icDEOGs) %>% 
  column_to_rownames('OG')

cicDEGO_id<- induced_cDEOGs %>%
  #select(!`S. habrochaites`) %>% 
  mutate(sum=rowSums(.)) %>%
  filter(sum>4) %>% 
  mutate(OG=rownames(.)) %>%
  select(OG) %>% 
  clipr::write_clip(.)


# get expression of those 

#expr <- read.csv("DeSeq/data/old/combined_ALL_inf_mock_LFCs.csv") 
expr <- read.csv("DeSeq/data/DeSeq_OUT/combined_INF_res_sus_DEGs.csv")

schil <- read.csv("OGs/data/genid2goid_schil.csv") %>% 
  select(OG, GeneID) %>% 
  filter(GeneID != "" & OG %in% cicDEGO_id$OG) %>% 
  left_join(expr, by=c("GeneID"="GeneID")) %>% 
  drop_na(log2FoldChange) %>% 
  filter(species=="S. chilense") %>%
  group_by(OG) %>%
  #slice_min(order_by = abs(padj), n = 1) %>% 
  slice_max(order_by = abs(log2FoldChange), n = 1) %>% 
  #reframe(log2FoldChange=mean(log2FoldChange)) %>%
  unique() %>% 
  rename( "S. chilense"=log2FoldChange)

shabro <- read.csv("OGs/data/genid2goid_shabro.csv")%>% 
  select(OG, GeneID) %>% 
  filter(GeneID != ""& OG %in% cicDEGO_id$OG) %>% 
  left_join(expr, by=c("GeneID"="GeneID")) %>% 
  drop_na(log2FoldChange) %>% 
  group_by(OG) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1) %>% 
  #reframe(log2FoldChange=mean(log2FoldChange)) %>%
  unique() %>% 
  rename( "S. habrochaites"=log2FoldChange)

slyco <- read.csv("OGs/data/genid2goid_slyco.csv") %>% 
  select(OG, GeneID) %>% 
  filter(GeneID != ""& OG %in% cicDEGO_id$OG)%>% 
  left_join(expr, by=c("GeneID"="GeneID")
  ) %>% 
  drop_na(log2FoldChange) %>% 
  group_by(OG) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1) %>% 
  #reframe(log2FoldChange=mean(log2FoldChange)) %>%
  unique() %>% 
  rename( "S. lycopersicoides"=log2FoldChange)

spen <- read.csv("OGs/data/genid2goid_spen.csv") %>% 
  select(OG, GeneID) %>% 
  filter(GeneID != ""&OG %in% cicDEGO_id$OG)%>% 
  left_join(expr, by=c("GeneID"="GeneID")) %>% 
  drop_na(log2FoldChange) %>% 
  group_by(OG) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1) %>% 
  #reframe(log2FoldChange=mean(log2FoldChange)) %>%
  unique() %>% 
  rename( "S. pennellii"=log2FoldChange)

spimp <- read.csv("OGs/data/genid2goid_spimp.csv") %>% 
  select(OG, GeneID) %>% 
  filter(GeneID !="" &  OG %in% cicDEGO_id$OG)%>% 
  left_join(expr, by=c("GeneID"="GeneID")) %>% 
  drop_na(log2FoldChange) %>% 
  group_by(OG) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1) %>% 
  #reframe(log2FoldChange=mean(log2FoldChange)) %>%
  unique() %>% 
  rename( "S. pimpinellifolium"=log2FoldChange)

merged_cicDEOG_expression <- schil %>% 
  left_join(shabro, by="OG") %>%
  left_join(slyco, by="OG") %>%
  left_join(spen, by="OG") %>%
  left_join(spimp, by="OG") %>% 
  select(OG, "S. pimpinellifolium", "S. pennellii", "S. lycopersicoides", "S. habrochaites", "S. chilense") %>% 
  column_to_rownames(var = "OG") %>%
  t() %>% 
  as.matrix()

names <- lapply(rownames(merged_cicDEOG_expression), 
                function(x) bquote(italic(.(x))))

# Draw the heatmap
svg("figures/fig_4/cicDEOGs_expression.svg", 
    width = 12, height = 5)
png("figures/fig_4/cicDEOGs_expression.png", 
    width = 16.5/2, height = 8, unit="cm", res=960)
pheatmap::pheatmap(merged_cicDEOG_expression,
                   scale = "none",
         cluster_rows = T, 
         cluster_cols = T, 
         treeheight_row=10,
         treeheight_col=10,
         #color = colorRampPalette(c("blue", "white", "red"))(50),
         #main = "ciDEOGs log2FoldChange(RES-SUS)",
         labels_row =as.expression(names), 
         border_color = "black",
         fontsize_row = 7,
         fontsize_col = 6)

dev.off()
