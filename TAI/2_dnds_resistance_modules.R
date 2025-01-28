# QDR RNAseq Solanum species
# Extract dnds of interesting resistance modules
# Severin Einspanier


###############################################################
# Setup
###############################################################
pacman::p_load(myTAI, tidyverse)

 setwd("")
# (Uncomment/modify if you need a working directory set)

###############################################################
# 1) Read and clean dNdS data
###############################################################
slyco_dnds <- read.csv("TDI/data/2024_10_11_dnds_slyco_melo.csv") %>% 
  mutate(
    query_id = query_id %>%
      gsub("GeneExt~|\\.1|\\.2|\\.3|\\.p[1-9]|mRNA:", "", .) %>%
      tolower()
  )

spen_dnds <- read.csv("TDI/data/2024_10_11_dnds_spen_melo.csv") %>% 
  mutate(
    query_id = query_id %>%
      gsub("GeneExt~|\\.1|\\.2|\\.3|\\.p[1-9]|mRNA:", "", .) %>%
      tolower()
  )

###############################################################
# 2) Read and clean input data for spen and slyco
###############################################################

## 2.1) S. pennellii
spen_genid2goid <- read.csv("OGs/data/genid2goid_spen.csv", row.names = 1) %>% 
  filter(GeneID != "") %>%
  mutate(
    gene = GeneID %>%
      gsub("GeneExt~", "", .) %>%
      gsub("mRNA_", "", .) %>%
      gsub("t\\.peak", "g.peak", .) %>%
      gsub("t\\.minus", "g.minus", .) %>%
      gsub("t\\.plus", "g.plus", .) %>%
      gsub("\\.[1-9].*|\\.p[1-9].*", "", .)
  ) %>% 
  select(OG, gene) %>% 
  distinct()

spen_nodes <- read.delim("TOM_INFL_filtered_node_2024-12-17_filtered.tsv") %>% 
  left_join(spen_genid2goid, by = c("nodeName" = "gene")) %>% 
  distinct()

## 2.2) S. lycopersicoides
slyco_genid2goid <- read.csv("OGs/data/genid2goid_slyco.csv", row.names = 1) %>% 
  filter(GeneID != "") %>%
  mutate(
    gene = GeneID %>%
      gsub("GeneExt~", "", .) %>%
      gsub("mRNA_", "", .) %>%
      gsub("t\\.peak", "g.peak", .) %>%
      gsub("t\\.minus", "g.minus", .) %>%
      gsub("t\\.plus", "g.plus", .) %>%
      gsub("\\.[1-9].*|\\.p[1-9].*", "", .)
  ) %>% 
  select(OG, gene) %>% 
  distinct()

slyco_nodes <- read.delim("TOM_INFL_filtered_node_2024-12-19_filtered.tsv") %>% 
  left_join(slyco_genid2goid, by = c("nodeName" = "gene")) %>% 
  distinct()

###############################################################
# 3) Identify "resistance modules" in spen and slyco
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

###############################################################
# 4) dNdS of the entire modules
###############################################################
whole_red <- spen_dnds %>%
  filter(query_id %in% tolower(spen_resistance_module$nodeName)) %>% 
  mutate(level = "whole_module") %>%
  select(query_id, level, dNdS) %>%
  distinct()

whole_turquoise <- slyco_dnds %>% 
  filter(query_id %in% tolower(slyco_resistance_module$nodeName)) %>%
  mutate(level = "whole_module") %>%
  select(query_id, level, dNdS) %>%
  distinct()

###############################################################
# 5) dNdS of overlapping genes
###############################################################
overlapping_genes <- spen_resistance_module %>% 
  rename(spen_GENEID = nodeName) %>% 
  merge(slyco_resistance_module, by = "OG") %>% 
  rename(slyco_GENEID = nodeName)

spen_overlap_dnds <- spen_dnds %>%
  filter(query_id %in% tolower(overlapping_genes$spen_GENEID)) %>%
  mutate(level = "overlap") %>%
  select(query_id, level, dNdS) %>%
  distinct()

slyco_overlap_dnds <- slyco_dnds %>%
  filter(query_id %in% tolower(overlapping_genes$slyco_GENEID)) %>%
  mutate(level = "overlap") %>%
  select(query_id, level, dNdS) %>%
  distinct()

###############################################################
# 6) dNdS of non-overlapping genes
###############################################################
spen_non_overlapping_genes <- spen_resistance_module %>% 
  rename(spen_GENEID = nodeName) %>% 
  anti_join(overlapping_genes, by = "spen_GENEID") %>%
  select(OG, spen_GENEID)

spen_non_overlap_dnds <- spen_dnds %>%
  filter(query_id %in% tolower(spen_non_overlapping_genes$spen_GENEID)) %>%
  mutate(level = "unique") %>%
  select(query_id, level, dNdS) %>%
  distinct()

slyco_non_overlapping_genes <- slyco_resistance_module %>% 
  rename(slyco_GENEID = nodeName) %>% 
  anti_join(overlapping_genes, by = "slyco_GENEID") %>%
  select(OG, slyco_GENEID)

slyco_non_overlap_dnds <- slyco_dnds %>%
  filter(query_id %in% tolower(slyco_non_overlapping_genes$slyco_GENEID)) %>%
  mutate(level = "unique") %>%
  select(query_id, level, dNdS) %>%
  distinct()

###############################################################
# 7) Merge data frames & plot
###############################################################
# 7.1) S. lycopersicoides
merged_df_slyco <- bind_rows(whole_turquoise, slyco_overlap_dnds, slyco_non_overlap_dnds) %>%
  filter(is.finite(dNdS))

p_slyco <- ggplot(merged_df_slyco, aes(x = level, y = dNdS, fill = level)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  ylim(0, 3.2) +
  labs(
    x    = "Level",
    y    = "dNdS",
    fill = "Level"
    # color = "Level" # optional if you also want a color legend
  ) +
  theme_bw() +
  theme(
    axis.text        = element_text(size = 12, color = "black"),
    strip.text       = element_text(size = 15, color = "black", face = "italic"),
    axis.title       = element_text(size = 15, color = "black"),
    strip.background = element_rect(fill = "white")
  )

# 7.2) S. pennellii
merged_df_spen <- bind_rows(whole_red, spen_overlap_dnds, spen_non_overlap_dnds) %>%
  filter(is.finite(dNdS))

p_spen <- ggplot(merged_df_spen, aes(x = level, y = dNdS, fill = level)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  ylim(0, 3.2) +
  labs(
    x    = "Level",
    y    = "dNdS",
    fill = "Level"
  ) +
  theme_bw() +
  theme(
    axis.text        = element_text(size = 12, color = "black"),
    strip.text       = element_text(size = 15, color = "black", face = "italic"),
    axis.title       = element_text(size = 15, color = "black"),
    strip.background = element_rect(fill = "white")
  )

###############################################################
# 8) Combine and save
###############################################################
p_tot <- ggpubr::ggarrange(
  p_spen, p_slyco,
  ncol = 2, nrow = 1, common.legend = TRUE,
  labels = c("S. pennellii", "S. lycopersicoides")
)

ggsave(
  filename = "2025_01_02_dnds_resistance_modules.png",
  plot     = p_tot,
  width    = 10,
  height   = 5,
  dpi      = 900,
  bg       = "white"
)

# Optional: If you want to add significance testing on the plot,
# uncomment the ggpubr::stat_compare_means() lines and adjust as needed.
