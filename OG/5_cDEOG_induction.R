# QDR RNAseq Solanum species
# This script is used to assess how many of the cDEOGs are induced upon infection 
# -> icDEOGs for both resistant (res) and susceptible (sus) genes
# Severin Einspanier

# Load required libraries
rm(list=ls())
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)

# Set working directory
setwd("")

# Load data
DEGs_res <- read.csv("DeSeq/data/DeSeq_OUT/combined_RES_inf_mock_DEGs.csv")
DEGs_sus <- read.csv("DeSeq/data/DeSeq_OUT/combined_SUS_inf_mock_DEGs.csv")
core_DEOGs <- read.csv("OGs/data/OG_cDEOGs_IDs.csv", check.names = FALSE)

# Identify induced cDEOGs for resistant (res) and susceptible (sus) genes
induced_cDEOGs <- core_DEOGs %>% 
  group_by(species) %>% 
  mutate(
    induced_sus = ifelse(GeneID %in% na.omit(DEGs_sus[['GeneID']]), 1, 0),
    induced_res = ifelse(GeneID %in% na.omit(DEGs_res[['GeneID']]), 1, 0)
  ) %>% 
  select(OG, species, induced_sus, induced_res)

# Create UpSet plots for susceptible (sus) data
induced_cDEOGs_sus <- induced_cDEOGs %>% 
  select(OG, species, induced_sus) %>% 
  pivot_wider(names_from = species, values_from = induced_sus) %>% 
  column_to_rownames('OG')

selected_species_sus <- colnames(induced_cDEOGs_sus)
upset(
  induced_cDEOGs_sus, 
  nsets = ncol(induced_cDEOGs_sus), 
  sets = rev(selected_species_sus), 
  keep.order = TRUE,
  order.by = "freq",
  mb.ratio = c(0.6, 0.4),
  set_size.show = TRUE
)

# Create UpSet plots for resistant (res) data
induced_cDEOGs_res <- induced_cDEOGs %>% 
  select(OG, species, induced_res) %>% 
  pivot_wider(names_from = species, values_from = induced_res) %>% 
  column_to_rownames('OG')

selected_species_res <- colnames(induced_cDEOGs_res)
upset(
  induced_cDEOGs_res, 
  nsets = ncol(induced_cDEOGs_res), 
  sets = rev(selected_species_res), 
  keep.order = TRUE,
  order.by = "freq",
  mb.ratio = c(0.6, 0.4),
  set_size.show = TRUE
)

# Create a combined UpSet plot for both res and sus data
# An OG is considered induced if it is a DEG in either sus or res
combined_induced_cDEOGs <- core_DEOGs %>% 
  group_by(species) %>% 
  mutate(
    induced_sus = ifelse(GeneID %in% na.omit(DEGs_sus[['GeneID']]), 1, 0),
    induced_res = ifelse(GeneID %in% na.omit(DEGs_res[['GeneID']]), 1, 0)
  ) %>% 
  select(OG, species, induced_sus, induced_res) %>% 
  mutate(icDEOGs = ifelse(induced_sus + induced_res > 0, 1, 0)) %>% 
  select(OG, species, icDEOGs) %>% 
  pivot_wider(names_from = species, values_from = icDEOGs) %>% 
  column_to_rownames('OG')

selected_species_combined <- colnames(combined_induced_cDEOGs)

# Save the combined plot as an SVG
svg(paste0("figures/fig_3/", Sys.Date(), "_cDEOGs_induction.svg"),
    width = 8, height = 5)
upset(
  combined_induced_cDEOGs, 
  nsets = ncol(combined_induced_cDEOGs), 
  sets = rev(selected_species_combined), 
  keep.order = TRUE,
  order.by = "freq",
  mb.ratio = c(0.6, 0.4),
  set_size.show = TRUE
)
dev.off()

induced_cDEOGs <- core_DEOGs %>% 
  dplyr::group_by(species) %>% 
  mutate(induced_sus = ifelse(GeneID %in% na.omit(DEGs_sus[['GeneID']]), 1, 0), 
         induced_res = ifelse(GeneID %in% na.omit(DEGs_res[['GeneID']]), 1, 0)) %>% 
  select(OG, species, induced_sus, induced_res) %>% 
  mutate(icDEOGs=ifelse(induced_sus + induced_res > 0, 1, 0)) %>% 
  select(OG, species, icDEOGs) %>% 
  pivot_wider(names_from = species, values_from = icDEOGs) %>% 
  column_to_rownames('OG')


# get table 

tab <- induced_cDEOGs %>% 
  mutate(sum = rowSums(.)) %>% 
  group_by(sum) %>%
  summarise(n = n())
write.csv(tab, "data/cDEOGs_induction_summary.csv", row.names=FALSE)  

# Extract GO terms of the cicDEOGs
cicDEOGs <- induced_cDEOGs %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum == 5) %>% 
  mutate(OG = rownames(.)) %>% 
  select(OG) %>% 
  rownames(NULL) %>% 
  clipr::write_clip()

induced_cDEOGs %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum > 3) %>% 
  mutate(OG = rownames(.)) %>% 
  select(OG) %>% 
  rownames(NULL) %>% 
  clipr::write_clip()

# write ID of cicDEOGs 

cicDEOGs <- induced_cDEOGs %>% 
  mutate(sum = rowSums(.)) %>% 
  filter(sum > 4) %>% 
  mutate(OG = rownames(.)) %>% 
  select(OG) %>% 
  rownames(NULL)


write.csv(cicDEOGs, "data/cicDEOGs_IDs_5.csv", row.names=FALSE)
# OLD!
# Create UpSet plots using ComplexHeatmap
m1 <- make_comb_mat(induced_cDEOGs_res, min_set_size = 0)
m1norm <- normalize_comb_mat(m1, full_comb_sets = TRUE)
m2 <- make_comb_mat(induced_cDEOGs_sus, min_set_size = 0)
m2norm <- normalize_comb_mat(m2, full_comb_sets = TRUE)

png("sequencing/figures/cDEOGs_induction_res.png",
    width = 10, height = 4, units = 'in', res = 600)
ComplexHeatmap::UpSet(m1, comb_order = order(comb_size(m1))) 
dev.off()

png("sequencing/figures/cDEOGs_induction_sus.png",
    width = 10, height = 4, units = 'in', res = 600)
ComplexHeatmap::UpSet(m2, comb_order = order(comb_size(m2)))
dev.off()