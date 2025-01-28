# QDR RNAseq Solanum species
# Plot cDEOGs vs sDEOGs as upset plot
# Severin Einspanier

library(tidyverse)
library(UpSetR)
setwd("")

orthogroups_genid_wide <- read.csv("OGs/data/OG_with_DEGs_binary_wide.csv") %>% 
  column_to_rownames('X') %>% 
  rename('S. chilense'='S..chilense','S. pennellii'='S..pennellii', 'S. pimpinellifolium'='S..pimpinellifolium', 
         'S. habrochaites'='S..habrochaites', 'S. lycopersicoides'='S..lycopersicoides')


selected_species <- colnames(orthogroups_genid_wide)[1:(ncol(orthogroups_genid_wide))] 
selected_species

svg(paste0("figures/fig_3/", Sys.Date(), "_UpSet_DEOGs.svg"), 
    width=8, height=4, bg="transparent")
UpSetR::upset(orthogroups_genid_wide, 
              nsets = ncol(orthogroups_genid_wide),
              queries = list(list(query=intersects, params=selected_species, color="tomato", active=T)),
              keep.order = F,
              order.by = "freq"
)
dev.off()
