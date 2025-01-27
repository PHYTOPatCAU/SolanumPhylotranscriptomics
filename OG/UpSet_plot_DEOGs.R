# Plot cDEOGs vs sDEOGs

# plotting
library(tidyverse)
library(UpSetR)
setwd("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/")

orthogroups_genid_wide <- read.csv("OGs/data/OG_with_DEGs_binary_wide.csv") %>% 
  column_to_rownames('X') %>% 
  rename('S. chilense'='S..chilense','S. pennellii'='S..pennellii', 'S. pimpinellifolium'='S..pimpinellifolium', 
         'S. habrochaites'='S..habrochaites', 'S. lycopersicoides'='S..lycopersicoides')


selected_species <- colnames(orthogroups_genid_wide)[1:(ncol(orthogroups_genid_wide))] 
selected_species

#orthogroups_filtered_all_upset <- orthogroups_filtered_all[1:5]

svg(paste0("C:/Users/suaph281/Nextcloud/ResiDEvo/sequencing/figures/fig_3/", Sys.Date(), "_UpSet_DEOGs.svg"), 
    width=8, height=4, bg="transparent")
UpSetR::upset(orthogroups_genid_wide, 
              nsets = ncol(orthogroups_genid_wide),
              queries = list(list(query=intersects, params=selected_species, color="tomato", active=T)),
              keep.order = F,
              order.by = "freq"
)
dev.off()



# Check this subset with other Pentapetalea

core_DEOGs <- orthogroups_genid_wide %>%
  mutate(sum=rowSums(.)) %>%
  filter(sum==5) %>% 
  rownames_to_column("Orthogroup") %>% 
  select(Orthogroup)

orthogroup_counts <- read.csv("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/orthofinder/NOV_2024_ITAG_PANNZER/Orthogroups/Orthogroups.GeneCount.tsv",  header=T, sep="\t", stringsAsFactors = F) 
colnames(orthogroup_counts)[1]="Orthogroup"

core_DEOGs_joined <- core_DEOGs %>% 
  left_join(orthogroup_counts, by="Orthogroup") %>% 
  mutate(`A. thaliana` = if_else(Arabidopsis_thaliana_TAIR10_prot > 0, 1,0), 
         `B. vulgaris`= if_else(Beta_vulgaris_BETV122_prot>0,1,0),
         `H. annuus` = if_else(Helianthus_annuus_HanXRQr1_prot>0,1,0),
         `P. vulgaris`= if_else(Phaseolus_vulgaris_PVULGBAT93_prot>0,1,0),
         `R. communis`= if_else(Ricinus_communis_GCF000151685_prot>0,1,0),
         `S. lycopersicum`=if_else(Solanum_lycopersicum_ITAG32_prot>0,1,0)) %>% 
  select(`A. thaliana`,`B. vulgaris`, `H. annuus`, `P. vulgaris`, `R. communis`,`S. lycopersicum`, Orthogroup)

write.csv(core_DEOGs_joined, "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/OGs/data/coreDEOGs_joined_pentapetals.csv")


# and plot as well:

selected_species <- colnames(core_DEOGs_joined)[1:(ncol(core_DEOGs_joined))-1] 
selected_species

svg(paste0("C:/Users/suaph281/Nextcloud/ResiDEvo/sequencing/figures/fig_3/", Sys.Date(),"_UpSet_coreDEOGs_broad.svg"), 
    width=6, height=4, bg="transparent")
upset(core_DEOGs_joined, 
      nsets = ncol(core_DEOGs_joined),
      #queries=list(upset_query(intersect=c('schil_comp', 'spen_comp'), color='orange')),
      queries = list(list(query=intersects, params=selected_species, color="green4", active=T)),
      keep.order = F,
      order.by = "freq"
)
dev.off()
