# NON CORE DEOGs
# SPECIES SPECIFIC DEOGs


orthogroups <- read.table("C:/Users/suaph281/Desktop/nesh_local/orthofinder/Results_Sep06/Orthogroups.csv",  header=T, sep="\t", stringsAsFactors = F) 
orthogroups_genid_edited <- orthogroups %>% 
  # select only solanum
  select(!Arabidopsis_thaliana_TAIR10_prot& !Beta_vulgaris_BETV122_prot& !Helianthus_annuus_HanXRQr1_prot& !Phaseolus_vulgaris_PVULGBAT93_prot& !Ricinus_communis_GCF000151685_prot& !Solanum_lycopersicum_ITAG32_prot) %>% 
  summarize(
    # Split the strings and remove "GeneExt~" and other unwanted parts
    schil = str_split(Chilense_longest_orfs_pep, ",\\s*") %>%
      map(~str_replace_all(.x, "^GeneExt~", "") %>%
            str_replace_all(., "\\.p[1-9]", "")%>%    # Remove '.p1' to '.p9'
            str_replace_all("\\.[1-9]", "")%>%     # Remove '.1' to '.8'
            str_replace_all("t\\.", "g.")),         # Replace 't.' with 'g.'),        
    spen = str_split(pennellii_longest_orfs_pep, ",\\s*") %>%
      map(~str_replace_all(.x, "^GeneExt~", "") %>%
            str_replace_all(., "\\.p[1-9]", "")%>%    # Remove '.p1' to '.p9'
            str_replace_all("\\.[1-9]", "")%>%     # Remove '.1' to '.8'
            str_replace_all("t\\.", "g.")),          # Replace 't.' with 'g.'),        
    shabro = str_split(habrochaites_longest_orfs_pep, ",\\s*") %>%
      map(~str_replace_all(.x, "^GeneExt~", "") %>%
            str_replace_all(., "\\.p[1-9]", "")%>%    # Remove '.p1' to '.p9'
            str_replace_all("\\.[1-9]", "") %>%       # Remove '.1' to '.8  ,
            str_replace_all(., "^mRNA_", "")%>%     
            str_replace_all("t\\.", "g.")),          # Replace 't.' with 'g.'    ),
    slyco = str_split(lycopersicoides_longest_orfs_pep, ",\\s*") %>%
      map(~str_replace_all(.x, "^GeneExt~", "") %>%
            str_replace_all(., "\\.p[1-9]", "")%>%    # Remove '.p1' to '.p9'
            str_replace_all("\\.[1-9]", "")%>%     # Remove '.1' to '.8'
            str_replace_all("t\\.", "g.")),          # Replace 't.' with 'g.'),        
    spimp = str_split(pimpinellifolium_longest_orfs_pep, ",\\s*") %>%
      map(~str_replace_all(.x, "^GeneExt~", "") %>%
            str_replace_all(., "\\.p[1-9]", "")%>%    # Remove '.p1' to '.p9'
            str_replace_all("\\.[1-9]", "")%>%     # Remove '.1' to '.8'
            str_replace_all("t\\.", "g.")),          # Replace 't.' with 'g.'),        
    Orthogroup = X
  )
orthogroups_genid_edited_long <- orthogroups_genid_edited %>% 
  rename("S. chilense"=schil, "S. pennellii"=spen, "S. habrochaites"=shabro, "S. lycopersicoides"=slyco, "S. pimpinellifolium"=spimp) %>% 
  pivot_longer(cols=c("S. chilense", "S. pennellii", "S. habrochaites", "S. lycopersicoides", "S. pimpinellifolium"), names_to="species", values_to="genelist") %>% 
  unnest()



write.csv(orthogroups_genid_edited_long, "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/DeSeq/data/2024_09_09_ORTHOGORUPS_name_edited.csv", row.names = F)

