##
## Merge the OG-GOs of the five species 
## Approach: 
## 1. Infer gene ID and OG ID
## 2. Read the GO files of the five species
## 3. Filter the GO terms based on PPV > .4
## 4. Merge the GO terms of the five species and keep the shared GO terms (? might get complicated)
## 5. Write the merged GO terms to a file



rm(list=ls())
library(tidyverse)


##### 1. Infer gene ID and OG ID #####

sng_cp_ogs <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/OGs/data/orthogroups_singlecopy_goid.csv")

# S. lycopersicoides

# I will adjust the geneID later, 'cause then I can directly join with GO terms

OG_data_slyco <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/OGs/data/genid2goid_slyco.csv") %>% 
  filter(OG %in%sng_cp_ogs$OG) %>% 
  filter(GeneID !="") %>% 
  mutate(gene=gsub("GeneExt~", "", GeneID))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))  


#OG_code_slyco <- OG_data_slyco %>% 
#  mutate(Orthogroup=str_replace_all(V2, "/mnt/c/Users/suaph281/Desktop/nesh_local/LDT_RNAseq/orthofinder/single_copy_OG/sequences//", "") %>% 
#           str_replace_all(., ".fasta", ""), 
#         species="S. lycopersicoides") %>%
#  rename(GeneID=V1) %>%
#  dplyr::select( !V2)
#rm(OG_data_slyco)         
           
# S. pennellii
OG_data_spen <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/OGs/data/genid2goid_spen.csv") %>% 
  filter(OG %in%sng_cp_ogs$OG) %>% 
  filter(GeneID !="") %>% 
  mutate(gene=gsub("GeneExt~", "", GeneID))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))  

#OG_code_spen <- OG_data_spen %>% 
#  mutate(Orthogroup=str_replace_all(V2, "/mnt/c/Users/suaph281/Desktop/nesh_local/LDT_RNAseq/orthofinder/single_copy_OG/sequences//", "") %>% 
#           str_replace_all(., ".fasta", ""), 
#         species="S. pennellii") %>%
#  rename(GeneID=V1) %>%
#  dplyr::select( !V2)
#rm(OG_data_spen)


# S. habrochaites

OG_data_shabro <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/OGs/data/genid2goid_shabro.csv") %>% 
  filter(OG %in%sng_cp_ogs$OG) %>% 
  filter(GeneID !="") %>% 
  mutate(gene=gsub("GeneExt~", "", GeneID))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))  

#OG_code_shab <- OG_data_shab %>% 
#  mutate(Orthogroup=str_replace_all(V2, "/mnt/c/Users/suaph281/Desktop/nesh_local/LDT_RNAseq/orthofinder/single_copy_OG/sequences//", "") %>% 
#           str_replace_all(., ".fasta", ""), 
#         species="S. habrochaites") %>%
#  rename(GeneID=V1) %>%
#  dplyr::select( !V2)
#rm(OG_data_shab)

# S. pimpinellifolium

OG_data_spimp <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/OGs/data/genid2goid_spimp.csv") %>% 
  filter(OG %in%sng_cp_ogs$OG) %>% 
  filter(GeneID !="") %>% 
  mutate(gene=gsub("GeneExt~", "", GeneID))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))  

#OG_code_spimp <- OG_data_spimp %>% 
#  mutate(Orthogroup=str_replace_all(V2, "/mnt/c/Users/suaph281/Desktop/nesh_local/LDT_RNAseq/orthofinder/single_copy_OG/sequences//", "") %>% 
#           str_replace_all(., ".fasta", ""), 
#         species="S. pimpinellifolium") %>%
#  rename(GeneID=V1) %>%
#  dplyr::select( !V2)
#rm(OG_data_spimp)

# S. chilense

OG_data_schil <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/OGs/data/genid2goid_schil.csv") %>% 
  filter(OG %in%sng_cp_ogs$OG) %>% 
  filter(GeneID !="") %>% 
  mutate(gene=gsub("GeneExt~", "", GeneID))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))  

#OG_code_schi <- OG_data_schi %>% 
#  mutate(Orthogroup=str_replace_all(V2, "/mnt/c/Users/suaph281/Desktop/nesh_local/LDT_RNAseq/orthofinder/single_copy_OG/sequences//", "") %>% 
#           str_replace_all(., ".fasta", ""), 
#         species="S. chilense") %>%
#  rename(GeneID=V1) %>%
#  dplyr::select( !V2)
#rm(OG_data_schi)



#### 2, 3. Read the GO files of the five species and filter the GO terms based on PPV > .4 ####

# 2.1 S. lycopersicoides
GOterms_slyco <- read.table("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/pannzer2_out/manual_pannzer2/slyco/GO.out",  
                            header=T, sep="\t", 
                            colClasses = c("character", "character", "character", "character", "numeric", "numeric", "numeric", "numeric")) %>% 
  filter(ARGOT_PPV >0.4) %>% 
  select(qpid, goid) %>%
  mutate(gene=gsub("GeneExt~", "", qpid))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>% 
  right_join(OG_data_slyco, by=c("gene"="gene") )

summary(is.na(GOterms_slyco$OG))
# no NA --> all GO-terms have an OG ID

# 2.2 S. pennellii  
GOterms_spen <- read.table("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/pannzer2_out/manual_pannzer2/spen/GO.out", 
                           header=T, sep="\t", 
                           colClasses = c("character", "character", "character", "character", "numeric", "numeric", "numeric", "numeric")) %>% 
  filter(ARGOT_PPV >0.4) %>% 
  select(qpid, goid) %>%
  mutate(gene=gsub("GeneExt~", "", qpid))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>% 
  right_join(OG_data_spen, by=c("gene"="gene") )


summary(is.na(GOterms_spen$OG))

# 2.3 S. habrochaites

GOterms_shabro <- read.table("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/pannzer2_out/manual_pannzer2/shabr/GO.out", 
                           header=T, sep="\t", 
                           colClasses = c("character", "character", "character", "character", "numeric", "numeric", "numeric", "numeric")) %>% 
  filter(ARGOT_PPV >0.4) %>% 
  select(qpid, goid) %>%
  mutate(gene=gsub("GeneExt~", "", qpid))%>%  
  mutate(gene=gsub("mRNA:", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>% 
  right_join(OG_data_shabro, by=c("gene"="gene") )

summary(is.na(GOterms_shabro$OG))

# 2.4 S. pimpinellifolium

GOterms_spimp <- read.table("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/pannzer2_out/manual_pannzer2/spimp/GO.out", 
                            header=T, sep="\t", colClasses = c("character", "character", "character", "character", "numeric", "numeric", "numeric", "numeric")) %>% 
  filter(ARGOT_PPV >0.4) %>% 
  select(qpid, goid) %>%
  mutate(gene=gsub("GeneExt~", "", qpid))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>% 
  right_join(OG_data_spimp, by=c("gene"="gene") )

summary(is.na(GOterms_spimp$OG))


# 2.5 S. chilense

GOterms_schil <- read.table("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/pannzer2_out/manual_pannzer2/schil/GO.out", 
                           header=T, sep="\t", 
                           colClasses = c("character", "character", "character", "character", "numeric", "numeric", "numeric", "numeric")) %>% 
  filter(ARGOT_PPV >0.4) %>% 
  select(qpid, goid) %>%
  mutate(gene=gsub("GeneExt~", "", qpid))%>%  
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>% 
  right_join(OG_data_schil, by=c("gene"="gene") )


summary(is.na(GOterms_schil$OG))


#### 4. Merge the GO terms of the five species and keep the shared GO terms ####
# and filter for PPV > .4, and GO term must be assigned to each OG at least n times

# per OG
# get the goid's which are represented at least twice. 

GOterms_all <- GOterms_slyco %>% 
  bind_rows(GOterms_spen) %>% 
  bind_rows(GOterms_shabro) %>% 
  bind_rows(GOterms_spimp) %>% 
  bind_rows(GOterms_schil) %>% 
  dplyr::select(OG, goid) %>% 
  group_by(OG, goid) %>%
  drop_na(goid) %>% 
  summarise(count=n(),
            goid=goid
            ) %>% 
  unique() %>% 
  filter(count > 2)





# How many GO terms has each OG? 

(count_GO_per_OG <- GOterms_all %>% 
  group_by(OG) %>% 
  summarise(n_GO=n()) %>% 
  arrange(desc(n_GO)) %>% 
  ggplot(aes(x=n_GO))+
  geom_histogram(binwidth=1)
)

# is this normal to have so many GO terms per OG sometimes?
# This might need some filtering here. 

 
### 5. Write the merged GO terms to a file ####


# Select the required columns and modify the goid column
out_file <- GOterms_all %>% 
  select(OG, goid) %>%
  mutate(goid = paste0("= ", goid))

# Define the output file path
output_file_path <- "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/OGs/data/OG_GO_terms.txt"

# Write the header lines to the file
writeLines(c("(species=SolanumMix) (type=Consortium) (curator=GO)"), con = output_file_path)

# Append the selected data to the file
write.table(out_file, file = output_file_path, append = TRUE, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
