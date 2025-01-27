# Are sDEOGs induced?

rm(list=ls())


library(tidyverse)
# RES

setwd("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/")

DEGs_res <- read.csv("DeSeq/data/DeSeq_OUT/combined_RES_inf_mock_DEGs.csv") %>% 
  select(log2FoldChange, padj, species, GeneID) %>% 
  rename(Res_log2FoldChange = log2FoldChange, Res_padj = padj)
DEGs_sus <- read.csv("DeSeq/data/DeSeq_OUT/combined_SUS_inf_mock_DEGs.csv") %>% 
  select(log2FoldChange, padj, species, GeneID) %>% 
  rename(Sus_log2FoldChange = log2FoldChange, Sus_padj = padj)

sDEOGs <- read.csv("OGs/data/OG_sDEOGs_IDs.csv", check.names = F)

# if GeneID and species match in DEGs_res set isDEOGs-value = 1, else=0

isDEOGs <- sDEOGs %>% 
  dplyr::group_by(species) %>% 
  left_join(DEGs_res, by = c("GeneID" = "GeneID", "species")) %>% 
  left_join(DEGs_sus, by = c("GeneID" = "GeneID", "species")) %>% 
  select(OG, species, GeneID, Res_log2FoldChange, Res_padj, Sus_log2FoldChange, Sus_padj) %>%
  mutate(isDEOGs_res = ifelse(!is.na(Res_log2FoldChange), 1, 0),
         isDEOGs_sus = ifelse(!is.na(Sus_log2FoldChange), 1, 0)) %>%
  mutate(sum=isDEOGs_res+isDEOGs_sus)

summary(isDEOGs$sum)

plot_isDEOGs <- isDEOGs %>% 
  select(OG, species, sum) %>%
  group_by(OG, species) %>% 
  # select one gene per OG. Always the max.
  filter(sum == max(sum)) %>%
  mutate(sum=ifelse(sum>0, 1, 0)) %>%
  group_by(OG, species) %>% 
  unique() %>% 
  group_by(species) %>%
  mutate(tot_sum=sum(sum), 
         n=n(),)

svg(paste0("C:/Users/suaph281/Nextcloud/ResiDEvo/sequencing/figures/fig_3/", Sys.Date(), "_induction_sDEOGs.svg"), 
    width=4, height=5, bg="transparent")

(p1 <- plot_isDEOGs %>% 
    summarize( species=species, share=sum(sum)/n(), 
               tot_sum=tot_sum) %>% 
    unique() %>% 
    ggplot(aes(x=fct_reorder(as.factor(species), tot_sum), 
               y=share*100, fill=species))+
    geom_col(show.legend = F)+
    labs(x="Species", 
         y="Share of infection-induced \nsDEOGs", 
         guides="Species")+
    ylim(0,100)+
    scale_fill_manual(values=(c("S. habrochaites"="#ffffff", 
                               "S. chilense" = "#c6c9d8",
                               "S. pimpinellifolium"="#7b86c4",
                               "S. pennellii" = "#3043af", 
                               "S. lycopersicoides"="#0f1e70"))) + # Scale of blues
    theme_bw()+
    #scale_fill_manual(values = c("LA1282" = "#007F94", "LA1809" = "#EED78D", "LA1941" = "#C22B26")) +
    #scale_y_continuous(limits=c(0,6.4), breaks = c(2,4,6))+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(vjust=1, face="italic",
                                     size=11, color="black",
                                     angle=45, hjust=1),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm'),
          panel.grid.major.x =  element_blank()
  )  
)
dev.off()


# output for GO terms 

go_isDEOGs <- isDEOGs %>% 
  select(OG, species, sum) %>%
  group_by(OG, species) %>% 
  # select one gene per OG. Always the max.
  filter(sum == max(sum)) %>%
  mutate(sum=ifelse(sum>0, 1, 0)) %>%
  group_by(OG, species) %>% 
  unique() %>% 
  filter(sum>0)

go_isDEOGs %>% 
  filter(species=="S. pimpinellifolium") %>% 
  ungroup() %>% 
  select(OG) %>% 
  clipr::write_clip()

go_nisDEOGs <- isDEOGs %>% 
  select(OG, species, sum) %>%
  group_by(OG, species) %>% 
  # select one gene per OG. Always the max.
  filter(sum == max(sum)) %>%
  mutate(sum=ifelse(sum>0, 1, 0)) %>%
  group_by(OG, species) %>% 
  unique() %>% 
  filter(sum==0)

go_nisDEOGs %>% 
  filter(species=="S. pimpinellifolium") %>% 
  ungroup() %>%
  select(OG) %>% 
  clipr::write_clip()
