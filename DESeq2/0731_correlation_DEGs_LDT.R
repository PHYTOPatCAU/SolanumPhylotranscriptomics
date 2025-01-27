## 
## Get Corelation Gene Expression after infection and infection score/LDT

rm(list=ls())
library(tidyverse)

DEGs_susceptible <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/DeSeq/data/DeSeq_OUT/combined_SUS_inf_mock_DEGs.csv")
DEGs_resistant <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/DeSeq/data/DeSeq_OUT/combined_RES_inf_mock_DEGs.csv")
stats_ldt <- read.csv("C:/Users/suaph281/Nextcloud/papers/2024_Phenotyping_Solanum_Sclerotinia/tables/2023_10_24_ldt_Populations_Means.csv", 
                      header = T, row.names = 1) 


# schil

data_schil <- read.csv("Z:/00_navautron/EXPERIMENTS/2024_07_25_schil_1980/0/slopes_2.txt", header=F, sep="\t")
head(data_schil)
colnames(data_schil) <- c("ID", "lag", "slope", "comm", "LDT")

id <- data.frame(
  ID=seq(1:79), 
  genotype=c(rep("LA3111-3", 29), rep("LA4107-10", 39), rep("C32", 11))
)


data_schil_edited <- data_schil %>% 
  left_join(id, by="ID") %>% 
  filter(comm=="ok" & genotype %in% c("LA3111-3", "LA4107-10")) %>% 
  mutate(start_date="2024_07_25", 
         inoculum="ss", 
         species="S. chilense", 
         SampleID=paste0("2024_07_25_schil_", genotype, ID, inoculum)) %>% 
  group_by(genotype) %>% 
  summarize(lsmean=mean(LDT), 
            genotype=genotype) %>% 
  unique()

stats_ldt <- stats_ldt %>% 
  bind_rows(data_schil_edited)

accessions <- data.frame(
  species=c("S. chilense", "S. chilense", "S. pennellii","S. pennellii", "S. habrochaites", "S. habrochaites", "S. pimpinellifolium", "S. pimpinellifolium", "S. lycopersicoides",  "S. lycopersicoides"),
  accession=c("LA3111-3", "LA4107-10", "LA1303", "LA2963", "LA2167", "LA1721", "LA2347", "LA1332", "LA2777", "LA2951"),
  level=c("res", "sus", "sus", "res", "sus", "res", "sus", "res", "sus", "res")
)


LDT_sus <- accessions %>% 
  filter(level=="sus") %>%
  left_join(stats_ldt, by = c("accession" = "genotype"))

LDT_res <- accessions %>% 
  filter(level=="res") %>%
  left_join(stats_ldt, by = c("accession" = "genotype"))


DEGs_susceptible_edit <- DEGs_susceptible %>% 
  left_join(LDT_sus, by = c("species" = "species")) %>% 
  group_by(species) %>%
  mutate(UP=ifelse(log2FoldChange>0,1,0)) %>%
  summarise(n_DEGs = n(),
            up=sum(UP),
            down=n()-sum(UP),
            lsmean = lsmean, 
            level=level,
            genotype=accession) %>% 
  unique() 

DEGs_all <- DEGs_resistant %>% 
  left_join(LDT_res, by = c("species" = "species")) %>% 
  group_by(species) %>%
  mutate(UP=ifelse(log2FoldChange>0,1,0)) %>%
  summarise(n_DEGs = n(),
            up=sum(UP),
            down=n()-sum(UP),
            lsmean = lsmean, 
            level=level,
            genotype=accession) %>% 
  unique() %>% 
  bind_rows(DEGs_susceptible_edit)

svg(paste0("C:/Users/suaph281/Nextcloud/ResiDEvo/sequencing/figures/fig_2/",
           Sys.Date(),"_corr_LFC_LDT.svg"), 
    width = 6, height = 4)

# Load necessary libraries
library(ggplot2)
library(ggpubr) # For stat_cor()

# Your plot with correlation analysis
p1 <- ggplot(DEGs_all, aes(x = lsmean, y = n_DEGs)) +
  geom_smooth(method = "lm", se = TRUE, linetype = "solid", size = 1, alpha = 0.3) +  
  stat_cor(method = "pearson", label.x = 6.5, label.y = 12000, 
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), size = 5) + 
  geom_point(aes(col = species),size = 3) +
  theme_bw() +
  labs(x = "lsmean(LDT)",
       y = "# DEGs \n(INF vs. MOCK)", 
       col = "Species") +
  scale_x_continuous(
    limit = c(5.5, 8.5), 
    breaks = c(5.88, 6.291569, 6.984716, 7.677864, 8.371011),
    labels = c(6, 9, 18, 36, 72)
  ) +
  theme(
    axis.text = element_text(size = 11, color = "Black"),
    axis.title = element_text(size = 11, color = "Black"),
    legend.title = element_text(size = 13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.text = element_text(size = 11), 
    legend.spacing.x = unit(.5, 'cm'),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(size = .4, color = "grey50", linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    #axis.text.x = element_text(hjust = 0.5, face = "italic", color = "black"),
    axis.title.x = element_blank()
  )

# Print the plot
png(paste0("C:/Users/suaph281/Nextcloud/ResiDEvo/sequencing/figures/fig_2/",
           Sys.Date(),"_corr_LFC_LDT.png"), 
    width = 7, height = 4, units = "in", res = 600)
svg(paste0("C:/Users/suaph281/Nextcloud/ResiDEvo/sequencing/figures/fig_2/",
           Sys.Date(),"_corr_LFC_LDT.svg"), 
    width = 7, height = 4)
print(p1)
dev.off()

ggplot(DEGs_all, aes(x=lsmean, y=up, col=species))+
  geom_point()
ggplot(DEGs_all, aes(x=lsmean, y=down, col=species))+
  geom_point()



# Cor with schil




ldt_wth_schil <- ldt_phenotypes %>% 
  filter(genotype%in%accessions$accession) %>%
  filter(comm=="ok") %>% 
  bind_rows(data_schil_edited)



accessions_order <- c("LA4107-10", "LA3111-3", "LA2167", "LA1721",
                      "LA2777", "LA2951", "LA1303", "LA2963", "LA2347", "LA1332")



