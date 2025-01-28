# QDR RNAseq Solanum species
# Get Corelation Gene Expression after infection and infection score/LDT
# Severin Einspanier 
# 2024_10_13

# Clear the workspace
rm(list=ls())

# Load necessary libraries
library(tidyverse)

# Read in the data for susceptible and resistant DEGs
DEGs_susceptible <- read.csv("DeSeq_OUT/combined_SUS_inf_mock_DEGs.csv")
DEGs_resistant <- read.csv("DeSeq_OUT/combined_RES_inf_mock_DEGs.csv")

# Read in the LDT statistics
stats_ldt <- read.csv("2023_10_24_ldt_Populations_Means.csv", 
                      header = T, row.names = 1) 

# Read in the schil data
data_schil <- read.csv("slopes_2.txt", header=F, sep="\t")
head(data_schil)
colnames(data_schil) <- c("ID", "lag", "slope", "comm", "LDT")

# Create a dataframe for genotype IDs
id <- data.frame(
  ID=seq(1:79), 
  genotype=c(rep("LA3111-3", 29), rep("LA4107-10", 39), rep("C32", 11))
)

# Edit the schil data
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

# Combine the schil data with the LDT statistics
stats_ldt <- stats_ldt %>% 
  bind_rows(data_schil_edited)

# Create a dataframe for accessions
accessions <- data.frame(
  species=c("S. chilense", "S. chilense", "S. pennellii","S. pennellii", "S. habrochaites", "S. habrochaites", "S. pimpinellifolium", "S. pimpinellifolium", "S. lycopersicoides",  "S. lycopersicoides"),
  accession=c("LA3111-3", "LA4107-10", "LA1303", "LA2963", "LA2167", "LA1721", "LA2347", "LA1332", "LA2777", "LA2951"),
  level=c("res", "sus", "sus", "res", "sus", "res", "sus", "res", "sus", "res")
)

# Filter and join the LDT data for susceptible and resistant accessions
LDT_sus <- accessions %>% 
  filter(level=="sus") %>%
  left_join(stats_ldt, by = c("accession" = "genotype"))

LDT_res <- accessions %>% 
  filter(level=="res") %>%
  left_join(stats_ldt, by = c("accession" = "genotype"))

# Edit the susceptible DEGs data
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

# Combine the resistant DEGs data with the edited susceptible DEGs data
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

# Create the plot
svg(paste0("figures/fig_2/",
           Sys.Date(),"_corr_LFC_LDT.svg"), 
    width = 6, height = 4)

# Load necessary libraries for plotting
library(ggplot2)
library(ggpubr) # For stat_cor()

# Plot with correlation analysis
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

# Print the plot to PNG and SVG files
png(paste0("figures/fig_2/",
           Sys.Date(),"_corr_LFC_LDT.png"), 
    width = 7, height = 4, units = "in", res = 600)
svg(paste0("figures/fig_2/",
           Sys.Date(),"_corr_LFC_LDT.svg"), 
    width = 7, height = 4)
print(p1)
dev.off()