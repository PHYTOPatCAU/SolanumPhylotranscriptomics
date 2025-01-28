# Fisher-testing of resistance module preservation spen and slyco
# Severin Einspanier
# 2025_01_01
# 2025_01_17 adjusted for more defence-modules

###############################################################
# 0) Setup
###############################################################
rm(list=ls())
pacman::p_load(tidyverse)

setwd("")

###############################################################
# 1) Read and clean input data
###############################################################

## 1.1) S. pennellii
spen_genid2goid <- read.csv("OGs/data/genid2goid_spen.csv", row.names = 1) %>% 
  filter(GeneID != "") %>% 
  mutate(gene = GeneID %>%
           gsub("GeneExt~", "", .) %>%
           gsub("mRNA_", "", .) %>%
           gsub("t\\.peak", "g.peak", .) %>%
           gsub("t\\.minus", "g.minus", .) %>%
           gsub("t\\.plus", "g.plus", .) %>%
           gsub("\\.[1-9].*|\\.p[1-9].*", "", .)) %>% 
  select(OG, gene) %>% 
  distinct()

spen_nodes <- read.delim("spen/TOM_INFL_filtered_node_2024-12-17_filtered.tsv") %>% 
  left_join(spen_genid2goid, by = c("nodeName" = "gene")) %>% 
  distinct()

## 1.2) S. lycopersicoides
slyco_genid2goid <- read.csv("OGs/data/genid2goid_slyco.csv", row.names = 1) %>% 
  filter(GeneID != "") %>% 
  mutate(gene = GeneID %>%
           gsub("GeneExt~", "", .) %>%
           gsub("mRNA_", "", .) %>%
           gsub("t\\.peak", "g.peak", .) %>%
           gsub("t\\.minus", "g.minus", .) %>%
           gsub("t\\.plus", "g.plus", .) %>%
           gsub("\\.[1-9].*|\\.p[1-9].*", "", .)) %>% 
  select(OG, gene) %>% 
  distinct()

slyco_nodes <- read.delim("slyco/TOM_INFL_filtered_node_2024-12-19_filtered.tsv") %>% 
  left_join(slyco_genid2goid, by = c("nodeName" = "gene")) %>% 
  distinct()

###############################################################
# 2) FORWARD CONTRAST:
#    Are spen "red" OGs enriched in slyco modules?
###############################################################

# 2.1) Identify spen "red" OGs, as well as pink & blue 
resistance_module_spen <- c("red", "blue", "pink")
spen_resistance_module <- spen_nodes %>%
  filter(nodeAttr.nodesPresent... %in% resistance_module_spen)

spen_red_OGs <- spen_resistance_module$OG %>%
  unique()

# 2.2) Define slyco background
slyco_all_OGs <- slyco_nodes$OG %>%
  unique()

# M = # of spen "red" OGs that appear in slyco
M <- length(intersect(spen_red_OGs, slyco_all_OGs))
# N = total # of unique OGs in slyco
N <- length(slyco_all_OGs)

# 2.3) Fisher’s test for each slyco module
slyco_modules <- slyco_nodes %>%
  filter(!is.na(nodeAttr.nodesPresent...)) %>%
  distinct(nodeAttr.nodesPresent...) %>%
  pull(nodeAttr.nodesPresent...)

fisher_results <- data.frame()

for (mod in slyco_modules) {
  # OGs in this slyco module
  slyco_mod_OGs <- slyco_nodes %>%
    filter(nodeAttr.nodesPresent... == mod) %>%
    pull(OG) %>%
    unique()
  
  # k = # of OGs in slyco module
  k <- length(slyco_mod_OGs)
  
  # x = overlap with spen "red"
  x <- length(intersect(spen_red_OGs, slyco_mod_OGs))
  
  # Contingency table
  contingency <- matrix(c(x,             # in spen_red & slyco module
                          k - x,         # not spen_red but in slyco module
                          M - x,         # spen_red but not in this slyco module
                          N - k - (M - x)),
                        nrow=2, byrow=TRUE)
  
  ft <- fisher.test(contingency, alternative="greater")
  
  fisher_results <- rbind(
    fisher_results,
    data.frame(
      slyco_module      = mod,
      p_value           = ft$p.value,
      odds_ratio        = ft$estimate,
      overlap           = x,
      slyco_module_size = k,
      spen_red_size     = M,
      background        = N,
      stringsAsFactors  = FALSE
    )
  )
}

# 2.4) Adjust p-values & plot
fisher_results <- fisher_results %>%
  mutate(
    padj         = p.adjust(p_value, method = "fdr"),
    neg_log10_p  = -log10(padj),
    significance = if_else(padj < 0.05, "significant", "not_significant")
  )

(spen_red <- ggplot(fisher_results, 
       aes(x = slyco_module, 
           y = neg_log10_p, 
           size = overlap, 
           color = significance)) +
  geom_point(alpha = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_discrete(guide = "none")+
  labs(
    x      = "S. lycopersicoides Module",
    y      = "-log10(p-adjust.)",
    color  = "Significance",
    size   = "Overlap Count"#,
    #title  = "Enrichment of 'red' S. pennellii genes \nin S. lycopersicoides modules"
  ) +    
     scale_size_area(max_size = 10, 
                     breaks=c(250,500,1000,1500,2000))+
  theme_bw() +
  theme(axis.text.x      = element_text(angle=45, hjust=1, color="black"),
        axis.text        = element_text(color="black"),
        strip.text       = element_text(size=12, color="black", face="italic"),
        axis.title       = element_text(size=12, color="black"),
        text             = element_text(family="serif"),
        strip.background = element_rect(fill="white"))
)

###############################################################
# 3) REVERSE CONTRAST:
#    Are slyco "magenta" OGs enriched in spen modules?
###############################################################

# 3.1) slyco "magenta" OGs
resistance_module_slyco <- c("turquoise", "green", "pink")
slyco_resistance_module <- slyco_nodes %>%
  filter(nodeAttr.nodesPresent... %in% resistance_module_slyco)

slyco_magenta_OGs <- slyco_resistance_module$OG %>%
  unique()

# 3.2) spen background
spen_all_OGs <- spen_nodes$OG %>%
  unique()

# M = # of slyco "magenta" OGs appearing in spen
M <- length(intersect(slyco_magenta_OGs, spen_all_OGs))
# N = total # of unique OGs in spen
N <- length(spen_all_OGs)

# 3.3) Fisher’s test for each spen module
spen_modules <- spen_nodes %>%
  filter(!is.na(nodeAttr.nodesPresent...)) %>%
  distinct(nodeAttr.nodesPresent...) %>%
  pull(nodeAttr.nodesPresent...)

fisher_results_rev <- data.frame()

for (mod in spen_modules) {
  spen_mod_OGs <- spen_nodes %>%
    filter(nodeAttr.nodesPresent... == mod) %>%
    pull(OG) %>%
    unique()
  
  k <- length(spen_mod_OGs)  # # of OGs in this spen module
  x <- length(intersect(slyco_magenta_OGs, spen_mod_OGs))  # overlap
  
  contingency <- matrix(c(x,
                          k - x,
                          M - x,
                          N - k - (M - x)),
                        nrow=2, byrow=TRUE)
  
  ft <- fisher.test(contingency, alternative="greater")
  
  fisher_results_rev <- rbind(
    fisher_results_rev,
    data.frame(
      spen_module       = mod,
      p_value           = ft$p.value,
      odds_ratio        = ft$estimate,
      overlap           = x,
      spen_module_size  = k,
      slyco_magenta_size= M,
      background        = N,
      stringsAsFactors  = FALSE
    )
  )
}

# 3.4) Adjust p-values & plot
fisher_results_rev <- fisher_results_rev %>%
  mutate(
    padj         = p.adjust(p_value, method = "fdr"),
    neg_log10_p  = -log10(padj),
    significance = if_else(padj < 0.05, "significant", "not_significant")
  )

(slyco_magenta <- ggplot(fisher_results_rev, 
       aes(x = spen_module, 
           y = neg_log10_p, 
           size = overlap, 
           color = significance)) +
  geom_point(alpha = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_discrete(guide = "none")+
    #ylim(-10, 170)+
  labs(
    x      = "S. pennellii Module",
    y      = "-log10(p-adjust.)",
    color  = "Significance",
    size   = "Overlap Count"#,
    #title  = "Enrichment of S.lycopersicoides \n'turquoise' in S. pennellii modules"
  ) + 
    scale_size_area(max_size = 10, 
                    breaks=c(250,500,1000,1500,2000)) +
  theme_bw() +
  theme(axis.text.x      = element_text(angle=45, hjust=1, color="black"),
        axis.text        = element_text(color="black"),
        strip.text       = element_text(size=12, color="black", face="italic"),
        axis.title       = element_text(size=12, color="black"),
        text             = element_text(family="serif"),
        strip.background = element_rect(fill="white"))
)


######
# Viz both 
######

svg(paste0("figures/fig_6/", Sys.Date(),
           "_resistance_genes_fisher.svg"), width = 6.496063, height = 3.14)
(p_tot <- ggpubr::ggarrange(spen_red, slyco_magenta,
                            ncol=2, nrow=1, common.legend = T, 
                            labels = c("A)", "B)")))
dev.off()

ggsave("figures/fig_6/resistance_genes_fisher_test.png", 
       p_tot, width = 16.5, height = 8, units = "cm", dpi = 900,
       bg="white")
