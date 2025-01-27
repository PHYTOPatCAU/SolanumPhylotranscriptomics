# Module correlations 
# but with a more complex model
# Severin Einspanier
# as of 2024_11_29
# S. pimpinellifolium
rm(list=ls())
pacman::p_load(WGCNA, lme4, tidyverse, pheatmap, car, effectsize)

options(stringsAsFactors = FALSE)
setwd("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/")

ids <- read.delim("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/PROTEOME/spimp_curated_proteome_OG_pannzer_dedub_ids.txt",
                  header=FALSE) %>%
  mutate(gene=gsub(">", "", V1),
         gene=gsub("GeneExt~|mRNA_", "", gene),
         gene=gsub("t\\.peak", "g.peak", gene),
         gene=gsub("t\\.minus", "g.minus", gene),
         gene=gsub("t\\.plus", "g.plus", gene),
         gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene)) %>%
  select(gene)

# Filter the gene expression:
datExpr <- read.csv("DeSeq/data/norm_counts_all_rlog.csv") %>%
  filter(species=="S. pimpinellifolium" & Geneid %in% ids$gene) %>%
  select(!species & !source) %>%
  pivot_wider(names_from = Geneid, values_from = normalized_count) %>%
  column_to_rownames("sample")

datTraits <- read.csv2("DeSeq/data/sample_infos.csv") %>% 
  filter(species=="spimp") %>% 
  column_to_rownames("X") %>% 
  dplyr::select(infection, genotype, rep)

moduleColors <- read.delim("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/WGCNA/spimp/TOM_final_moduleColors_genids.txt", 
                           sep=" ", header=T)

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors$ModuleColor)$eigengenes
MEs = orderMEs(MEs0)
MEs <- MEs[, !colnames(MEs) %in% "MEgrey"]

datModules <- as.data.frame(MEs)
datModules$SampleID <- rownames(datModules)
longModules <- pivot_longer(
  datModules,
  cols = starts_with("ME"),
  names_to = "Module",
  values_to = "Eigengene"
)

mergedData <- merge(longModules, datTraits, by.x = "SampleID", by.y = "row.names")

# Standardize the Eigengene per module
mergedData <- mergedData %>%
  group_by(Module) %>%
  mutate(Eigengene_scaled = scale(Eigengene)) %>%
  ungroup()

# Ensure factors
mergedData$infection <- as.factor(mergedData$infection)
mergedData$genotype <- as.factor(mergedData$genotype)

# Initialize results
eta2_df <- data.frame(Module=character(),
                      Infection_eta2=double(),
                      Genotype_eta2=double(),
                      Interaction_eta2=double(),
                      stringsAsFactors=FALSE)

p_value_df <- data.frame()

# Loop through each module
for (module in unique(mergedData$Module)) {
  module_data <- subset(mergedData, Module == module)
  
  # Fit model with interaction
  model_full <- lm(Eigengene_scaled ~ infection * genotype, data = module_data)
  
  # Type III ANOVA for p-values
  aov_tab <- Anova(model_full, type = "III")
  
  # Extract p-values for infection, genotype, and interaction
  terms_of_interest <- c("infection", "genotype", "infection:genotype")
  aov_interest <- aov_tab[rownames(aov_tab) %in% terms_of_interest, , drop=FALSE]
  
  # Store p-values
  p_value_df <- rbind(
    p_value_df,
    data.frame(
      Module = module,
      Effect = rownames(aov_interest),
      p_value = aov_interest[,"Pr(>F)"],
      stringsAsFactors = FALSE
    )
  )
  
  # Compute non-partial eta²
  eta_res <- eta_squared(model_full, partial = FALSE)
  
  infection_eta     <- eta_res$Eta2[eta_res$Parameter == "infection"]
  genotype_eta      <- eta_res$Eta2[eta_res$Parameter == "genotype"]
  interaction_eta   <- eta_res$Eta2[eta_res$Parameter == "infection:genotype"]
  
  eta2_df <- rbind(eta2_df, data.frame(
    Module = module,
    Infection_eta2 = ifelse(length(infection_eta)>0, infection_eta, NA),
    Genotype_eta2 = ifelse(length(genotype_eta)>0, genotype_eta, NA),
    Interaction_eta2 = ifelse(length(interaction_eta)>0, interaction_eta, NA),
    stringsAsFactors = FALSE
  ))
}

# Correct p-values (Benjamini-Hochberg)
p_value_df <- p_value_df %>%
  group_by(Effect) %>%
  mutate(corrected_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Reshape eta2 for the heatmap
eta2_long <- eta2_df %>%
  pivot_longer(cols = c("Infection_eta2", "Genotype_eta2", "Interaction_eta2"), 
               names_to = "Effect", 
               values_to = "Eta2") %>%
  mutate(Effect = case_when(
    Effect == "Infection_eta2" ~ "infection",
    Effect == "Genotype_eta2" ~ "genotype",
    Effect == "Interaction_eta2" ~ "infection:genotype"
  ))

# Merge eta² with p-values
merged_annot <- p_value_df %>%
  left_join(eta2_long, by = c("Module", "Effect")) %>%
  mutate(annotation = paste0("η²=", round(Eta2,2), "\np=", signif(corrected_p_value, 2)))

# Create wide matrix for Eta²
eta2_wide <- eta2_long %>%
  pivot_wider(names_from = Effect, values_from = Eta2)

mat_eta2 <- as.matrix(eta2_wide[, -1])
rownames(mat_eta2) <- eta2_wide$Module

# Create wide annotation matrix
annot_wide <- merged_annot %>%
  select(Module, Effect, annotation) %>%
  pivot_wider(names_from = Effect, values_from = annotation)

display_mat <- as.matrix(annot_wide[, -1])
rownames(display_mat) <- annot_wide$Module

colnames(mat_eta2) <- c("Infection", "Genotype", "G:I")
colnames(display_mat) <- c("Infection", "Genotype", "G:I")

# Plot heatmap
png(paste0("WGCNA/documentation/pics/spimp/", Sys.Date(),
           "_heatmap_module_trait.png"),
    width = 18/2, height = 10, units = "cm", res = 900)
svg(paste0("WGCNA/documentation/pics/spimp/", Sys.Date(),
           "_heatmap_module_trait.svg"),
    width = 3.248031, height = 3.14961)
pheatmap(
  mat_eta2,
  cluster_rows = T,
  treeheight_row=15,
  legend=F,
  cluster_cols = F,
  display_numbers = display_mat,
  fontsize_number = 8,
  number_color = "black",
  angle_col = 0,
  fontsize = 11,
  border_color = "black",
  color = colorRampPalette(c("white", "#C22B26"))(100)
)
dev.off()

