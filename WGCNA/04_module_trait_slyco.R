rm(list=ls())
pacman::p_load(WGCNA, lme4, tidyverse, pheatmap, car, effectsize)
options(stringsAsFactors = FALSE)

setwd("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/")

# Read gene IDs
ids <- read.delim("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/PROTEOME/slyco_curated_proteome_OG_pannzer_dedub_ids.txt",
                  header=FALSE) %>%
  mutate(gene = gsub(">", "", V1),
         gene = gsub("GeneExt~|mRNA_", "", gene),
         gene = gsub("t\\.peak", "g.peak", gene),
         gene = gsub("t\\.minus", "g.minus", gene),
         gene = gsub("t\\.plus", "g.plus", gene),
         gene = gsub("\\.[1-9].*|\\.p[1-9].*", "", gene)) %>%
  select(gene)

# Read expression data
datExpr <- read.csv("DeSeq/data/norm_counts_all_rlog.csv") %>%
  filter(species == "S. lycopersicoides" & Geneid %in% ids$gene) %>%
  select(!species & !source) %>%
  pivot_wider(names_from = Geneid, values_from = normalized_count) %>%
  column_to_rownames("sample")

# Read sample traits
datTraits <- read.csv2("DeSeq/data/sample_infos.csv") %>%
  filter(species == "slyco") %>%
  column_to_rownames("X") %>%
  select(infection, genotype, rep)

# Read module colors
moduleColors <- read.table("WGCNA/data/colors/slyco_module_colors_TOM_filtered_genids.txt", header=TRUE) %>%
  select(ModuleColor)

# Calculate module eigengenes
MEs0 = moduleEigengenes(datExpr, moduleColors$ModuleColor)$eigengenes
MEs = orderMEs(MEs0)
MEs <- MEs[, !colnames(MEs) %in% "MEgrey"]

# Prepare data in long format
datModules <- as.data.frame(MEs) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(cols = starts_with("ME"), names_to = "Module", values_to = "Eigengene")

mergedData <- merge(datModules, datTraits, by.x = "SampleID", by.y = "row.names")

# Ensure factors
mergedData$infection <- as.factor(mergedData$infection)
mergedData$genotype <- as.factor(mergedData$genotype)

# Initialize result data frames
eta2_df <- data.frame(Module=character(),
                      Infection_eta2=double(),
                      Genotype_eta2=double(),
                      Interaction_eta2=double(),
                      stringsAsFactors=FALSE)

p_value_df <- data.frame()

# Loop through each module
for (module in unique(mergedData$Module)) {
  # Subset data for the current module
  module_data <- subset(mergedData, Module == module)
  
  # Fit model with interaction
  model_full <- lm(Eigengene ~ infection * genotype, data = module_data)
  
  # Type III ANOVA to get p-values
  aov_tab <- Anova(model_full, type = "III")
  
  # Terms of interest
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
  
  # Compute non-partial eta² using effectsize
  eta_res <- eta_squared(model_full, partial = FALSE)
  
  # Extract eta² for infection, genotype, and interaction
  infection_eta  <- eta_res$Eta2[eta_res$Parameter == "infection"]
  genotype_eta   <- eta_res$Eta2[eta_res$Parameter == "genotype"]
  interaction_eta <- eta_res$Eta2[eta_res$Parameter == "infection:genotype"]
  
  # Append results (use NA if the parameter doesn't exist)
  eta2_df <- rbind(eta2_df, data.frame(
    Module = module,
    Infection_eta2 = ifelse(length(infection_eta)>0, infection_eta, NA),
    Genotype_eta2 = ifelse(length(genotype_eta)>0, genotype_eta, NA),
    Interaction_eta2 = ifelse(length(interaction_eta)>0, interaction_eta, NA),
    stringsAsFactors = FALSE
  ))
}

# Correct p-values
p_value_df <- p_value_df %>%
  group_by(Effect) %>%
  mutate(corrected_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup()



# Reshape eta² data to long format
eta2_long <- eta2_df %>%
  pivot_longer(cols = c("Infection_eta2", "Genotype_eta2", "Interaction_eta2"),
               names_to = "Effect", 
               values_to = "Eta2") %>%
  mutate(Effect = case_when(
    Effect == "Infection_eta2" ~ "infection",
    Effect == "Genotype_eta2" ~ "genotype",
    Effect == "Interaction_eta2" ~ "infection:genotype"
  ))

# Combine p-values with eta²
merged_annot <- p_value_df %>%
  left_join(eta2_long, by = c("Module", "Effect")) %>%
  mutate(annotation = paste0("η²=", round(Eta2,2),
                             "\np=", signif(corrected_p_value, 2)))

# Create a wide matrix of Eta2 for pheatmap
eta2_wide <- eta2_long %>%
  pivot_wider(names_from = Effect, values_from = Eta2)

mat_eta2 <- as.matrix(eta2_wide[, -1]) # numeric matrix of Eta² values
rownames(mat_eta2) <- eta2_wide$Module

# Create a wide matrix of annotations
annot_wide <- merged_annot %>%
  select(Module, Effect, annotation) %>%
  pivot_wider(names_from = Effect, values_from = annotation)

display_mat <- as.matrix(annot_wide[, -1])
rownames(display_mat) <- annot_wide$Module

# Rename columns for neat heatmap labels
colnames(mat_eta2) <- c("Infection", "Genotype", "G:I")
colnames(display_mat) <- c("Infection", "Genotype", "G:I")

# Plot heatmap
png(paste0("WGCNA/documentation/pics/slyco/", Sys.Date(), 
           "_heatmap_module_trait_eta2_nonpartial.png"),
    width = 16.5/2, height = 8, units = "cm", res = 900)

svg(paste0("WGCNA/documentation/pics/slyco/", Sys.Date(),
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
