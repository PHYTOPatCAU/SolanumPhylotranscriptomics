# Module Trait correlations 
# using a linear model with scaled eigengenes
# FDR correction applied to p-values
# Severin Einspanier
# as of 2024_12_09
# Orthogroups

rm(list=ls()) 
pacman::p_load(tidyverse, WGCNA, reshape2, pheatmap)

options(stringsAsFactors = FALSE)
setwd("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/")

# --- Data Loading and Preparation ---
datExpr <- read.csv("OGs/data/orthogroup_expression_data_rlog_all_species.csv", header=TRUE) %>% 
  column_to_rownames("X") %>%
  t() %>% 
  as.data.frame()

ranked_genotypes <- read.csv2("DeSeq/data/sample_infos.csv") %>%
  group_by(accession) %>%
  summarise(mean_LDT = mean(lsmean.LDT.)) %>%
  arrange(mean_LDT) %>%
  mutate(rank = row_number())

datTraits <- read.csv2("DeSeq/data/sample_infos.csv") %>% 
  left_join(ranked_genotypes, by = "accession") %>%
  column_to_rownames("X") %>%
  select(infection, species, rank, rep)

datTraits$infection <- as.factor(datTraits$infection)
datTraits$species <- relevel(as.factor(datTraits$species), ref="slyco")

# Set rank as a factor with "1" as reference
datTraits$rank <- relevel(factor(datTraits$rank, levels = as.character(1:10)), ref="10")
datTraits$rep <- as.factor(datTraits$rep)
datTraits$infection <- relevel(datTraits$infection, ref="0")

moduleColors <- read.delim("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/WGCNA/OGs/TOM_sft_9_ds0_mch35_moduleColors_genid.txt", 
                           header=TRUE, sep=" ")

# Calculate module eigengenes
MEs0 <- moduleEigengenes(datExpr, moduleColors$ModuleColor, impute=TRUE, excludeGrey=TRUE)$eigengenes
MEs <- orderMEs(MEs0)

datModules <- as.data.frame(MEs)
datModules$SampleID <- rownames(datModules)

longModules <- pivot_longer(
  datModules,
  cols = starts_with("ME"),
  names_to = "Module",
  values_to = "Eigengene"
)

mergedData <- merge(longModules, datTraits, by.x = "SampleID", by.y = "row.names")

# Scale the response (Eigengene)
mergedData <- mergedData %>%
  group_by(Module) %>%
  mutate(Eigengene_scaled = scale(Eigengene)) %>%
  ungroup()

# Fit simple linear models with scaled Eigengene
model_results <- list()

# Data frame to store partial R² and full R²
r2_df <- data.frame(Module = character(),
                    Full_R2 = double(),
                    Infection_partial_R2 = double(),
                    Rank_partial_R2 = double(),
                    stringsAsFactors = FALSE)

for (module in unique(mergedData$Module)) {
  module_data <- subset(mergedData, Module == module)
  
  # Full model
  model_full <- lm(Eigengene_scaled ~ infection + rank, data = module_data)
  R2_full <- summary(model_full)$r.squared
  
  # No infection
  model_no_infection <- lm(Eigengene_scaled ~ rank, data = module_data)
  R2_no_infection <- summary(model_no_infection)$r.squared
  
  # No rank
  model_no_rank <- lm(Eigengene_scaled ~ infection, data = module_data)
  R2_no_rank <- summary(model_no_rank)$r.squared
  
  infection_partial_R2 <- R2_full - R2_no_infection
  rank_partial_R2 <- R2_full - R2_no_rank
  
  model_results[[module]] <- model_full
  
  r2_df <- rbind(r2_df, data.frame(
    Module = module,
    Full_R2 = R2_full,
    Infection_partial_R2 = infection_partial_R2,
    Rank_partial_R2 = rank_partial_R2
  ))
}

# Extract fixed effects for all modules
fixed_effects_df <- data.frame()

for (module in names(model_results)) {
  model <- model_results[[module]]
  fixed <- as.data.frame(coef(summary(model)))
  fixed$Effect <- rownames(fixed)
  fixed$Module <- module
  fixed_effects_df <- rbind(fixed_effects_df, fixed)
}

colnames(fixed_effects_df) <- c("Estimate","StdError","tValue","Pr(>|t|)","Effect","Module")

df <- fixed_effects_df %>% 
  filter(Effect != "(Intercept)")

# Apply FDR correction
df$Padj <- p.adjust(df$`Pr(>|t|)`, method = "BH")

# Create text annotation with Estimates and adjusted p-values
#df$Text <- paste0(round(df$Estimate, 2), "\n(p=", signif(df$Padj, 2), ")")
df$Text <- ifelse(df$Padj < 0.01, paste0(round(df$Estimate, 2), "\n(p < 0.01)"), paste0(round(df$Estimate, 2), "\n(p=", signif(df$Padj, 2), ")"))
# Reshape for heatmap
heatmap_data <- dcast(df, Module ~ Effect, value.var = "Estimate") 
text_data <- dcast(df, Module ~ Effect, value.var = "Text")

rownames(heatmap_data) <- heatmap_data$Module
heatmap_data$Module <- NULL
rownames(text_data) <- text_data$Module
text_data$Module <- NULL

# Order columns for heatmap
desired_order <- c("infection1", paste0("rank", 1:9))  # Specify desired column order
current_columns <- colnames(heatmap_data)

# Reorder the columns if they exist in the current dataset
ordered_columns <- intersect(desired_order, current_columns)  # Keep only existing columns
heatmap_data <- heatmap_data[, ordered_columns, drop = FALSE]
text_data <- text_data[, ordered_columns, drop = FALSE]

colnames(heatmap_data)=c("Infected", "Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Rank8", "Rank9")
colnames(text_data)=c("Infected", "Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Rank8", "Rank9")

# Generate heatmap with Padj and Estimate
png(paste0("WGCNA/documentation/pics/OG/", Sys.Date(), "_heatmap_module_trait_correlation_scaled_padj_estimate.png"),     
    width = 15, height = 8, units = "cm", res = 900)
#svg(paste0("C://Users/suaph281/Nextcloud/ResiDEvo/sequencing/figures/fig_5/", Sys.Date(), "_heatmap_module_trait_correlation_scaled_padj_estimate.svg"),
 #   width = 15, height = 8)
pheatmap(
  as.matrix(heatmap_data),
  legend = F,
  treeheight_row=15,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = as.matrix(text_data),
  #main = "OG-Network Module-Trait Relationships \n(Scaled Estimates and FDR-corrected p-values)",
  fontsize_number = 8,
  number_color="black",
  angle_col=0,
  fontsize = 11,
  border_color = "black"
)
dev.off()

# Save R² table
r2_summary_file <- paste0("WGCNA/documentation/pics/OG/", Sys.Date(), "_R2_summary_table.csv")
write.csv(r2_df, file = r2_summary_file, row.names = FALSE)

# Preview saved R² table
head(r2_df)
