# QDR RNAseq Solanum species
# Severin Einspanier
# Perform basic DeSeq Analysis on all samples
# 2024_09

rm(list=ls())
library(DESeq2)
library(tidyverse)


### FUNCTIONS

perform_deseq_analysis <- function(counts_path, coldata_path, colnames_vec_path, contrast_name) {
  # Check if all required variables are provided
  if (missing(counts_path) || missing(coldata_path) || missing(contrast_name)) {
    stop("All arguments (counts_path, coldata_path, contrast_name) must be provided.")
  }
  # Load the counts table and sample information
  cts <- as.matrix(read.delim(counts_path, sep = "\t", row.names = "Geneid"))
  coldata <- read.csv2(coldata_path, row.names = 1)
  
  # Set the column names of the counts table
  colnames(cts) <- gsub(".*\\.([A-H][0-9]{2})_.*", "\\1", colnames(cts))
  
  
  # Ensure that the row names of the sample information match the column names of the counts table
  if (!all(rownames(coldata) == colnames(cts))) {
    stop("Row names of the sample information do not match the column names of the counts table.")
  }
  
  # Convert genotype and treatment to factors
  coldata$genotype <- as.factor(coldata$genotype)
  coldata$treatment <- as.factor(coldata$treatment)
  
  # Create DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~genotype * treatment)
  # Set the reference level for treatment
  dds$treatment <- factor(dds$treatment, levels = c("mock", "sclero"))
  
  dds$contrasts <- factor(paste0(dds$genotype, dds$treatment))
  design(dds) <- ~0+contrasts
  dds <- DESeq2::DESeq(dds)
  # Prefilter the dataset
  smallestGroupSize <- 4
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize 
  dds <- dds[keep, ]
  DESeq2::resultsNames(dds)
  
  # ORDER: contrasts1mock contrasts1sclero contrasts2mock contrasts2sclero
  # 2 meaning resistant genotype!
  
  TreatmentMatrix <- rbind(
    "infected_comp_genotypes" = c(0, -1, 0, 1),
    "mock_comp_genotypes" = c(-1, 0, 1, 0),
    "both_comp_infection" = c(-0.5, 0.5, -0.5, 0.5),
    "Resistant_comp_infection" = c(0,0,-1,1),
    "Susceptible_comp_infection" = c(-1,1,0,0) 
    
  )
  # Check if the contrast_name exists in the TreatmentMatrix
  if (!contrast_name %in% rownames(TreatmentMatrix)) {
    stop("Invalid contrast name provided.")
  }
  
  # Return the selected contrast
  selected_contrast <- TreatmentMatrix[contrast_name, ]
  
  
  # Treatments
  return_contrast <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type="ashr")
  
  return(return_contrast)
}


# Get actual DEGs information 

get_DEGs <- function(df, species){
  key <- intersect(rownames(df)[which(abs(df$log2FoldChange) > 1)] , rownames(df)[which(df$padj<=0.05)])
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  results$species <- species
  return(results)
}

setwd("C:/Users/suaph281/Desktop/")
setwd("C:/Users/Sever//Desktop/")

# Define the paths for each species
species_paths <- data.frame(
  species = c("S. pennellii", "S. chilense", "S. pimpinellifolium", "S. habrochaites", "S. lycopersicoides"),
  counts_table_path = c(
    "counts_table_matrix_spen",
    "counts_table_matrix_schil",
    "counts_table_matrix_spim",
    "counts_table_matrix_shabr",
    "Gcounts_table_matrix_slyco"
  ),
  sample_info_path = c(
    "sample_infos_spen.csv",
    "sample_infos_schil.csv",
    "sample_infos_spimp.csv",
    "sample_infos_shabro.csv",
    "sample_infos_slyco.csv"
  ),
  stringsAsFactors = FALSE
)

# Initialize lists to store results
results_DEGs_list <- list()
DEGs_all_list <- list()
results_list <- list()

# Iterate over each species and perform DESeq analysis
for (i in 1:nrow(species_paths)) {
  species <- species_paths$species[i]
  counts_table_path <- species_paths$counts_table_path[i]
  sample_info_path <- species_paths$sample_info_path[i]
  colnames_vec_path <- species_paths$colnames_vec_path[i]
  
  #infected_comp_genotypes, mock_comp_genotypes, both_comp_infection, Resistant_comp_infection, Susceptible_comp_infection  
  
  result <- perform_deseq_analysis(counts_table_path, 
                                   sample_info_path, 
                                   colnames_vec_path, 
                                   "Susceptible_comp_infection")

  DEGs_all_list[[species]] <- get_DEGs(result, species)
  result_out <- as.data.frame(result)
  result_out$species = species
  results_list[[species]] <- result_out
}
#combined_results_DEGs <- do.call(rbind, DEGs_all_list)

combined_results_DEGs <- do.call(rbind, lapply(DEGs_all_list, function(df) {
  df$GeneID <- rownames(df)
  rownames(df) <- NULL
  return(df)
}))

combined_results_all_genes <- do.call(rbind, lapply(results_list, function(df) {
  df$GeneID <- rownames(df)
  rownames(df) <- NULL
  return(df)
}))


# Access results for each species
results_DEGs_pen <- results_DEGs_list[["S. pennellii"]]
DEG_names_spen_up <- DEG_names_up_list[["S. pennellii"]]
DEG_names_spen_down <- DEG_names_down_list[["S. pennellii"]]
spen_DEGs_all <- DEGs_all_list[["S. pennellii"]]

results_DEGs_schil <- results_DEGs_list[["S. chilense"]]
DEG_names_schil_up <- DEG_names_up_list[["S. chilense"]]
DEG_names_schil_down <- DEG_names_down_list[["S. chilense"]]
schil_DEGs_all <- DEGs_all_list[["S. chilense"]]

results_DEGs_spimp <- results_DEGs_list[["S. pimpinellifolium"]]
DEG_names_spimp_up <- DEG_names_up_list[["S. pimpinellifolium"]]
DEG_names_spimp_down <- DEG_names_down_list[["S. pimpinellifolium"]]
spim_DEGs_all <- DEGs_all_list[["S. pimpinellifolium"]]

results_DEGs_shabro <- results_DEGs_list[["S. habrochaites"]]
DEG_names_shabro_up <- DEG_names_up_list[["S. habrochaites"]]
DEG_names_shabro_down <- DEG_names_down_list[["S. habrochaites"]]
shabr_DEGs_all <- DEGs_all_list[["S. habrochaites"]]

results_DEGs_slyco <- results_DEGs_list[["S. lycopersicoides"]]
DEG_names_slyco_up <- DEG_names_up_list[["S. lycopersicoides"]]
DEG_names_slyco_down <- DEG_names_down_list[["S. lycopersicoides"]]
slyco_DEGs_all <- DEGs_all_list[["S. lycopersicoides"]]


## write RESULTS

write.csv(combined_results_all_genes, "DeSeq_OUT/combined_SUS_inf_mock.csv")
write.csv(combined_results_DEGs, "DeSeq_OUT/combined_SUS_inf_mock_DEGs.csv")
