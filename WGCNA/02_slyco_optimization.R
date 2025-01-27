# Find networksettings for S. lycopersicoides
# module load R/4.3.1 gcc/12.3.0
# This is the script to define networks for spimp on the nesh. 

rm(list=ls())
library(tidyverse)
library(WGCNA)
enableWGCNAThreads(n=30)

setwd("/gxfs_home/cau/suaph281/2024_solanum_ldt_rnaseq/")
options(stringsAsFactors = FALSE)

ids <- read.delim("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/PROTEOME/slyco_curated_proteome_OG_pannzer_dedub_ids.txt",
    header=F) %>%
    mutate(gene=gsub(">", "", V1)) %>%
    mutate(gene=gsub("GeneExt~", "", gene))%>% 
    mutate(gene=gsub("mRNA_", "", gene))%>%  
    mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
    mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
    mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
    mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))%>%
    select(gene)


# filter the gene expression:

datExpr <- read.csv("DeSeq/data/norm_counts_all_rlog.csv") %>%
  dplyr::filter(species=="S. lycopersicoides" & Geneid %in% ids$gene) %>% 
  select(!species & !source) %>% 
  pivot_wider(names_from = Geneid, values_from = normalized_count) %>%
  column_to_rownames("sample")

# Get Outliers

# sampleTree = hclust(dist(datExpr), method = "average");
# # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# # The user should change the dimensions if the window is too large or too small.
# sizeGrWindow(12,9)
# png(file = "WGCNA/documentation/pics/slyco/filtered_wgcna/sampleClustering.png", width = 1200, height = 900);
# par(cex = 0.6);
# par(mar = c(5, 4, 2, 1))  # Adjusted margins
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
#     cex.axis = 1.5, cex.main = 2)
# dev.off()

# # # Compute SFT 
# powers = c(c(1:10), seq(from = 12, to=30, by=1))
# # Call the network topology analysis function
# sft = pickSoftThreshold(datExpr, removeFirst=F,
#                         networkType = "signed hybrid",
#                         blockSize = ncol(datExpr), 
#                         powerVector = powers, verbose = 5)
# # Plot the results:
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# png(file = "WGCNA/documentation/pics/slyco/filtered_wgcna/sft_combined.png", width = 1800, height = 900);
# par(mfrow = c(1, 2))  # Set up the plotting area to have 1 row and 2 columns

# # First plot: Scale independence
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n",
#     main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#     labels=powers, cex=cex1, col="red");
# abline(h=0.90, col="red");

# # Second plot: Mean connectivity
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
#     main = paste("Mean connectivity"));
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red");

# dev.off()

# # get network 


## 1. DeepSplit

# I might want to consider using multithreading for this.
# enableWGCNAThreads(n=30)

 findDeepSplit <- function(datExpr, deepSplit_values) {
   results <- list()
   colors_list <- list()
 
   for (deep_split in deepSplit_values) {
     result <- blockwiseModules(datExpr, checkMissingData = T, 
                                 maxBlockSize = ncol(datExpr), 
                                networkType = "signed hybrid",
                                power =9, 
                                TOMType = "signed", 
                                minModuleSize = 30,
                                corType = "bicor", 
                                reassignThreshold = 0.01, 
                                mergeCutHeight = 0.2, 
                                deepSplit = deep_split,
                                detectCutHeight = 0.99,
                                numericLabels = TRUE,
                                saveTOMs = FALSE, 
                                verbose = 3)
     key <- paste0("DeepSplit_", deep_split)
     results[[key]] <- result
     colors_list[[key]] <- labels2colors(result$colors)
   }
  
   return(list(results = results, colors_list = colors_list))
 }

# ## Example usage
 deepSplit_values <- c(0, 1, 2, 3, 4)
 results <- findDeepSplit(datExpr, deepSplit_values)

results_col <- results$colors_list
results_dend <- results$results

# Plot the dendrogram and module colors for the first combination
first_key <- names(results_dend)[1]

png(file = "WGCNA/documentation/pics/slyco/filtered_wgcna/DeepSplit_sft_9.png", width = 10000, height = 8000, 
     res=600);
 plotDendroAndColors(results_dend[[first_key]]$dendrograms[[1]], 
                     do.call(cbind, results_col), 
                     c(names(results_col)), 
                     dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05, 
                     autoColorHeight = FALSE,
                     colorHeight = 0.65)
 dev.off()


#stop("End of DeepSplit-Determination. Evaluate output and preceed with mergeCutheight")

# Take DS=1

## 2. Find mergeCutHeight

find_mergeCutHeight <- function(datExpr, variable_values) {
  results <- list()
  colors_list <- list()
  
  for (merge_cut_height in variable_values) {
    result <- blockwiseModules(datExpr, checkMissingData = FALSE, replaceMissingAdjacencies = TRUE, 
                               maxBlockSize = ncol(datExpr), 
                               networkType = "signed hybrid",
                               power = 9, 
                               TOMType = "signed", 
                               minModuleSize = 30,
                               corType = "bicor", 
                               reassignThreshold = 0.01, 
                               mergeCutHeight = merge_cut_height, 
                               deepSplit = 1,
                               detectCutHeight = 0.99,
                               numericLabels = TRUE,
                               saveTOMs = FALSE, 
                               #saveTOMFileBase = "C:/Users/suaph281/Desktop/nesh_local/LDT_RNAseq/WGCNA/TOM_SE_OG_sft_4",
                               verbose = 3)
    results[[paste0("Variable_set_", merge_cut_height)]] <- result
    colors_list[[paste0("Variable_set_", merge_cut_height)]] <- labels2colors(result$colors)
  }

  
  return(list(results = results, colors_list = colors_list))
}

# # Example usage
variable_values <- c(0.1, 0.15, 0.2 , 0.25, 0.3)
# #datExpr <- read.csv("path_to_your_data.csv")  # Replace with your actual data
results <- find_mergeCutHeight(datExpr, variable_values)

results_col<- results$colors_list
results_dend <- results$results

png(file = "WGCNA/documentation/pics/slyco/filtered_wgcna/MergeCut_sft_9_ds1.png", width = 10000, height = 8000, 
     res=600);

#   # Plot the dendrogram and module colors
 plotDendroAndColors(results_dend[[1]]$dendrograms[[1]], 
                       do.call(cbind, results_col), 
                       c(paste0("Variable_set_", variable_values)), 
                       dendroLabels = FALSE, hang = 0.03,
                       addGuide = TRUE, guideHang = 0.05, 
                     autoColorHeight = F,
                     colorHeight = .65)

 dev.off()
# stop("eval output and make final dendrogram")

# # Extract module labels (numeric)
#module_labels <- result$colors
# head(module_labels)  # View first few labels

# # Convert numeric labels to colors
# module_colors <- labels2colors(module_labels)
# head(module_colors)  # View first few colors

# # Create a data frame linking gene IDs to module colors
# gene_module_df <- data.frame(
#   GeneID = colnames(datExpr),  # Replace with your gene IDs
#   ModuleColor = module_colors
# )
# head(gene_module_df)
# write.table(gene_module_df,"/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/slyco/TOM/module_colors_TOM_filtered_genids.txt")


# Final PLots:

result <- blockwiseModules(datExpr, checkMissingData = T, 
  replaceMissingAdjacencies = F, 
                               maxBlockSize = ncol(datExpr), 
                               networkType = "signed hybrid",
                               power =9, 
                               TOMType = "signed", 
                               minModuleSize = 30,
                               corType = "bicor", 
                               reassignThreshold = 0.01, 
                               mergeCutHeight = 0.2, 
                               deepSplit = 1,
                               detectCutHeight = 0.95,
                               numericLabels = TRUE,
                               saveTOMs = TRUE, 
                               saveTOMFileBase = "/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/slyco/TOM/TOM_filtered_sft9",
                               verbose = 3)
  
colors=labels2colors(result$colors)


# Plot final figure:

png(file = "WGCNA/documentation/pub/slyco_network_filtered_sft9.png", width = 2000, height = 2000, 
    res=600);

svg(file = "WGCNA/documentation/pub/slyco_network_filtered_sft9.svg", width = 7, height = 5)

# Plot the dendrogram and module colors
plotDendroAndColors(result$dendrograms[[1]], 
                      colors[result$blockGenes[[1]]], 
                      c(""), 
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = F, guideHang = 0.05, 
                    autoColorHeight = T,
                    colorHeight = .2)

dev.off()


moduleColors=labels2colors(result$colors)

# Extract module labels (numeric)
module_labels <- result$colors
head(module_labels)  # View first few labels

# Convert numeric labels to colors
module_colors <- labels2colors(module_labels)
head(module_colors)  # View first few colors

# Create a data frame linking gene IDs to module colors
gene_module_df <- data.frame(
  GeneID = colnames(datExpr),  # Replace with your gene IDs
  ModuleColor = module_colors
)
head(gene_module_df)
write.table(gene_module_df,"/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/slyco/TOM/module_colors_sft9_TOM_filtered_genids.txt")
