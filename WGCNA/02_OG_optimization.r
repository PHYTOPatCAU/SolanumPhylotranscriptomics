# OG Optimization 

# This is the script to define networks for OGs on the nesh. 

rm(list=ls())
library(tidyverse)
library(WGCNA)

setwd("/gxfs_home/cau/suaph281/2024_solanum_ldt_rnaseq/")
options(stringsAsFactors = FALSE)
datExpr <- read.csv("OGs/data/orthogroup_expression_data_rlog_all_species.csv")

datExpr <-  datExpr %>% 
  column_to_rownames("X")

datExpr <- as.data.frame(t(datExpr))

datTraits <- read.csv2("DeSeq/data/sample_infos.csv")

datTraits_cor <- datTraits %>% 
  column_to_rownames("X") %>% 
  select(genotype, infection, species_short, lsmean.LDT.)

# Get Outliers
sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
png(file = "WGCNA/documentation/pics/OG/sampleClustering_newProt.png", width = 1200, height = 900);
par(cex = 0.6);
par(mar = c(5, 4, 2, 1))  # Adjusted margins
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
    cex.axis = 1.5, cex.main = 2)
dev.off()


# # Compute SFT 



powers = c(c(1:10), seq(from = 12, to=30, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, removeFirst=F,
                        networkType = "signed hybrid",
                        blockSize = ncol(datExpr), 
                        powerVector = powers, verbose = 5)
# Plot the results:
#sizeGrWindow(9, 5)
#par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
png(file = "WGCNA/documentation/pics/OG/sft_combined_newProt.png", width = 1800, height = 900);
par(mfrow = c(1, 2))  # Set up the plotting area to have 1 row and 2 columns

# First plot: Scale independence
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n",
    main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers, cex=cex1, col="red");
abline(h=0.90, col="red");

# Second plot: Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"));
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red");

dev.off()

# get network 

## 1. DeepSplit

# I might want to consider using multithreading for this.
enableWGCNAThreads(n=30)

findDeepSplit <- function(datExpr, deepSplit_values) {
  results <- list()
  colors_list <- list()
  
  for (deep_split in deepSplit_values) {
    result <- blockwiseModules(datExpr, checkMissingData = FALSE, replaceMissingAdjacencies = TRUE, 
                               maxBlockSize = ncol(datExpr), 
                               networkType = "signed hybrid",
                               power =9, 
                               TOMType = "signed", 
                               minModuleSize = 30,
                               corType = "bicor", 
                               reassignThreshold = 0.01, 
                               mergeCutHeight = 0.25, 
                               deepSplit = deep_split,
                               detectCutHeight = 0.99,
                               numericLabels = TRUE,
                               saveTOMs = FALSE, 
                               verbose = 1)
    key <- paste0("DeepSplit_", deep_split)
    results[[key]] <- result
    colors_list[[key]] <- labels2colors(result$colors)
  }
  
  return(list(results = results, colors_list = colors_list))
}

# Example usage
deepSplit_values <- c(0, 1, 2, 3, 4)
results <- findDeepSplit(datExpr, deepSplit_values)

results_col <- results$colors_list
results_dend <- results$results

# Plot the dendrogram and module colors for the first combination
first_key <- names(results_dend)[1]

png(file = "WGCNA/documentation/pics/OG/DeepSplit_sft_9_new.png", width = 10000, height = 8000, 
    res=600);
plotDendroAndColors(results_dend[[first_key]]$dendrograms[[1]], 
                    do.call(cbind, results_col), 
                    c(names(results_col)), 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    autoColorHeight = FALSE,
                    colorHeight = 0.65)

dev.off()


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
                               deepSplit = 0,
                               detectCutHeight = 0.99,
                               numericLabels = TRUE,
                               saveTOMs = FALSE, 
                               #saveTOMFileBase = "C:/Users/suaph281/Desktop/nesh_local/LDT_RNAseq/WGCNA/TOM_SE_OG_sft_4",
                               verbose = 1)
    results[[paste0("Variable_set_", merge_cut_height)]] <- result
    colors_list[[paste0("Variable_set_", merge_cut_height)]] <- labels2colors(result$colors)
  }

  
  return(list(results = results, colors_list = colors_list))
}

# Example usage
variable_values <- c(0.1, 0.15, 0.2 , 0.25, 0.3, 0.35, 0.4, 0.5)
#datExpr <- read.csv("path_to_your_data.csv")  # Replace with your actual data
results <- find_mergeCutHeight(datExpr, variable_values)

results_col<- results$colors_list
results_dend <- results$results

png(file = "WGCNA/documentation/pics/OG/MergeCut_sft_9_ds0_new.png", width = 10000, height = 8000, 
    res=600);

  # Plot the dendrogram and module colors
plotDendroAndColors(results_dend[[1]]$dendrograms[[1]], 
                      do.call(cbind, results_col), 
                      c(paste0("Variable_set_", variable_values)), 
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, 
                    autoColorHeight = F,
                    colorHeight = .65)

dev.off()

# 3. Get final network:


result <- blockwiseModules(datExpr, checkMissingData = FALSE, replaceMissingAdjacencies = TRUE, 
                           maxBlockSize = ncol(datExpr), 
                           networkType = "signed hybrid",
                           power = 9, 
                           TOMType = "signed", 
                           minModuleSize = 30,
                           corType = "bicor", 
                           reassignThreshold = 0.01, 
                           mergeCutHeight = 0.35, 
                           deepSplit = 0,
                           detectCutHeight = 0.95,
                           numericLabels = TRUE,
                           saveTOMs = F, 
                           saveTOMFileBase = "/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/OG/TOM_sft_9_ds0_cm3",
                           verbose = 2)

# Get the module colors
module_colors <- labels2colors(result$colors)

write.table(module_colors, "/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/OG/TOM_sft_9_ds0_mch35_moduleColors.txt", 
col.names = FALSE, row.names = TRUE, quote = FALSE, sep = "\t")


png(file = "WGCNA/documentation/pics/OG/network_sft9_ds0_ch35.png", width = 10000, height = 6000, 
    res=600); 

  # Plot the dendrogram and module colors
plotDendroAndColors(result$dendrograms[[1]], 
                      module_colors, 
                      c("SFT_6_DS1_CH2"), 
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, 
                    autoColorHeight = F,
                    colorHeight = .15)

dev.off()

# Correlating modules with traits


