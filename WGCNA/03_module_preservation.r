# Test Module precervation between spimp and spen 

library(WGCNA)
library(tidyverse)

setwd("/gxfs_home/cau/suaph281/2024_solanum_ldt_rnaseq/")
options(stringsAsFactors = FALSE)
datExpr <- read.csv("DeSeq/data/norm_counts_all_rlog.csv")


datExpr_slyco <-  datExpr %>% 
  dplyr::filter(species=="S. lycopersicoides") %>% 
  select(!species & !source) %>% 
  pivot_wider(names_from = Geneid, values_from = normalized_count) %>%
  column_to_rownames("sample")

datExpr_spimp<-  datExpr %>% 
  dplyr::filter(species=="S. pimpinellifolium") %>% 
  select(!species & !source) %>% 
  pivot_wider(names_from = Geneid, values_from = normalized_count) %>%
  column_to_rownames("sample")


setLabels = c("S. pimpinellifolium", "S. lycopersicoides")

multiExpr=list(`S. pimpinellifolium`=list(data=datExpr_spimp),
           `S. lycopersicoides`=list(data=datExpr_slyco))

moduleColorsSpim=module_colors

multiColor=list(`S. pimpinellifolium`=moduleColorsSpim)
# Perservation analysis
nPermutations1=10
set.seed(1)
system.time({
    mp = modulePreservation(multiExpr, multiColor,
    referenceNetworks = 1, nPermutations = nPermutations1,
    networkType = "signed hybrid", corFnc="bicor", 
    randomSeed = 1234, quickCor = 0, verbose = 4, 
    checkData=F)
})
# Save the results of the module preservation analysis
save(mp, file = "modulePreservation.RData")
# If needed, reload the data:
load(file = "modulePreservation.RData")

# specify the reference and the test networks
ref=1; test = 2

Obs.PreservationStats= mp$preservation$observed[[ref]][[test]]
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
# Look at the observed preservation statistics
Obs.PreservationStats
# Z statistics from the permutation test analysis
Z.PreservationStats


# Viz: 

modColors = rownames(Obs.PreservationStats)
moduleSize = Obs.PreservationStats$moduleSize
# we will omit the grey module (background genes)
# and the gold module (random sample of genes)
selectModules = !(modColors %in% c("grey", "gold"))
# Text labels for points
point.label = modColors[selectModules]

#Composite preservation statistics
medianRank=Obs.PreservationStats$medianRank.pres
Zsummary=Z.PreservationStats$Zsummary.pres
# Open a PNG device for the first plot
png("WGCNA/documentation/pics/preservation/medianRank_Preservation_checkOff.png", width = 800, height = 600)
par(mfrow=c(1,1), mar = c(4.5, 4.5, 2.5, 1)) # Reset to single plot layout

# Plot medianRank versus module size
plot(moduleSize[selectModules], medianRank[selectModules], col = 1,
     bg = modColors[selectModules], pch = 21, main = "medianRank Preservation",
     cex = 2, ylab = "medianRank", xlab = "Module size", log = "x")
labelPoints(moduleSize[selectModules], medianRank[selectModules], point.label, cex = 1, offs = 0.03)

# Close the PNG device
dev.off()

# Open a PNG device for the second plot
png("WGCNA/documentation/pics/preservation/Zsummary_Preservation_checkOff.png", width = 800, height = 600)
par(mfrow=c(1,1), mar = c(4.5, 4.5, 2.5, 1)) # Reset to single plot layout

# Plot Zsummary versus module size
plot(moduleSize[selectModules], Zsummary[selectModules], col = 1,
     bg = modColors[selectModules], pch = 21, main = "Zsummary Preservation",
     cex = 2, ylab = "Zsummary", xlab = "Module size", log = "x")
labelPoints(moduleSize[selectModules], Zsummary[selectModules], point.label, cex = 1, offs = 0.03)

# Add threshold lines for Zsummary
abline(h = 0)
abline(h = 2, col = "blue", lty = 2)
abline(h = 10, col = "red", lty = 2)

# Close the PNG device
dev.off()

#Not much perservation between the two species??
# this journey end here! 
# because identical geneIDs are required for this approach to work, 
# all relevant publications are mapping all species on one reference genome. 
# Next: 12.6: consensus-modules 

# number of networks used in the consensus network analysis:
nSets = 2
# Vector with descriptive names of the two sets.
setLabels = c("S. pimpinellifolium", "S. lycopersicoides")
shortLabels = c("spimp", "slyco")
# Define a list whose components contain the data
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = datExpr_spimp)
names(multiExpr[[1]]$data) = names(datExpr_spimp)
rownames(multiExpr[[1]]$data) = dimnames(datExpr_spimp)[[1]]

multiExpr[[2]] = list(data = datExpr_slyco)
names(multiExpr[[2]]$data) = names(datExpr_slyco)
rownames(multiExpr[[2]]$data) = dimnames(datExpr_slyco)[[1]]
# Check that the data has the correct format:
exprSize = checkSets(multiExpr)

# The journey ends here? The tool expects same length of genes in the two data sets. 

# The variable exprSize contains useful information about the sizes of all
# of the data sets now we run automatic module detection procedure
netConsensus = blockwiseConsensusModules(multiExpr, maxBlockSize = 5000, power = 7, 
    minModuleSize = 30, deepSplit = 2, pamRespectsDendro = FALSE, mergeCutHeight = 0.25, 
    numericLabels = TRUE, minKMEtoStay = 0, saveTOMs = TRUE)