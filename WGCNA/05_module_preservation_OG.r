# Preservation testing of modules to identify modules that are of high condifence
# Orthogroup network
# Severin Einspanier
# 2024_12_01

rm(list=ls())
library(tidyverse)
library(WGCNA)

setwd("")
options(stringsAsFactors = FALSE)
datExpr <- read.csv("OGs/data/orthogroup_expression_data_rlog_all_species.csv")

datExpr <-  datExpr %>% 
  column_to_rownames("X")

datExpr <- as.data.frame(t(datExpr))

datTraits <- read.csv2("DeSeq/data/sample_infos.csv")

datTraits_cor <- datTraits %>% 
  column_to_rownames("X")

datTraits_cor

enableWGCNAThreads(n=30)

# Sample subsets for genotype and treatment combinations
spen_mock <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$species_short == 4 & datTraits$treatment == "mock"],]%>% as.matrix() 
spen_sclero <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$species_short == 4 & datTraits$treatment == "sclero"],]%>% as.matrix() 
schil_mock <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$species_short == 1 & datTraits$treatment == "mock"],]%>% as.matrix() 
schil_sclero <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$species_short == 1 & datTraits$treatment == "sclero"],]%>% as.matrix() 
slyco_mock <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$species_short == 3 & datTraits$treatment == "mock"],]%>% as.matrix() 
slyco_sclero <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$species_short == 3 & datTraits$treatment == "sclero"],]%>% as.matrix() 
shabro_mock <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$species_short == 2 & datTraits$treatment == "mock"],]%>% as.matrix() 
shabro_sclero <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$species_short == 2 & datTraits$treatment == "sclero"],]%>% as.matrix() 
spimp_mock <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$species_short == 5 & datTraits$treatment == "mock"],]%>% as.matrix() 
spimp_sclero <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$species_short == 5 & datTraits$treatment == "sclero"],]%>% as.matrix() 


bad_genes=""

# Filter the input matrix and subsets
input_in <- datExpr %>% as.matrix()


modules <- read.table("data/WGCNA/OG/TOM/TOM_sft_9_ds0_mch35_moduleColors_genid.txt", 
    sep=" ", row.names=NULL) %>%
    column_to_rownames("GeneID")

clrs_in <- modules[!(rownames(modules) %in% bad_genes), 2]
names(clrs_in) <- rownames(modules[!(rownames(modules) %in% bad_genes), ])

mp <- modulePreservation(
  multiData = list(
    all = list(data = input_in),
    spen_Mock = list(data = spen_mock),
    spen_Sclero = list(data = spen_sclero),
    schil_Mock = list(data = schil_mock),
    schil_Sclero = list(data = schil_sclero),
    slyco_Mock = list(data = slyco_mock),
    slyco_Sclero = list(data = slyco_sclero),
    shabro_Mock = list(data = shabro_mock),
    shabro_Sclero = list(data = shabro_sclero),
    spimp_Mock = list(data = spimp_mock),
    spimp_Sclero = list(data = spimp_sclero)
  ),
  multiColor = list(all = clrs_in),   # Adjust module colors as needed
  dataIsExpr = TRUE,
  referenceNetworks = 1,  # Use all as reference
  nPermutations = 10,
  maxModuleSize = 5000, 
  randomSeed = 1,
  quickCor = 0,
  parallelCalculation = TRUE,
  verbose = 3,
  networkType = "signed hybrid",
  corFnc = "bicor", 
  checkData=F
)

# Calculate Z statistics for all subsets
zstat <- bind_rows(
  spen_Mock = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.spen_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "spen_Mock", 
           species = "Solanum pennellii"),
  spen_Sclero = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.spen_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "spen_Sclero",
            species = "Solanum pennellii"),
  schil_Mock = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.schil_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "schil_Mock", 
    species="Solanum chilense"),
  schil_Sclero = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.schil_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "schil_Sclero", 
    species="Solanum chilense"),
  slyco_Mock = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.slyco_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "slyco_Mock", 
    species="Solanum lycopersicoides"),
  slyco_Sclero = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.slyco_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "slyco_Sclero", 
    species="Solanum lycopersicoides"),
  shabro_Mock = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.shabro_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "shabro_Mock", 
    species="Solanum habrochaites"),
  shabro_Sclero = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.shabro_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "shabro_Sclero", 
    species="Solanum habrochaites"),
  spimp_Mock = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.spimp_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "spimp_Mock",
    species="Solanum pimpinellifolium"),
  spimp_Sclero = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.spimp_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "spimp_Sclero",
species = "Solanum pimpinellifolium")
) 

write.table(zstat, "/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/OG/zstat/2024_12_01_zstat.tsv", sep="\t", row.names = FALSE)

png("WGCNA/documentation/pics/OG/presevation_zstats.png", width = 800, height = 400)

# Plot Z summary preservation statistics
ggplot(zstat, aes(x = moduleSize, y = Zsummary.pres,
                  shape = species, fill = moduleColor)) +
  geom_point(size = 3) +
  ylab("Preservation (Zsummary)") +
  scale_fill_manual(values = as.character(zstat$moduleColor),
                    guide = "none") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)) +
  geom_hline(yintercept = 10, col = "red") +
  geom_hline(yintercept = 2, col = "darkgreen") +
  theme_bw() +
  labs(title = "Module Preservation (Zsummary)", x = "Module Size")
dev.off()

# Calculate median rank statistics for all subsets
med_rank_stat <- bind_rows(
  spen_Mock = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.spen_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "spen_Mock", 
    species = "Solanum pennellii"),
  spen_Sclero = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.spen_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "spen_Sclero", 
    species = "Solanum pennellii"),
  schil_Mock = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.schil_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "schil_Mock", 
    species = "Solanum chilense"),
  schil_Sclero = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.schil_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "schil_Sclero", 
    species = "Solanum chilense"),
  slyco_Mock = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.slyco_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "slyco_Mock", 
    species = "Solanum lycopersicoides"),
  slyco_Sclero = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.slyco_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "slyco_Sclero", 
    species = "Solanum lycopersicoides"),
  shabro_Mock = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.shabro_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "shabro_Mock", 
    species = "Solanum habrochaites"),
  shabro_Sclero = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.shabro_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "shabro_Sclero", 
    species = "Solanum habrochaites"),
  spimp_Mock = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.spimp_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "spimp_Mock", 
    species = "Solanum pimpinellifolium"),
  spimp_Sclero = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.spimp_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "spimp_Sclero", 
    species = "Solanum pimpinellifolium")
) 

# Plot median rank preservation statistics
png("WGCNA/documentation/pics/OG/presevation_median_rank.png", width = 800, height = 400)

ggplot(med_rank_stat, aes(x = moduleSize, y = medianRank.pres,
                          shape = species, fill = moduleColor)) +
  geom_point(size = 3) +
  scale_y_reverse() +
  ylab("Preservation (Median rank)") +
  scale_fill_manual(values = as.character(med_rank_stat$moduleColor),
                    guide = "none") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)) +
  theme_bw() +
  labs(title = "Module Preservation (Median Rank)", x = "Module Size")

dev.off()
