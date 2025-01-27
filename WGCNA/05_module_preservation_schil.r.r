# Preservation testing of modules to identify modules that are of high condifence
# Solanum chilense
# Severin Einspanier
# 2024_12_01

rm(list=ls())
library(tidyverse)
library(WGCNA)

setwd("/gxfs_home/cau/suaph281/2024_solanum_ldt_rnaseq/")
options(stringsAsFactors = FALSE)
datExpr <- read.csv("DeSeq/data/norm_counts_all_rlog.csv") %>% 
  dplyr::filter(species=="S. chilense") %>% 
  select(!species & !source) %>% 
  pivot_wider(names_from = Geneid, values_from = normalized_count) %>%
  column_to_rownames("sample")

datTraits <- read.csv2("DeSeq/data/sample_infos.csv") %>% 
  filter(species=="schil") 

enableWGCNAThreads(n=30)

# Sample subsets for genotype and treatment combinations
genotype1_mock <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$genotype == 1 & datTraits$treatment == "mock"],]
genotype1_sclero <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$genotype == 1 & datTraits$treatment == "sclero"],]
genotype2_mock <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$genotype == 2 & datTraits$treatment == "mock"],]
genotype2_sclero <- datExpr[rownames(datExpr) %in% datTraits$X[datTraits$genotype == 2 & datTraits$treatment == "sclero"],]


bad_genes <- unique(c(
  rownames(genotype1_mock[!goodGenes(t(genotype1_mock)), ]),
  rownames(genotype1_sclero[!goodGenes(t(genotype1_sclero)), ]),
  rownames(genotype2_mock[!goodGenes(t(genotype2_mock)), ]),
  rownames(genotype2_sclero[!goodGenes(t(genotype2_sclero)), ])
))

# Filter the input matrix and subsets
input_in <- datExpr[!(rownames(datExpr) %in% bad_genes), ] %>% as.matrix() 
genotype1_mock_in <- genotype1_mock[!(rownames(genotype1_mock) %in% bad_genes), ] %>% as.matrix()
genotype1_sclero_in <- genotype1_sclero[!(rownames(genotype1_sclero) %in% bad_genes), ] %>% as.matrix()
genotype2_mock_in <- genotype2_mock[!(rownames(genotype2_mock) %in% bad_genes), ] %>% as.matrix()
genotype2_sclero_in <- genotype2_sclero[!(rownames(genotype2_sclero) %in% bad_genes), ] %>% as.matrix()

modules <- read.table("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/schil/TOM/moduleColorsTOM_ds2_mch03_genids.txt", 
    sep=" ", row.names=NULL) %>%
    column_to_rownames("GeneID")


#

clrs_in <- modules[!(rownames(modules) %in% bad_genes), 2]
names(clrs_in) <- rownames(modules[!(rownames(modules) %in% bad_genes), ])



mp <- modulePreservation(
  multiData = list(
    all = list(data = input_in),
    genotype1_Mock = list(data = genotype1_mock_in),
    genotype1_Sclero = list(data = genotype1_sclero_in),
    genotype2_Mock = list(data = genotype2_mock_in),
    genotype2_Sclero = list(data = genotype2_sclero_in)
  ),
  multiColor = list(all = clrs_in),   # Adjust module colors as needed
  dataIsExpr = TRUE,
  referenceNetworks = 1,  # Use all as reference
  nPermutations = 100,
  maxModuleSize = 5000, 
  randomSeed = 1,
  quickCor = 0,
  parallelCalculation = TRUE,
  verbose = 3,
  networkType = "signed hybrid",
  corFnc = "bicor"
)

# Calculate Z statistics for all subsets
zstat <- bind_rows(
  genotype1_Mock = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.genotype1_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "genotype1_Mock"),
  genotype1_Sclero = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.genotype1_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "genotype1_Sclero"),
  genotype2_Mock = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.genotype2_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "genotype2_Mock"),
  genotype2_Sclero = mp$preservation$Z$ref.all$inColumnsAlsoPresentIn.genotype2_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "genotype2_Sclero")
) 

write.table(zstat, "/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/schil/zstat/2024_12_01_zstat.tsv", sep="\t", row.names = FALSE)

png("WGCNA/documentation/pics/schil/presevation_zstats.png", width = 800, height = 400)

# Plot Z summary preservation statistics
ggplot(zstat, aes(x = moduleSize, y = Zsummary.pres,
                  shape = dataset, fill = moduleColor)) +
  geom_point(size = 3) +
  ylab("Preservation (Zsummary)") +
  scale_fill_manual(values = as.character(zstat$moduleColor),
                    guide = "none") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  geom_hline(yintercept = 10, col = "red") +
  geom_hline(yintercept = 2, col = "darkgreen") +
  theme_bw() +
  labs(title = "Module Preservation (Zsummary)", x = "Module Size")
dev.off()



# Calculate median rank statistics for all subsets
med_rank_stat <- bind_rows(
  genotype1_Mock = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.genotype1_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "genotype1_Mock"),
  genotype1_Sclero = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.genotype1_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "genotype1_Sclero"),
  genotype2_Mock = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.genotype2_Mock %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "genotype2_Mock"),
  genotype2_Sclero = mp$preservation$observed$ref.all$inColumnsAlsoPresentIn.genotype2_Sclero %>%
    rownames_to_column("moduleColor")%>%
    mutate(dataset = "genotype2_Sclero")
) 


# Plot median rank preservation statistics

png("WGCNA/documentation/pics/schil/presevation_median_rank.png", width = 800, height = 400)

ggplot(med_rank_stat, aes(x = moduleSize, y = medianRank.pres,
                          shape = dataset, fill = moduleColor)) +
  geom_point(size = 3) +
  scale_y_reverse() +
  ylab("Preservation (Median rank)") +
  scale_fill_manual(values = as.character(med_rank_stat$moduleColor),
                    guide = "none") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  theme_bw() +
  labs(title = "Module Preservation (Median Rank)", x = "Module Size")

dev.off()
