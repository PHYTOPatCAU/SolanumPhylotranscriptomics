###################################
# Vizualization of Resistance-TFs 
# of S. lycopersicoides in all other 
# Solanum species
# 
# Severin Einspanier
# 2025_01_13
###################################

# Load packages
pacman::p_load(tidyverse, ComplexHeatmap, circlize)

# Clean workspace
rm(list = ls())

# Set working directory
setwd("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/")

#==========================
# 1) LOAD ORTHOGROUP DATA
#==========================
spenID2OG   <- read.csv("OGs/data/genid2goid_spen.csv",   row.names=1)
slycoID2OG  <- read.csv("OGs/data/genid2goid_slyco.csv",  row.names=1)
schilID2OG  <- read.csv("OGs/data/genid2goid_schil.csv",  row.names=1)
shabroID2OG <- read.csv("OGs/data/genid2goid_shabro.csv", row.names=1)
spimpID2OG  <- read.csv("OGs/data/genid2goid_spimp.csv",  row.names=1)

# define the genes. Upregulated in the resistant genotype
slyco_genes <- c("Solyd03g050610", "Solyd02g064330", "Solyd09g071070")
# Which Orthogroups do these genes belong to?

OG_of_interest <- slycoID2OG %>% 
  filter(GeneID %in% slyco_genes)

OG_info <- read_tsv("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/orthofinder/NOV_2024_ITAG_PANNZER/Orthogroups/Orthogroups.GeneCount.tsv"
) %>% 
  filter(Orthogroup %in% OG_of_interest$OG)
# found in all species

#==========================
# 2) RESISTANCE MODULES
#==========================
res_modules <- read.table("WGCNA/data/resistance_modules_all_species.txt", header = TRUE) %>% 
  mutate(is_resistance = 1)

res_modules_spimp <- res_modules %>% filter(Species == "spimp")
res_modules_schil <- res_modules %>% filter(Species == "schil")
res_modules_shabro <- res_modules %>% filter(Species == "shabro")
res_modules_slyco <- res_modules %>% filter(Species == "slyco")
res_modules_spen  <- res_modules %>% filter(Species == "spen")


#===========================
# 3) LOAD EXPRESSION DATA
#===========================
# We'll define a little function that:
#   (a) reads combined_RES_inf_mock.csv
#   (b) reads combined_SUS_inf_mock.csv
#   (c) merges them, merges OG, filters OG_of_interest
#   (d) returns a combined data.frame
# so we don’t repeat similar code for each species.

load_expression_data <- function(species_label, og_map) {
  
  # Resist
  df_res <- read.csv("DeSeq/data/DeSeq_OUT/combined_RES_inf_mock.csv") %>%
    filter(species == species_label) %>%
    mutate(res = log2FoldChange) %>%
    dplyr::select(GeneID, res)
  
  # Susceptible
  df_sus <- read.csv("DeSeq/data/DeSeq_OUT/combined_SUS_inf_mock.csv") %>%
    filter(species == species_label) %>%
    mutate(sus = log2FoldChange) %>%
    dplyr::select(GeneID, sus)
  
  # Merge them
  df_combined <- df_sus %>%
    left_join(df_res, by = "GeneID") %>%
    left_join(og_map, by = "GeneID") %>%
    filter(OG %in% OG_of_interest$OG) %>%
    mutate(species = species_label)
  
  return(df_combined)
}

# Apply function to each species
schil_genes_all <- load_expression_data("S. chilense",        schilID2OG) %>%
  # for consistent naming
  mutate(GeneID = str_replace(GeneID, "chilense_Hirise_Scaffold_", "Schil"),
         GeneID = str_replace(GeneID, "__1_contigs__length_", "c"),
         GeneID = str_replace(GeneID, "__2_contigs__length_", "c"),
         GeneID = str_replace(GeneID, "__3_contigs__length_", "c"))

shabro_genes_all <- load_expression_data("S. habrochaites",    shabroID2OG)
spimp_genes_all  <- load_expression_data("S. pimpinellifolium", spimpID2OG)
slyco_genes_all  <- load_expression_data("S. lycopersicoides",  slycoID2OG)
spen_genes_all   <- load_expression_data("S. pennellii",        spenID2OG)

# Merge all expression data
merged <- bind_rows(
  spen_genes_all,
  slyco_genes_all,
  shabro_genes_all,
  spimp_genes_all,
  schil_genes_all
) %>%
  group_by(species)


#==============================
# 4) WHICH OGs TO VISUALIZE?
#==============================

OGs <- c("OG0000116", "OG0003738", "OG0009560")


# Filter your orthogroup mappings to these OGs only
slycoID2OG_f  <- slycoID2OG  %>% filter(OG %in% OGs)
schilID2OG_f  <- schilID2OG  %>% filter(OG %in% OGs)
shabroID2OG_f <- shabroID2OG %>% filter(OG %in% OGs)
spimpID2OG_f  <- spimpID2OG  %>% filter(OG %in% OGs)
spenID2OG_f   <- spenID2OG   %>% filter(OG %in% OGs)


#=============================
# 5) MODULE MEMBERSHIP
#=============================
# We'll read the gene->moduleColor file for each species
# and tag each gene as "Resistance" or "Non-Resistance"
# depending on membership in the modules from res_modules_*

schil_module <- read.table("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/WGCNA/schil/moduleColorsTOM_ds2_mch25_genids.txt",
                           header = TRUE) %>%
  right_join(schilID2OG_f, by = "GeneID") %>%
  mutate(WGCNA = if_else(ModuleColor %in% res_modules_schil$Module,
                         "Resistance","Non-Resistance")) %>%
  mutate(GeneID = str_replace(GeneID, "chilense_Hirise_Scaffold_", "Schil"),
         GeneID = str_replace(GeneID, "__1_contigs__length_", "c"),
         GeneID = str_replace(GeneID, "__2_contigs__length_", "c"),
         GeneID = str_replace(GeneID, "__3_contigs__length_", "c"))

slyco_module <- read.table("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/WGCNA/slyco/slyco_module_colors_TOM_filtered_genids.txt",
                           header = TRUE) %>%
  right_join(slycoID2OG_f, by = "GeneID") %>%
  mutate(WGCNA = if_else(ModuleColor %in% res_modules_slyco$Module,
                         "Resistance", "Non-Resistance"))

shabro_module <- read.table("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/WGCNA/shabro/TOMfinal_module_colors_geneids.txt",
                            header = TRUE) %>%
  right_join(shabroID2OG_f, by = "GeneID") %>%
  mutate(WGCNA = if_else(ModuleColor %in% res_modules_shabro$Module,
                         "Resistance", "Non-Resistance"))

spimp_module <- read.table("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/WGCNA/spimp/TOM_final_moduleColors_genids.txt",
                           header = TRUE) %>%
  right_join(spimpID2OG_f, by = "GeneID") %>%
  mutate(WGCNA = if_else(ModuleColor %in% res_modules_spimp$Module,
                         "Resistance", "Non-Resistance"))

spen_module <- read.table("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/WGCNA/spen/module_colors_TOM_filtered_genids.txt",
                          header = TRUE) %>%
  right_join(spenID2OG_f, by = "GeneID") %>%
  mutate(WGCNA = if_else(ModuleColor %in% res_modules_spen$Module,
                         "Resistance", "Non-Resistance"))

# Combine all module membership info
all_modules <- bind_rows(
  schil_module,
  slyco_module,
  shabro_module,
  spimp_module,
  spen_module
)

#=========================
# 6) CREATE HEATMAP DATA
#=========================
# For each OG, we build:
#   - annotation matrix (WGCNA + species)
#   - expression matrix (res + sus columns)

# Utility function to produce annotation + expression matrix
make_anno_and_mat <- function(og_id) {
  
  ann_df <- merged %>%
    filter(OG == og_id) %>%
    left_join(all_modules, by = c("species", "GeneID", "OG")) %>%
    unique() %>% 
    column_to_rownames("GeneID") %>%
    dplyr::select(WGCNA, species)
  
  mat_df <- merged %>%
    filter(OG == og_id) %>%
    unique() %>% 
    column_to_rownames("GeneID") %>%
    dplyr::select(!c(OG, species))  # keep 'res' and 'sus' columns
  
  return(list(
    ann = as.matrix(ann_df),
    mat = as.matrix(mat_df)
  ))
}

# Generate data for each OG
og1_data <- make_anno_and_mat(OGs[1])
og2_data <- make_anno_and_mat(OGs[2])
og3_data <- make_anno_and_mat(OGs[3])

ann_1 <- og1_data$ann
mat1 <- og1_data$mat

ann_2 <- og2_data$ann
mat2 <- og2_data$mat

ann_3 <- og3_data$ann
mat3 <- og3_data$mat

#===============================
# 7) HEATMAP COLORS & ANNOTS
#===============================
# Define color gradients for each OG
col_fun1 <- circlize::colorRamp2(c(-2,0, 8),   c("darkblue", "white", "darkred"))
col_fun2 <- circlize::colorRamp2(c(-2,0, 8),   c("darkblue", "white", "darkred"))
col_fun3 <- circlize::colorRamp2(c(-2,0, 8),   c("darkblue", "white", "darkred"))

# Define annotation colors
annot_colors <- list(
  "WGCNA"   = c("Resistance" = "#b6a710", "Non-Resistance" = "#eae1b5"),
  "species" = c("S. chilense"         = "#2B2B2B",
                "S. habrochaites"     = "#820024",
                "S. lycopersicoides"  = "#00B8A3",
                "S. pennellii"        = "#F2EFDC",
                "S. pimpinellifolium" = "#FF9607")
)

ha1 <- rowAnnotation(
  df = ann_1, 
  col = annot_colors, 
  show_annotation_name = FALSE,
  show_legend = FALSE, 
  gp= gpar(col="black", lwd = 2)  # disable automatic annotation legends
)
ha2 <- rowAnnotation(
  df = ann_2, 
  col = annot_colors, 
  show_annotation_name = FALSE,
  show_legend = FALSE, 
  gp= gpar(col="black", lwd = 2)
)
ha3 <- rowAnnotation(
  df = ann_3, 
  col = annot_colors, 
  show_annotation_name = T,
  show_legend = FALSE, 
  gp= gpar(col="black", lwd = 2)
)

#=========================
# 8) CONSTRUCT HEATMAPS
#=========================
h1 <- Heatmap(
  mat1,
  left_annotation   = ha1,
  width             = ncol(mat1) * unit(5, "mm"), 
  height            = nrow(mat1) * unit(5, "mm"),
  cluster_rows      = TRUE,
  cluster_columns   = FALSE,
  name              = "log2FC \nOG0000116",
  row_title         = "OG0000116",
  col               = col_fun1,
  row_dend_width    = unit(1, "cm"),
  rect_gp           = gpar(col = "black", lwd = 2),
  border            = "black",
  na_col            = "grey",
  show_heatmap_legend = FALSE  # turn off color legend for this heatmap
)

h2 <- Heatmap(
  mat2,
  left_annotation   = ha2,
  width             = ncol(mat2) * unit(5, "mm"),
  height            = nrow(mat2) * unit(5, "mm"),
  cluster_rows      = TRUE,
  cluster_columns   = FALSE,
  row_title         = "OG0003738",
  name              = "log2FC \nOG0003738",
  col               = col_fun2,
  row_dend_width    = unit(1, "cm"),
  rect_gp           = gpar(col = "black", lwd = 2),
  border            = "black",
  na_col            = "grey",
  show_heatmap_legend = FALSE
)

h3 <- Heatmap(
  mat3,
  left_annotation   = ha3,
  width             = ncol(mat3) * unit(5, "mm"),
  height            = nrow(mat3) * unit(5, "mm"),
  cluster_rows      = TRUE,
  cluster_columns   = FALSE,
  row_title         = "OG0009560",
  name              = "log2FC \nOG0009560",
  col               = col_fun3,
  row_dend_width    = unit(1, "cm"),
  rect_gp           = gpar(col = "black", lwd = 2),
  border            = "black",
  na_col            = "grey",
  show_heatmap_legend = FALSE
)

# Combine heatmaps vertically
hlist <- h1 %v% h2 %v% h3

#============================
# 9) MANUAL LEGENDS
#============================
# Heatmap color legends (vertical)
lgd1 <- Legend(title = "log2FC \nOG0000116", col_fun = col_fun1, at=c(-2, 0, 2,4,6, 8))
lgd2 <- Legend(title = "log2FC \nOG0003738", col_fun = col_fun2, at=c(-2, 0, 2,4,6, 8))
lgd3 <- Legend(title = "log2FC \nOG0009560", col_fun = col_fun3, at=c(-2, 0, 2,4,6, 8))

#combined_legend <- packLegend(lgd1, lgd2, lgd3, direction = "vertical")
combined_legend <- packLegend(lgd1, direction = "vertical")
# Annotation legends (horizontal)
species_legend <- Legend(
  labels    = names(annot_colors$species),
  title     = "Species",
  legend_gp = gpar(fill = annot_colors$species)
)

WGCNA_legend <- Legend(
  labels    = names(annot_colors$WGCNA),
  title     = "WGCNA",
  legend_gp = gpar(fill = annot_colors$WGCNA)
)

annotation_legends <- packLegend(
  species_legend, 
  WGCNA_legend, 
  direction = "horizontal"
)

#============================
# 10) DRAW & SAVE PLOT
#============================
png(
  filename = paste0(
    "C:/Users/suaph281/Nextcloud/ResiDEvo/sequencing/figures/fig_7/",
    Sys.Date(),
    "_slyco_resis_TFs_expression_across.png"
  ),
  width  = 16.5, 
  height = 27, 
  units  = "cm", 
  res    = 900
)

draw(
  hlist, 
  annotation_legend_list = annotation_legends, 
  heatmap_legend_list    = combined_legend,
  gap                    = unit(10, "mm"),
  heatmap_legend_side    = "right",
  annotation_legend_side = "top"
)

dev.off()
