# Overlap analysis on the key-TFs

#key_tfs <- c("Sopen09g001470", "Sopen10g006210")
# conserved?

# which OG? 

spenID2OG <- read.csv("OGs/data/genid2goid_spen.csv", row.names=1) 
slycoID2OG <- read.csv("OGs/data/genid2goid_slyco.csv", row.names=1) 
schilID2OG <- read.csv("OGs/data/genid2goid_schil.csv", row.names=1)
shabroID2OG <- read.csv("OGs/data/genid2goid_shabro.csv", row.names=1)
spimpID2OG <- read.csv("OGs/data/genid2goid_spimp.csv", row.names=1)

# get general OG-info 

OG_info <- read_tsv("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/orthofinder/NOV_2024_ITAG_PANNZER/Orthogroups/Orthogroups.GeneCount.tsv"
) %>% 
  filter(Orthogroup %in% spenID2OG$OG)

# OG0015518 MADS box transcription factor

# get heatmap for the OG-TFs and see whether the regulation after infection is 
# always same?


# get list of resistance TFs:
tfs_slyco <- read.csv("GRN/data/slyco_resistance_TFs_annotation.csv", row.names = 1) %>% 
  left_join(slycoID2OG, by=c("GeneID"="GeneID")) 

tfs_spen <- read.csv("GRN/data/spen_resistance_TFs_annotation.csv", row.names = 1) %>% 
  left_join(spenID2OG, by=c("GeneID"="GeneID"))

tfs_merged <- tfs_slyco %>% 
  bind_rows(tfs_spen) 
summary(as.factor(tfs_merged$OG))

# get only the overlap 
overlap <- tfs_merged %>%
  unique() %>% 
  group_by(OG) %>%
  filter(n()>1) 
# 4 OGs overlap

# only Sopen05g003630 / Solyd05g052890 / OG0005445 relevantly overlapping


OG_info <- read_tsv("C:/Users/suaph281/Desktop/nesh_local/LDT_rnaseq/orthofinder/NOV_2024_ITAG_PANNZER/Orthogroups/Orthogroups.GeneCount.tsv"
) %>% 
  filter(Orthogroup == "OG0005445")
# very interesting

# Draw heatmap for all species here. 
# for the overlappers 
resistance_OGs <- c("OG0005445", "OG0001338")

# import expressiondata of all others 


schil_genes <- read.csv("DeSeq/data/DeSeq_OUT/combined_RES_inf_mock.csv")%>% 
  filter(species=="S. chilense") %>% 
  mutate(res=log2FoldChange) %>% 
  dplyr::select(GeneID, res)

schil_genes_all <- read.csv("DeSeq/data/DeSeq_OUT/combined_SUS_inf_mock.csv") %>% 
  filter(species=="S. chilense") %>% 
  mutate(sus=log2FoldChange) %>% 
  dplyr::select(GeneID, sus)%>% 
  left_join(schil_genes) %>% 
  left_join(schilID2OG, by=c("GeneID"="GeneID")) %>% 
  filter(OG %in% resistance_OGs) 
  
# shabro

shabro_genes <- read.csv("DeSeq/data/DeSeq_OUT/combined_RES_inf_mock.csv")%>% 
  filter(species=="S. habrochaites") %>% 
  mutate(res=log2FoldChange) %>% 
  dplyr::select(GeneID, res)

shabro_genes_all <- read.csv("DeSeq/data/DeSeq_OUT/combined_SUS_inf_mock.csv") %>%
  filter(species=="S. habrochaites") %>% 
  mutate(sus=log2FoldChange) %>% 
  dplyr::select(GeneID, sus)%>% 
  left_join(shabro_genes) %>% 
  left_join(shabroID2OG, by=c("GeneID"="GeneID")) %>% 
  filter(OG %in% resistance_OGs)

# spimp

spimp_genes <- read.csv("DeSeq/data/DeSeq_OUT/combined_RES_inf_mock.csv")%>% 
  filter(species=="S. pimpinellifolium") %>% 
  mutate(res=log2FoldChange) %>% 
  dplyr::select(GeneID, res)

spimp_genes_all <- read.csv("DeSeq/data/DeSeq_OUT/combined_SUS_inf_mock.csv") %>%
  filter(species=="S. pimpinellifolium") %>% 
  mutate(sus=log2FoldChange) %>% 
  dplyr::select(GeneID, sus)%>% 
  left_join(spimp_genes) %>% 
  left_join(spimpID2OG, by=c("GeneID"="GeneID")) %>% 
  filter(OG %in% resistance_OGs)

# slyco

slyco_genes <- read.csv("DeSeq/data/DeSeq_OUT/combined_RES_inf_mock.csv")%>% 
  filter(species=="S. lycopersicoides") %>% 
  mutate(res=log2FoldChange) %>% 
  dplyr::select(GeneID, res)

slyco_genes_all <- read.csv("DeSeq/data/DeSeq_OUT/combined_SUS_inf_mock.csv") %>%
  filter(species=="S. lycopersicoides") %>% 
  mutate(sus=log2FoldChange) %>% 
  dplyr::select(GeneID, sus)%>% 
  left_join(slyco_genes) %>% 
  left_join(slycoID2OG, by=c("GeneID"="GeneID")) %>% 
  filter(OG %in% resistance_OGs)

# spen

spen_genes <- read.csv("DeSeq/data/DeSeq_OUT/combined_RES_inf_mock.csv")%>% 
  filter(species=="S. pennellii") %>% 
  mutate(res=log2FoldChange) %>% 
  dplyr::select(GeneID, res)

spen_genes_all <- read.csv("DeSeq/data/DeSeq_OUT/combined_SUS_inf_mock.csv") %>%
  filter(species=="S. pennellii") %>% 
  mutate(sus=log2FoldChange) %>% 
  dplyr::select(GeneID, sus)%>% 
  left_join(spen_genes) %>% 
  left_join(spenID2OG, by=c("GeneID"="GeneID")) %>% 
  filter(OG %in% resistance_OGs) 

merged <- spen_genes_all %>% 
  bind_rows(slyco_genes_all) %>% 
  bind_rows(shabro_genes_all) %>% 
  bind_rows(spimp_genes_all) %>% 
  bind_rows(schil_genes_all) %>% 
  group_by(OG, species) %>%
  reframe(sus=mean(sus),
            res=mean(res)) %>% 
  unique()
OGs_1 <- merged %>% 
  filter(OG=="OG0005445") %>% 
  select(!OG) %>% 
  column_to_rownames("species") 
OGs_2 <- merged %>% 
  filter(OG=="OG0001338") %>% 
  select(!OG) %>% 
  column_to_rownames("species")


pheatmap::pheatmap(OGs_1, cluster_rows = F, cluster_cols = F, main="OG0005445")
pheatmap::pheatmap(OGs_2, cluster_rows = F, cluster_cols = F, main="OG0001338")


### complex graph

# Install if needed
# install.packages("ComplexHeatmap")
# install.packages("circlize")  # for color mapping
# install.packages("dplyr")
pacman::p_load(ComplexHeatmap, circlize, dplyr)
merged <- spen_genes_all %>% 
  bind_rows(slyco_genes_all) %>% 
  bind_rows(shabro_genes_all) %>% 
  bind_rows(spimp_genes_all) %>% 
  bind_rows(schil_genes_all)
# Extract just the numeric columns for the heatmap
mat <- as.matrix(merged[, c("sus", "res")])

# Assign rownames to identify each row by GeneID
rownames(mat) <- merged$GeneID

# A simple annotation data frame
row_anno_merged <- data.frame(
  OG = merged$OG,
  Species = merged$species
)
rownames(row_anno_merged) <- merged$GeneID

# Create the ComplexHeatmap row annotation
row_anno <- rowAnnotation(
  df = row_anno_merged,    # specify it's a data frame
  col = list(
    OG = c("OG0001338" = "gold", "OG0005445" = "skyblue"),
    Species = c("S. pennellii" = "#1b9e77", 
                "S. lycopersicoides" = "#d95f02",
                "S. habrochaites" = "#7570b3",
                "S. pimpinellifolium" = "#e7298a",
                "S. chilense" = "#66a61e")
  )
)

# Use circlize::colorRamp2 to define breaks
# For example, if your log-fold changes range approx. -2 to +4:
col_fun <- colorRamp2(
  breaks = c(-2, 0, 4),          # adapt to your min/mid/max
  colors = c("blue", "white", "red")
)
Heatmap(
  matrix          = mat,
  name            = "Value",          # legend title
  col             = col_fun,          # the color function we defined
  row_split       = merged$OG,            # split rows by OG
  cluster_columns = FALSE,            # usually only 2 columns (sus vs. res), so often no real cluster there
  show_row_names  = TRUE,             # show your GeneIDs
  row_names_gp    = gpar(fontsize=8), # text size for row labels
  top_annotation  = HeatmapAnnotation( # optional: label the columns
    Condition = c("sus","res"),
    col = list(Condition = c("sus"="lightgreen", "res"="orange"))
  ),
  left_annotation = row_anno          # species + OG annotation on the left
)

