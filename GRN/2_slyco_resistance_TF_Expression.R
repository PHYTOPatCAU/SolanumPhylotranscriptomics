# QDR RNAseq Solanum species
# Check the resistance modules 
# Slycopersicoides
# Severin Einspanier

rm(list=ls())
library(tidyverse)
library(ComplexHeatmap)
setwd("/")
GRN_data <- read.delim("slyco_grn_network_nodes.txt")
TFs <- read.delim("GRN/data/slyco_TFs/slyco_TFs.txt", header=F)%>%
  mutate(gene=gsub(">", "", V1)) %>%
  mutate(gene=gsub("GeneExt~", "", gene))%>% 
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))%>%
  dplyr::select(gene)
###############################################################
# 2) Get TAI of the whole resistance-modules
###############################################################


# 2.1) Identify 
resistance_module_slyco <- c("magenta","turquoise", "green", "red")

slyco_resistance_module <- GRN_data %>%
  filter(ModuleColor %in% resistance_module_slyco) %>%
  filter(hub=="hub" & hub_grn=="hub" & genes %in% TFs$gene) %>% 
  filter(protein !="Calcium-dependent protein kinase")

#genes:

slyco_resistance_module$genes

# from those we should make a directed subnetwork

# plot their regulation:

DEGs_res <- read.csv("DeSeq/data/DeSeq_OUT/combined_RES_inf_mock.csv")%>% 
  filter(species=="S. lycopersicoides") %>% 
  mutate(res=log2FoldChange) %>% 
  dplyr::select(GeneID, res)

DEGs <- read.csv("DeSeq/data/DeSeq_OUT/combined_SUS_inf_mock.csv") %>% 
  filter(species=="S. lycopersicoides") %>% 
  mutate(sus=log2FoldChange) %>% 
  dplyr::select(GeneID, sus)%>% 
  left_join(DEGs_res) %>% 
  filter(GeneID %in% slyco_resistance_module$genes) %>% 
  column_to_rownames("GeneID")

# write.csv(DEGs, "GRN/data/slyco_resistance_module_DEGs.csv")

DEGs_res_ann <- read.csv("DeSeq/data/DeSeq_OUT/combined_RES_inf_mock.csv")%>% 
  filter(species=="S. lycopersicoides") %>% 
  mutate(res=log2FoldChange,p_res=padj) %>% 
  dplyr::select(GeneID, res, p_res)

# consider adding information on: DEG res-sus to the heatmap too 

susresDEG <- read.csv("DeSeq/data/DeSeq_OUT/combined_INF_res_sus_DEGs.csv") %>% 
  filter(species =="S. lycopersicoides")

annotation <- read.csv("DeSeq/data/DeSeq_OUT/combined_SUS_inf_mock.csv") %>% 
  filter(species=="S. lycopersicoides") %>% 
  mutate(sus=log2FoldChange,p_sus=padj) %>% 
  dplyr::select(GeneID, sus, p_sus)%>% 
  left_join(DEGs_res_ann) %>%   
  right_join(slyco_resistance_module, by=c("GeneID"="genes"))  %>% 
  mutate(DEG_res = ifelse((abs(res) > 1 & p_res < 0.05), "DEG", "noDEG"),
         DEG_sus = if_else((abs(sus) > 1 & p_sus < 0.05 ), "DEG","noDEG"), 
         DEG_susVres = if_else(GeneID %in%susresDEG$GeneID,"DEG", "noDEG")) %>% 
  column_to_rownames("GeneID") %>% 
  dplyr::select(ModuleColor, DEG_sus, DEG_res, DEG_susVres)


# Define annotation colors
annot_colors <- list(ModuleColor = c(blue = "blue", red = "red", pink = "pink", 
                                     turquoise="turquoise", magenta="magenta",
                                     green="green"),
                     DEG_sus = c(DEG = "#fc7481", noDEG = "#7a7f82"),
                     DEG_res = c(DEG = "#bbef77", noDEG = "#7a7f82"),
                     DEG_susVres= c(DEG = "#bdb9e5", noDEG = "#7a7f82"))

# Create heatmap annotation
ha <- rowAnnotation(df = annotation, col = annot_colors, 
                    gp= gpar(col="black", lwd = 2), 
                    gap=unit(1, "mm"))
col_fun <- circlize::colorRamp2(c(-6, 0, 6), c("darkblue", "white", "darkred"))

# Plot heatmap
png(paste0("figures/fig_7/", Sys.Date(), "slyco_resistance_TFs.png"), 
    width = 16.5, height = 16.5, unit = "cm", res = 960)
Heatmap(as.matrix(DEGs), name = "log2FoldChange", 
        left_annotation = ha, 
        #annotation = lgd,
        cluster_columns = FALSE, 
        show_row_names = TRUE,
        width = ncol(DEGs) * unit(5, "mm"), 
        height = nrow(DEGs) * unit(5, "mm"),
        cluster_rows = TRUE, 
        row_dend_width = unit(1, "cm"),
        rect_gp = gpar(col = "black", lwd = 2),
        border = "black",
        na_col = "grey",
        col=col_fun,
        show_column_names = TRUE,
        heatmap_legend_param = list(
          at=c(-6, -3,0,3,6),
          title = "log2FoldChange",
          legend_height = unit(6, "cm"), 
          direction = "horizontal")
)
dev.off()

########
# write csv

annotation_print <- annotation %>% 
  rownames_to_column("GeneID") 

write.csv(annotation_print,"GRN/data/slyco_resistance_TFs_annotation.csv")
