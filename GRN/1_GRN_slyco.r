# QDR RNAseq Solanum species
# Perform GRN analysis on S lycopersicoides
# Severin Einspanier
# 2025_01_10

rm(list=ls())
pacman::p_load(tidyverse, GENIE3, segmented, pheatmap)

setwd("")
options(stringsAsFactors = FALSE)

ids <- read.delim("slyco_curated_proteome_OG_pannzer_dedub_ids.txt", header=F) %>%
    mutate(gene=gsub(">", "", V1)) %>%
    mutate(gene=gsub("GeneExt~", "", gene))%>% 
    mutate(gene=gsub("mRNA_", "", gene))%>%  
    mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
    mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
    mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
    mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))%>%
    dplyr::select(gene)

modules <- read.delim("slyco_module_colors_TOM_filtered_genids.txt", row.names=1, sep=" ")
GO_TERM <- read.csv("filtered_protein_names_spen.csv", header=T, row.names=1)
wgcna_hub <- read.csv("2025_01_08_slyco_hub.csv", row.names=1)

# filter the gene expression:

datExpr <- read.csv("norm_counts_all_rlog.csv") %>%
  dplyr::filter(species=="S. lycopersicoides" & Geneid %in% ids$gene) %>% 
  dplyr::select(!species & !source) %>%
  pivot_wider(names_from=sample, values_from=normalized_count) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

TF <- read.delim("slyco_TFs.txt", header=F)

colnames(TF) = c("gene", "TF")

TF_ids <- TF %>%
    mutate(gene=gsub(">", "", gene)) %>%
    mutate(gene=gsub("GeneExt~", "", gene))%>% 
    mutate(gene=gsub("mRNA_", "", gene))%>%  
    mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
    mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
    mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
    mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))%>%
    dplyr::select(gene) %>%
    unique()

head(TF)

head(TF_ids)
dim(TF_ids)

expr_TFs <- datExpr %>% 
  dplyr::filter(!rownames(.) %in% TF_ids$gene) 
dim(expr_TFs)

# which TFs are not findable?

missing <- TF_ids %>% 
  filter(!gene %in% rownames(datExpr))
dim(missing)
head(missing)
# but seems like editing is not the issue:
# cant find them in the rlogs ... 

expressed_TFs <- TF_ids %>% 
  filter(gene %in% rownames(datExpr))
dim(missing)

weightMat <- GENIE3(datExpr, regulators=expressed_TFs$gene, nCores=30)

# Created "linked list" from a weighted adjacency_matrix

wam_linked_list <- getLinkList(weightMatrix = weightMat) %>%
    mutate_if(is.factor, as.character)

# calculate edge weight threshold

brkpnt_fun <- function(x, y, output_file = "breakpoint_plot.png") {
  if (length(x) > 10000) {
    set.seed(54321)
    x <- sample(x, size = 10000)
  }
  df <- data.frame(  # generate data frame
    vec = sort(x), num = seq_along(sort(x))
  )
  glm_ <- glm(vec ~ num, data = df)  # create linear model
  seg <- segmented(glm_, seg.Z =  ~ num, npsi = y)  # calculate breakpoint(s)
  brkpnt <- seg[["psi"]][, 2] |> round() # return breakpoint(s)
  cat(  # print breakpoint(s)
    "Breakpoints are \nx:\t",
    paste0(brkpnt[1:y], collapse = ", "),
    "\ny:\t",
    paste0(df$vec[brkpnt[1:y]], collapse = ", "),
    "\n"
  )
  
  # Generate plot
  p <- ggplot() +  
    geom_point(data = df, aes(x = num, y = vec)) +
    geom_vline(xintercept = as.numeric(brkpnt), color='red') +
    geom_hline(yintercept = as.numeric(df$vec[brkpnt]), color='forestgreen') +
    theme_classic()
  
  # Save plot to file
  ggsave(output_file, plot = p)
}

breakpont <- brkpnt_fun(wam_linked_list$weight, 7, paste0("slyco_breakpoints.png"))

# picked 0.00396004960324319
# filter dataframe with eigencentrality threshold
grn_edges_all <- wam_linked_list[wam_linked_list$weight > 0.00396004960324319,]

genes <- append(wam_linked_list$regulatoryGene, wam_linked_list$targetGene) %>% as.character() %>% unique() %>% sort()

genes_names <- as.data.frame(genes) %>% 
    left_join(GO_TERM, by=c("genes"="gene"))%>%
    mutate(protein=ifelse(is.na(desc), "unknown", desc)) %>%
    dplyr::select(!desc)
dim(genes_names)

#### define hubs 

grn_nodes_all_eigen <- graph_from_data_frame(grn_edges_all)  # generate igraph object
grn_all_eigen <- eigen_centrality(grn_nodes_all_eigen)  # calculate eigenvector centrality

grn_all_eigen <- data.frame(egnvctr = unname(grn_all_eigen[["vector"]]),  # generate dataframe of genes
                             gene = names(grn_all_eigen[["vector"]]))     # and eigenvector centrality

breakpont <- brkpnt_fun((grn_all_eigen$egnvctr), 5, paste0("slyco_breakpoints_hubs.png"))
# selected 0.0516089702523873
hubs <- grn_all_eigen[grn_all_eigen$egnvctr > 0.0516089702523873,]

genes_names_wth_hubs <- genes_names %>% 
    left_join(hubs, by=c("genes"="gene")) %>%
    mutate(hub_grn=ifelse(is.na(egnvctr), "non-hub", "hub")) %>%
    dplyr::select(!egnvctr) %>% 
    left_join(modules, by=c("genes"="GeneID")) %>%
    left_join(wgcna_hub, by=c("genes"="gene")) 
    # also append WGRC hubs 

dim(genes_names_wth_hubs)
head(genes_names_wth_hubs)
write_delim(genes_names_wth_hubs, file = "slyco_grn_network_nodes.txt", delim = "\t")

# output GRN network
write_delim(grn_edges_all, file = "slyco_grn_network_edges.txt", delim = "\t")