# export networks 
# get hub genes
# spen 
# 2024_12_17 
# 
rm(list=ls())
library(tidyverse)
library(WGCNA)
library(segmented)

setwd("/gxfs_home/cau/suaph281/2024_solanum_ldt_rnaseq/")
options(stringsAsFactors = FALSE)

# take /gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/PROTEOME/spen_curated_proteome_OG_pannzer_dedub_ids.txt
# remove '>'
# remove 'GeneExt~'
# change 't.' to 'g.' 
# remove '.p1-9'

ids <- read.delim("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/PROTEOME/spen_curated_proteome_OG_pannzer_dedub_ids.txt",
    header=F) %>%
    mutate(gene=gsub(">", "", V1)) %>%
    mutate(gene=gsub("GeneExt~", "", gene))%>% 
    mutate(gene=gsub("mRNA_", "", gene))%>%  
    mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
    mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
    mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
    mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))%>%
    dplyr::select(gene)


# filter the gene expression:

datExpr <- read.csv("DeSeq/data/norm_counts_all_rlog.csv") %>%
  dplyr::filter(species=="S. pennellii" & Geneid %in% ids$gene) %>% 
  dplyr::select(!species & !source) %>% 
  pivot_wider(names_from = Geneid, values_from = normalized_count) %>%
  column_to_rownames("sample")

# load TOM

load("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/TOM/TOM_filtered-block.1.RData")
moduleColors = read.table("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/TOM/module_colors_TOM_filtered.txt", header = FALSE)
moduleColors  = as.character(moduleColors$V2)
head(moduleColors)



# Export to Cytoscape

out_network <- exportNetworkToCytoscape(
  TOM, 
   edgeFile = NULL,
   nodeFile = NULL,
   weighted = TRUE,
   threshold = 0.01,
   nodeNames = names(datExpr),
   altNodeNames = NULL,
   nodeAttr = moduleColors,
   includeColNames = TRUE)

summary(out_network$edgeData$weight)
unique(colors)

# Function to get the top N edges
get_top_edges <- function(edgeData, n) {
  edgeData %>%
    arrange(desc(weight)) %>%
    head(n) %>%
    tail(1) %>%
    dplyr::select(weight)
}

# Get top edges
top_1k_edges <- get_top_edges(out_network$edgeData, 1000)
top_10k_edges <- get_top_edges(out_network$edgeData, 10000)
top_100k_edges <- get_top_edges(out_network$edgeData, 100000)
top_1M_edges <- get_top_edges(out_network$edgeData, 1000000)

# Export networks
out_network_threshold_1k <- exportNetworkToCytoscape(
  TOM, 
  edgeFile = paste0("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE/TOM_1k_edge_", Sys.Date(), "_filtered.tsv"),
  nodeFile = paste0("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE/TOM_1k_node_", Sys.Date(), "_filtered.tsv"),
  weighted = TRUE,
  threshold = top_1k_edges$weight,
  nodeNames = names(datExpr),
  #altNodeNames = modGenes,
  nodeAttr = moduleColors,
  includeColNames = TRUE
)

out_network_threshold_10k <- exportNetworkToCytoscape(
  TOM, 
  edgeFile = paste0("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE/TOM_10k_edge_", Sys.Date(), "_filtered.tsv"),
  nodeFile = paste0("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE/TOM_10k_node_", Sys.Date(), "_filtered.tsv"),
  weighted = TRUE,
  threshold = top_10k_edges$weight,
  nodeNames = names(datExpr),
  #altNodeNames = modGenes,
  nodeAttr = moduleColors,
  includeColNames = TRUE
)

out_network_threshold_100k <- exportNetworkToCytoscape(
  TOM, 
  edgeFile = paste0("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE/TOM_100k_edge_", Sys.Date(), "_filtered.tsv"),
  nodeFile = paste0("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE/TOM_100k_node_", Sys.Date(), "_filtered.tsv"),
  weighted = TRUE,
  threshold = top_100k_edges$weight,
  nodeNames = names(datExpr),
  #altNodeNames = modGenes,
  nodeAttr = moduleColors,
  includeColNames = TRUE
)

out_network_threshold_1M <- exportNetworkToCytoscape(
  TOM, 
  edgeFile = paste0("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE/TOM_1M_edge_", Sys.Date(), "_filtered.tsv"),
  nodeFile = paste0("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE/TOM_1M_node_", Sys.Date(), "_filtered.tsv"),
  weighted = TRUE,
  threshold = top_1M_edges$weight,
  nodeNames = names(datExpr),
  #altNodeNames = modGenes,
  nodeAttr = moduleColors,
  includeColNames = TRUE
)

exportNetworkToCytoscape(
  TOM, 
   edgeFile = paste0("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE//TOM_full_edge_", Sys.Date(), "_filtered.tsv"),
   nodeFile = paste0("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE/TOM_full_node_", Sys.Date(), "_filtered.tsv"),
   weighted = TRUE,
   nodeNames = names(datExpr),
   threshold =0,
   #altNodeNames = modGenes,
   nodeAttr = moduleColors,
   includeColNames = TRUE)
   

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


breakpont <- brkpnt_fun(out_network$edgeData$weight,5,paste0("WGCNA/documentation/pics/spen/", Sys.Date(), "_spen_net_edges_filtered.tsv.png"))

# I'll select 0.0403409541741956

out_network_threshold <- exportNetworkToCytoscape(
  TOM, 
   edgeFile = paste0("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE/TOM_INFL_filtered_edge_", Sys.Date(), "_filtered.tsv"),
   nodeFile = paste0("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE/TOM_INFL_filtered_node_", Sys.Date(), "_filtered.tsv"),
   weighted = TRUE,
   threshold = 0.0403409541741956,
   nodeNames = names(datExpr),
   #altNodeNames = modGenes,
   nodeAttr = moduleColors,
   includeColNames = TRUE)


# Get Hub genes 

net_edges <- read.table("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE//TOM_INFL_filtered_edge_2024-12-17_filtered.tsv", header=T)

net_nodes <- read.table("/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/CYTOSCAPE//TOM_INFL_filtered_node_2024-12-17_filtered.tsv", header=T, sep = "\t") 


module_df <- net_nodes %>% 
  dplyr::select(nodeName, nodeAttr.nodesPresent...) %>%
  rename(gene_id = nodeName, colors = nodeAttr.nodesPresent...) 

# 
summary(as.factor(module_df$colors))

xtrxt_net <- function(x){  # put in desired module in quotes
  igraph::graph_from_data_frame(net_edges[((net_edges$fromNode %in% module_df$gene_id[module_df$colors == x]) &
                                     (net_edges$toNode %in% module_df$gene_id[module_df$colors == x])),],
                        vertices = net_nodes[net_nodes$nodeAttr.nodesPresent... == x,])
}

brkpnt_fun <- function(x, y) {
  # Ensure the input is numeric
  if (!is.numeric(x)) {
    stop("Input to brkpnt_fun must be numeric.")
  }
  
  df <- data.frame(vec = sort(x), num = seq_along(sort(x)))  # generate data frame
  glm_ <- glm(vec ~ num, data = df)  # create linear model
  
  # Try to calculate breakpoints and handle potential errors
  while (y > 0) {
    seg <- tryCatch({
      segmented::segmented(glm_, seg.Z = ~ num, npsi = y)  # calculate breakpoints
    }, warning = function(w) {
      warning(paste("Breakpoint estimation failed with", y, "breakpoints. Trying with", y - 1, "breakpoints."))
      return(NULL)
    }, error = function(e) {
      warning(paste("Breakpoint estimation failed with", y, "breakpoints. Trying with", y - 1, "breakpoints."))
      return(NULL)
    })
    
    if (!is.null(seg) && is.list(seg)) {
      brkpnt <- seg[["psi"]][, 2] %>% round()  # return the 2nd breakpoint
      print(paste0("Breakpoints are ", paste(brkpnt[1:y], collapse = ", ")))  # prints the breakpoints
      
      return(brkpnt[y])
    }
    
    y <- y - 1
  }
  
  # If all attempts fail, calculate the threshold for the top 5% biggest values
  warning("All attempts to estimate breakpoints failed. Returning top 5% threshold.")
  top_5_percent_index <- ceiling(0.95 * nrow(df))
  threshold <- df$num[top_5_percent_index]
  return(threshold)
}

process_module <- function(module) {
  # Filter for the current module
  net_col <- xtrxt_net(module)

  # Calculate eigencentrality
  egnvctr_col <- igraph::eigen_centrality(net_col)
  egnvctr_col <- data.frame(egnvctr = unname(egnvctr_col[["vector"]]),
                            gene = names(egnvctr_col[["vector"]]))

  # Sort and find breakpoint
  sorted <- egnvctr_col %>%
    dplyr::select(egnvctr) %>%
    arrange((egnvctr)) %>%
    rownames_to_column("position")
  brekpont <- brkpnt_fun(sorted$egnvctr,5)

  # Create hub variable and name the data frame
  hub_col <- egnvctr_col %>%
    arrange(egnvctr) %>%
    dplyr::select(gene, egnvctr) %>%
    rownames_to_column("position") %>%
    mutate(hub = ifelse(as.numeric(rownames(.)) > brekpont, "hub", "non-hub"))

  # Add a module column to the data frame
  hub_col$module <- module

  return(hub_col)
}


# Get unique modules
modules <- unique(module_df$colors)


combined_results <- modules %>%
  map_dfr(process_module)


# Plot the results
png(file = "WGCNA/documentation/pics/spen/filtered_wgcna/hub_genes_new.png", width = 10000, height = 6000, 
    res=600);

(p_1 <- combined_results %>%
  ggplot(aes(x = as.numeric(position), y = egnvctr, color = hub)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("black", "red"))+
    facet_wrap(~module, scales = "free_x") 
)

dev.off()

# manual adjustment:

combined_backup <- combined_results

combined_results$hub <- ifelse(combined_results$module == "pink" & as.numeric(combined_results$position) >= 262, "hub", combined_results$hub)
#combined_results$hub <- ifelse(combined_results$module == "black" & as.numeric(combined_results$position) >= 230, "hub", combined_results$hub)
combined_results$hub <- ifelse(combined_results$module == "red" & as.numeric(combined_results$position) >= 807, "hub", combined_results$hub)

png(file = "WGCNA/documentation/pics/spen/filtered_wgcna/hub_genes_new_edited.png", width = 10000, height = 6000, 
    res=600);

(p_1 <- combined_results %>%
  ggplot(aes(x = as.numeric(position), y = egnvctr, color = hub)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("black", "red"))+
    facet_wrap(~module, scales = "free_x") 
)

dev.off()

write.csv( combined_results,"/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/WGCNA/spen/HUB/2025_01_08_spen_hub.csv")

