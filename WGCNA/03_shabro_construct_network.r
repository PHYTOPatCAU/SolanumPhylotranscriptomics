# QDR RNAseq Solanum species
# Construct final network for S. habrochaites
# also filter network for Cytoscape
# Define Hub genes
# Severin Einspanier

rm(list=ls())
library(tidyverse)
library(WGCNA)

setwd("")
options(stringsAsFactors = FALSE)
ids <- read.delim("shabro_curated_proteome_OG_pannzer_dedub_ids.txt", header=F) %>%
    mutate(gene=gsub(">", "", V1)) %>%
    mutate(gene=gsub("GeneExt~", "", gene))%>% 
    mutate(gene=gsub("mRNA:", "", gene))%>%  
    mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
    mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
    mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
    mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))%>%
    dplyr::select(gene)

datExpr <- read.csv("DeSeq/data/norm_counts_all_rlog.csv")%>%
  dplyr::filter(species=="S. habrochaites"& Geneid %in% ids$gene) %>% 
  dplyr::select(!species & !source) %>% 
  pivot_wider(names_from = Geneid, values_from = normalized_count) %>%
  column_to_rownames("sample")

gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

# get network 

allowWGCNAThreads(n=30)

result <- blockwiseModules(datExpr, checkMissingData = FALSE, replaceMissingAdjacencies = TRUE, 
                               maxBlockSize = ncol(datExpr), 
                               networkType = "signed hybrid",
                               power =9, 
                               TOMType = "signed", 
                               minModuleSize = 30,
                               corType = "bicor", 
                               reassignThreshold = 0.01, 
                               mergeCutHeight = 0.2, 
                               deepSplit = 1,
                               detectCutHeight = 0.99,
                               numericLabels = TRUE,
                               saveTOMs = T, 
                               saveTOMFileBase = "shabr/TOM/TOM_final",
                               verbose = 3)

colors=labels2colors(result$colors)

png(file = "WGCNA/documentation/pub/shabro_network_filtered_sft9.png", width = 2000, height = 2000, 
    res=600);
plotDendroAndColors(result$dendrograms[[1]], 
                      colors[result$blockGenes[[1]]], 
                      c(""), 
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = F, guideHang = 0.05, 
                    autoColorHeight = T,
                    colorHeight = .2)
dev.off()

svg(file = "WGCNA/documentation/pub/shabro_network_filtered_sft9.svg", width = 7, height = 5)
plotDendroAndColors(result$dendrograms[[1]], 
                      colors[result$blockGenes[[1]]], 
                      c(""), 
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = F, guideHang = 0.05, 
                    autoColorHeight = T,
                    colorHeight = .2)

dev.off()

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
write.table(gene_module_df,"shabr/TOM/TOMfinal_module_colors_geneids.txt")

stop("stop here")

# Export to Cytoscape:
load("shabr/TOM/TOM_final-block.1.RData")
moduleColors = read.table("shabro/TOM/TOMfinal_module_colors_geneids.txt", sep=" ", row.names=1, header = T)
moduleColors  = as.character(moduleColors$ModuleColor)
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
  edgeFile = paste0("shabro/CYTOSCAPE/TOM_1k_edge_", Sys.Date(), "_filtered.tsv"),
  nodeFile = paste0("shabro/CYTOSCAPE/TOM_1k_node_", Sys.Date(), "_filtered.tsv"),
  weighted = TRUE,
  threshold = top_1k_edges$weight,
  nodeNames = names(datExpr),
  #altNodeNames = modGenes,
  nodeAttr = moduleColors,
  includeColNames = TRUE
)

out_network_threshold_10k <- exportNetworkToCytoscape(
  TOM, 
  edgeFile = paste0("shabro/CYTOSCAPE/TOM_10k_edge_", Sys.Date(), "_filtered.tsv"),
  nodeFile = paste0("shabro/CYTOSCAPE/TOM_10k_node_", Sys.Date(), "_filtered.tsv"),
  weighted = TRUE,
  threshold = top_10k_edges$weight,
  nodeNames = names(datExpr),
  #altNodeNames = modGenes,
  nodeAttr = moduleColors,
  includeColNames = TRUE
)

out_network_threshold_100k <- exportNetworkToCytoscape(
  TOM, 
  edgeFile = paste0("shabro/CYTOSCAPE/TOM_100k_edge_", Sys.Date(), "_filtered.tsv"),
  nodeFile = paste0("shabro/CYTOSCAPE/TOM_100k_node_", Sys.Date(), "_filtered.tsv"),
  weighted = TRUE,
  threshold = top_100k_edges$weight,
  nodeNames = names(datExpr),
  #altNodeNames = modGenes,
  nodeAttr = moduleColors,
  includeColNames = TRUE
)

out_network_threshold_1M <- exportNetworkToCytoscape(
  TOM, 
  edgeFile = paste0("shabro/CYTOSCAPE/TOM_1M_edge_", Sys.Date(), "_filtered.tsv"),
  nodeFile = paste0("shabro/CYTOSCAPE/TOM_1M_node_", Sys.Date(), "_filtered.tsv"),
  weighted = TRUE,
  threshold = top_1M_edges$weight,
  nodeNames = names(datExpr),
  #altNodeNames = modGenes,
  nodeAttr = moduleColors,
  includeColNames = TRUE
)

exportNetworkToCytoscape(
  TOM, 
   edgeFile = paste0("shabro/CYTOSCAPE//TOM_full_edge_", Sys.Date(), "_filtered.tsv"),
   nodeFile = paste0("shabro/CYTOSCAPE/TOM_full_node_", Sys.Date(), "_filtered.tsv"),
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

breakpont <- brkpnt_fun(out_network$edgeData$weight,5,paste0("WGCNA/documentation/pics/shabro/filtered_wgcna/", Sys.Date(), "_shabro_net_edges_filtered.tsv.png"))

# I'll select 0.0495303146022321

out_network_threshold <- exportNetworkToCytoscape(
  TOM, 
   edgeFile = paste0("shabro/CYTOSCAPE/TOM_INFL_filtered_edge_", Sys.Date(), "_filtered.tsv"),
   nodeFile = paste0("shabro/CYTOSCAPE/TOM_INFL_filtered_node_", Sys.Date(), "_filtered.tsv"),
   weighted = TRUE,
   threshold = 0.0495303146022321,
   nodeNames = names(datExpr),
   #altNodeNames = modGenes,
   nodeAttr = moduleColors,
   includeColNames = TRUE)


# Get Hub genes 

net_edges <- read.table("shabro/CYTOSCAPE/TOM_INFL_filtered_edge_2025-01-10_filtered.tsv", header=T)

net_nodes <- read.table("shabro/CYTOSCAPE//TOM_INFL_filtered_node_2025-01-10_filtered.tsv", header=T, sep = "\t") 

module_df <- net_nodes %>% 
  dplyr::select(nodeName, nodeAttr.nodesPresent...) %>%
  rename(gene_id = nodeName, colors = nodeAttr.nodesPresent...) 

head(module_df)
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

combined_results <- unique(module_df$colors) %>%
  map_dfr(process_module)

# Plot the results
png(file = "WGCNA/documentation/pics/shabro/filtered_wgcna/hub_genes_new.png", width = 10000, height = 6000, 
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

combined_results$hub <- ifelse(combined_results$module == "tan" & as.numeric(combined_results$position) >= 61, "hub", combined_results$hub)
combined_results$hub <- ifelse(combined_results$module == "salmon" & as.numeric(combined_results$position) >= 45, "hub", combined_results$hub)
combined_results$hub <- ifelse(combined_results$module == "purple" & as.numeric(combined_results$position) >= 108, "hub", combined_results$hub)
combined_results$hub <- ifelse(combined_results$module == "magenta" & as.numeric(combined_results$position) >= 150, "hub", combined_results$hub)

png(file = "WGCNA/documentation/pics/shabro/filtered_wgcna/hub_genes_new_edited.png", width = 10000, height = 6000, 
    res=600);

(p_1 <- combined_results %>%
  ggplot(aes(x = as.numeric(position), y = egnvctr, color = hub)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("black", "red"))+
    facet_wrap(~module, scales = "free_x") 
)

dev.off()

write.csv(combined_results, "shabro/HUB/2025_01_10_shabro_hub.csv")