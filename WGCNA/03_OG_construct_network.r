# QDR RNAseq Solanum species
# Construct final network for OGs
# also filter network for Cytoscape
# Define Hub genes
# Severin Einspanier

rm(list=ls())
library(tidyverse)
library(WGCNA)
allowWGCNAThreads(n=30)

options(stringsAsFactors = FALSE)
datExpr <- read.csv("OGs/data/orthogroup_expression_data_rlog_all_species.csv")%>%
  column_to_rownames("X")

datExpr <- as.data.frame(t(datExpr))

# Get Network 

allowWGCNAThreads(30)

result <- blockwiseModules(datExpr, checkMissingData = FALSE, replaceMissingAdjacencies = TRUE, 
                               maxBlockSize = ncol(datExpr), 
                               networkType = "signed hybrid",
                               power =9, 
                               TOMType = "signed", 
                               minModuleSize = 30,
                               corType = "bicor", 
                               reassignThreshold = 0.01, 
                               mergeCutHeight = 0.35, 
                               deepSplit = 0,
                               detectCutHeight = 0.95,
                               numericLabels = TRUE,
                               saveTOMs = F, 
                               saveTOMFileBase = "TOM_sft_9_ds0_cm3",
                               verbose = 3)

colors=labels2colors(result$colors)

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

write.table(gene_module_df,"OG/TOM/TOM_sft_9_ds0_mch35_moduleColors_genid.txt")


# Plot Dendrogram 

png(file=paste0("WGCNA/documentation/pics/OG/", Sys.Date(), "_final_network_pub.png"),
  width = 2000, height = 2000, res=600);

svg(file=paste0("WGCNA/documentation/pics/OG/", Sys.Date(), "_final_network_pub.svg") , 
width = 7, height = 5)
# Plot the dendrogram and module colors
plotDendroAndColors(result$dendrograms[[1]], 
                      colors[result$blockGenes[[1]]], 
                      c(""), 
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = F, guideHang = 0.05, 
                    autoColorHeight = T,
                    colorHeight = .2)

dev.off()

moduleColors=labels2colors(result$colors)
write.table(moduleColors, "OG/TOM_sft_9_ds0_mch35_moduleColors.txt", 
    col.names = FALSE, row.names = TRUE, quote = FALSE, sep = "\t")


load("OG/TOM_sft_9_ds0_cm3-block.1.RData")
moduleColors = read.table("OG/TOM_sft_9_ds0_mch35_moduleColors.txt", header = FALSE)
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

top_1k_edges <- out_network$edgeData %>%
  arrange(desc(weight)) %>%
  head(1000) %>% 
  tail(1) %>% 
  select(weight)

top_10k_edges <- out_network$edgeData %>%
  arrange(desc(weight)) %>%
  head(10000)%>% 
  tail(1) %>% 
  select(weight)

top_100k_edges <- out_network$edgeData %>%
  arrange(desc(weight)) %>%
  head(100000)%>% 
  tail(1) %>% 
  select(weight)

top_1M_edges <- out_network$edgeData %>%
  arrange(desc(weight)) %>%
  head(1000000)%>% 
  tail(1) %>% 
  select(weight)



  # EXPORT 

out_network_threshold <- exportNetworkToCytoscape(
  TOM, 
   edgeFile = paste0("OG/CYTOSCAPE/TOM_1k_edge_", Sys.Date(), ".tsv"),
   nodeFile = paste0("OG/CYTOSCAPE/TOM_1k_node_", Sys.Date(), ".tsv"),
   weighted = TRUE,
   threshold = top_1k_edges$weight,
   nodeNames = names(datExpr),
   #altNodeNames = modGenes,
   nodeAttr = moduleColors,
   includeColNames = TRUE)

out_network_threshold_10k <- exportNetworkToCytoscape(
  TOM, 
   edgeFile = paste0("OG/CYTOSCAPE/TOM_10k_edge_", Sys.Date(), ".tsv"),
   nodeFile = paste0("OG/CYTOSCAPE//TOM_10k_node_", Sys.Date(), ".tsv"),
   weighted = TRUE,
   threshold = top_10k_edges$weight,
   nodeNames = names(datExpr),
   #altNodeNames = modGenes,
   nodeAttr = moduleColors,
   includeColNames = TRUE)

out_network_threshold_100k <- exportNetworkToCytoscape(
  TOM, 
   edgeFile = paste0("OG/CYTOSCAPE//TOM_100k_edge_", Sys.Date(), ".tsv"),
   nodeFile = paste0("OG/CYTOSCAPE//TOM_100k_node_", Sys.Date(), ".tsv"),
   weighted = TRUE,
   threshold = top_100k_edges$weight,
   nodeNames = names(datExpr),
   #altNodeNames = modGenes,
   nodeAttr = moduleColors,
   includeColNames = TRUE)

out_network_threshold_1M <- exportNetworkToCytoscape(
  TOM, 
   edgeFile = paste0("OG/CYTOSCAPE//TOM_1M_edge_", Sys.Date(), ".tsv"),
   nodeFile = paste0("OG/CYTOSCAPE/TOM_1M_node_", Sys.Date(), ".tsv"),
   weighted = TRUE,
   threshold = top_1M_edges$weight,
   nodeNames = names(datExpr),
   #altNodeNames = modGenes,
   nodeAttr = moduleColors,
   includeColNames = TRUE)

exportNetworkToCytoscape(
  TOM, 
   edgeFile = paste0("OG/CYTOSCAPE//TOM_full_edge_", Sys.Date(), ".tsv"),
   nodeFile = paste0("OG/CYTOSCAPE/TOM_full_node_", Sys.Date(), ".tsv"),
   weighted = TRUE,
   nodeNames = names(datExpr),
   threshold =0,
   #altNodeNames = modGenes,
   nodeAttr = moduleColors,
   includeColNames = TRUE)
   

# Get Hub genes 

net_edges <- read.table("OG/CYTOSCAPE//TOM_full_edge_2024-11-28.tsv", header=T)

net_nodes <- read.table("OG/CYTOSCAPE//TOM_full_node_2024-11-28.tsv", header=T, sep = "\t") 


module_df <- net_nodes %>% 
  select(nodeName, nodeAttr.nodesPresent...) %>%
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
    select(egnvctr) %>%
    arrange((egnvctr)) %>%
    rownames_to_column("position")
  brekpont <- brkpnt_fun(sorted$egnvctr,5)

  # Create hub variable and name the data frame
  hub_col <- egnvctr_col %>%
    arrange(egnvctr) %>%
    select(gene, egnvctr) %>%
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
png(file = "WGCNA/documentation/pics/OG/hub_genes_sft9_ds0_mch035.png", width = 10000, height = 6000, 
    res=600);

(p_1 <- combined_results %>%
  ggplot(aes(x = as.numeric(position), y = egnvctr, color = hub)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("black", "red"))+
    facet_wrap(~module, scales = "free_x") 
)

dev.off()


combined_backup <- combined_results

combined_results$hub <- ifelse(combined_results$module == "brown" & as.numeric(combined_results$position) >= 450, "hub", combined_results$hub)
combined_results$hub <- ifelse(combined_results$module == "greenyellow" & as.numeric(combined_results$position) >= 60, "hub", combined_results$hub)
combined_results$hub <- ifelse(combined_results$module == "magenta" & as.numeric(combined_results$position) >= 74, "hub", combined_results$hub)
combined_results$hub <- ifelse(combined_results$module == "pink" & as.numeric(combined_results$position) >= 77, "hub", combined_results$hub)
combined_results$hub <- ifelse(combined_results$module == "purple" & as.numeric(combined_results$position) >= 70, "hub", combined_results$hub)


png(file = "WGCNA/documentation/pics/OG/hub_genes_sft9_ds0_mch035_edited.png", width = 10000, height = 6000, 
    res=600);

(p_1 <- combined_results %>%
  ggplot(aes(x = as.numeric(position), y = egnvctr, color = hub)) +
  geom_point() +
  theme_classic() +out <- read.csv("WGCNA/data/hubs/2024_11_18_slyco_hub.csv") %>% 
  .[grep(temp, ., ignore.case = T), ]


  scale_color_manual(values = c("black", "red"))+
    facet_wrap(~module, scales = "free_x") 
)

dev.off()

# write list of genes
write.csv( combined_results,"OG/2024_11_28_OG_hub.csv")
