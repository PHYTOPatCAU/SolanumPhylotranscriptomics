# QDR RNAseq Solanum species
# Extract geneIDs and merge into one consensus
# Severin Einspanier 

pacman::p_load(tidyverse)
setwd("")
# List all files matching the pattern
files <- list.files("Ogs/data", pattern = "^genid2goid_", full.names = TRUE)

# Initialize an empty list to store processed data frames
processed_list <- list()

# Loop over each file
for (f in files) {
  # Read the data; adjust read.table arguments to your file format
  df <- read.csv(f, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  
  # Process the data
  processed_df <- df %>%
    filter(GeneID != "") %>%
    mutate(
      gene = gsub("GeneExt~", "", GeneID),
      gene = gsub("mRNA_", "", gene),
      gene = gsub("t\\.peak", "g.peak", gene),
      gene = gsub("t\\.minus", "g.minus", gene),
      gene = gsub("t\\.plus", "g.plus", gene),
      gene = gsub("\\.[1-9].*|\\.p[1-9].*", "", gene)
    ) %>%
    select(OG, gene, species) %>%
    unique()

  # Store the processed df in the list
  processed_list[[f]] <- processed_df
}

# Combine all processed data frames into one
merged_df <- bind_rows(processed_list)

# If you need to ensure uniqueness across all files, you can do so here:
merged_df <- unique(merged_df)

# Write the merged data to a file
write.table(merged_df, file = paste0("data/",
  Sys.Date(), "_genid2goid_all_species.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  