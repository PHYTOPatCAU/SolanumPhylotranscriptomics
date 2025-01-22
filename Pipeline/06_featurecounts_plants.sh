#!/bin/bash
# Perform featureCounts on the mapped reads and curated GTF
# Severin Einspanier

# Load necessary modules and activate conda environment
module load gcc12-env/12.3.0
module load miniconda3/23.5.2
conda activate featurecounts

# Define directories for each species
species_dirs=("schil" "spim" "spen" "shabr" "slyco")

# Loop through each species directory and run featureCounts
for species in "${species_dirs[@]}"; do
    echo "${species^^}"  # Print species name in uppercase

    # Set paths for BAM files, output feature counts, and GTF file
    BAM="data/mapped/clean/$species/"
    FEATURE_out="data/feature_counts/$species/counts_table_GeneExtMM"
    GTF="data/GeneExt/$species/${species}_GeneExt.gtf"

    # Run featureCounts
    featureCounts -s 1 -T 16 -M -t exon -g gene_id -a "$GTF" -o "$FEATURE_out" "${BAM}"*_chr_rename.bam

    # Process the output to create a matrix
    grep -v '#' "$FEATURE_out" | cut -d$'\t' -f1,7- > "${FEATURE_out}_matrix"
done

# IMPORTANT: multimapping allowed increased featureCounts to 70%!