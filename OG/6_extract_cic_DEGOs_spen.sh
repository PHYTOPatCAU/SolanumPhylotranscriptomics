# QDR RNAseq Solanum species
# get cic DEGs for functional enrichment using spen data
# Severin Einspanier

#!/bin/bash

# Load necessary modules
module load gcc12-env
module load miniconda3

# Activate the seqkit Conda environment
conda activate seqkit

# Path to the CSV file
CSV_FILE="OG_cDEOGs_IDs_spen.csv"

# Path to the species-proteome FASTA file
PROTEOME_FILE="spen_curated_proteome_OG_pannzer_dedub.fasta"

# Path to the output FASTA file
OUTPUT_FILE="all_extracted_sequences.fasta"

# Initialize the output file
> "$OUTPUT_FILE"

# Read each sequence_id from the CSV and extract corresponding sequences
while IFS= read -r sequence_id
do
    echo "$sequence_id"
    # Check if the line is not empty
    if [[ -n "$sequence_id" ]]; then
        # Extract the sequence using seqkit with partial match and append to the output file
        seqkit grep -rp "$sequence_id" "$PROTEOME_FILE" >> "$OUTPUT_FILE"
    fi
done < "$CSV_FILE"