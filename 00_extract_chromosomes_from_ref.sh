#!/bin/bash
# Script to extract only the chromosomes from the reference-genomes
# Severin Einspanier
# Function to extract sequences based on sequence ID
extract_sequence() {
    local id="$1"
    awk -v id="$id" -v found=0 '/^>/{if($0~id) found=1; else found=0} found' "$input_file"
}

# Process S. habrochaites chromosomes
# Define the input and output file names
input_file="GWHBJTH00000000.genome.fasta"
output_file="GWHBJTH00000000.genome_top12.fasta"

# Define the sequence IDs as a single string
sequence_ids="GWHBJTH00000001
GWHBJTH00000002
GWHBJTH00000003
GWHBJTH00000004
GWHBJTH00000005
GWHBJTH00000006
GWHBJTH00000007
GWHBJTH00000008
GWHBJTH00000009
GWHBJTH00000010
GWHBJTH00000011
GWHBJTH00000012"

# Convert the sequence IDs string into an array
IFS=$'\n' read -r -d '' -a sequence_ids_array <<< "$sequence_ids"

# Loop through each sequence ID and extract the sequences
for id in "${sequence_ids_array[@]}"; do
    extract_sequence "$id" >> "$output_file"
done

# Process S. lycopersicoides chromosomes
# Define the input and output file names
input_file="GCA_022817965.1_SlydLA2951_v2.0_genomic.fna"
output_file="GCA_022817965.1_SlydLA2951_v2.0_genomic_top12.fasta"

# Define the sequence IDs as a single string
sequence_ids="CM040618.1
CM040619.1
CM040620.1
CM040621.1
CM040622.1
CM040623.1
CM040624.1
CM040625.1
CM040626.1
CM040627.1
CM040628.1
CM040629.1"

# Convert the sequence IDs string into an array
IFS=$'\n' read -r -d '' -a sequence_ids_array <<< "$sequence_ids"

# Loop through each sequence ID and extract the sequences
for id in "${sequence_ids_array[@]}"; do
    extract_sequence "$id" >> "$output_file"
done

# Process S. pennellii chromosomes
# Define the input and output file names
input_file="Spenn.fasta"
output_file="Spenn_top12.fasta"

# Define the sequence IDs as a single string
sequence_ids="Spenn-ch01
Spenn-ch02
Spenn-ch03
Spenn-ch04
Spenn-ch05
Spenn-ch06
Spenn-ch07
Spenn-ch08
Spenn-ch09
Spenn-ch10
Spenn-ch11
Spenn-ch12"

# Convert the sequence IDs string into an array
IFS=$'\n' read -r -d '' -a sequence_ids_array <<< "$sequence_ids"

# Loop through each sequence ID and extract the sequences
for id in "${sequence_ids_array[@]}"; do
    extract_sequence "$id" >> "$output_file"
done