#!/bin/bash
#Create rRNA-ome reference file for all species
#Severin Einspanier
#2024_06

# Load necessary modules and activate conda environment
module load gcc12-env/12.3.0
module load miniconda3/23.5.2
conda activate samtools

# List of folders to process
folder=(
    "schil"
    "shabr"
    "slyco"
    "spen"
    "spim"
    )

# Loop through each folder
for folder in "${folder[@]}"; do
  # Move into the folder (assuming subfolders are named with trailing slash)
  cd "$folder"

  # Remove existing rRNA sequence files if they exist
  if [ -e "rRNA_sequences.bed" ]; then
    rm "rRNA_sequences.bed"
    rm "rRNA_sequences_ordered.bed"
    rm "rRNA_sequences_*"
  fi
  
  # Create a new rRNA sequences file
  touch rRNA_sequences.bed

  # Loop through each .txt file in the current folder
  for file in *tab.txt; do
    # Read the first line of the file
    first_line=$(head -n 1 "$file")

    # Extract title (3rd column), start (10th column), and end (11th column)
    title=$(echo "$first_line" | awk -F'\t' '{print $2 }')
    start=$(echo "$first_line" | awk -F'\t' '{print $10 }')
    end=$(echo "$first_line" | awk -F'\t' '{print $11 }')

    # Write data to the rRNA sequences file with tab separation
    echo -e "$title\t$start\t$end" >> rRNA_sequences.bed
  done

  # Adjust coordinates and create ordered rRNA sequences file
  awk -F'\t' '{if ($3 == 61731525) {print $1"\t"($2-100)"\t"($3)} else {{if ($2 > $3) {print $1"\t"($3-100)"\t"($2+100)} else {print $1"\t"($2-100)"\t"($3+100)}}}}' rRNA_sequences.bed > rRNA_sequences_ordered.bed
  
  # Move back to the parent directory
  cd ..
done

# Extract sequences from fasta using coordinates
bedtools getfasta -fi GCA_001406875.2_SPENNV200_genomic_top12.fasta  -bed rRNA_sequences_ordered.bed -fo rRNA_sequences_spen.fa
bedtools getfasta -fi S_chilense_Hirise_top12.fasta  -bed rRNA_sequences_ordered -fo rRNA_sequences_schil.fa
bedtools getfasta -fi GWHBJTH00000000.genome_top12.fasta  -bed rRNA_sequences_ordered -fo rRNA_sequences_habro.fa
bedtools getfasta -fi GCA_022817965.1_SlydLA2951_v2.0_genomic_top12.fasta -bed rRNA_sequences_ordered -fo rRNA_sequences_lyco.fa
bedtools getfasta -fi LA2093_genome_top12.fasta -bed rRNA_sequences_ordered -fo rRNA_sequences_pimp.fa