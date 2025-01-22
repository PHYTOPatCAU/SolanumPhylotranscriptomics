#!/bin/bash

module load gcc12-env/12.3.0
module load blast-plus

# Database folders (modify these paths to your actual data)
databases=(
  "S_chilense_Hirise_top12.fasta"
  "GWHBJTH00000000.genome_top12.fasta"
  "GCA_022817965.1_SlydLA2951_v2.0_genomic_top12.fasta"
  "GCA_001406875.2_SPENNV200_genomic_top12.fasta"
  "LA2093_genome_top12.fasta"
)

# Query sequences (modify these paths to your actual data)
queries=(
  "5_8S_AB373816.fa"
  "25S_OK073662.fa"
  "18S_OK073663.fa"
)

# BLAST parameters (adjust these as needed)
evalue=1e-6
num_threads=16

# Loop through each database and query combination
for database in "${databases[@]}"; do
  for query in "${queries[@]}"; do
    # Extract filename from database path (assuming filename is at the end)
    basename_db=$(basename "$database")
    dirname_db=$(dirname "$database")

    # Construct output filename with desired format
    output_filename_read="$dirname_db"/blast_output/"$(basename "$query")_blast_${basename_db%.fasta}_read.txt"

    # Run BLAST command for each combination
    blastn -db "$database" -query "$query" -evalue "$evalue" -num_threads "$num_threads" -out "$output_filename_read"
    output_filename="$dirname_db"/blast_output/"$(basename "$query")_blast_${basename_db%.fasta}_tab.txt"

    blastn -db "$database" -query "$query" -evalue "$evalue" -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs 5 -num_threads "$num_threads" -out "$output_filename"

    echo "Running BLAST:  $database vs $query output: $output_filename"
  done
done

echo "All BLAST searches completed."