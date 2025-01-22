#!/bin/bash
# Perform Mapping of rRNA-depleted reads on the reference genomes
# Severin Einspanier
# 2024_06

# Set the number of threads to use for STAR
THREADS=16*2

# Load necessary module
module load gcc12-env/12.3.0

# List of directories to process
dirs=("schil" "shabr" "slyco" "spen" "spim")

# Loop through each directory
for folder in "${dirs[@]}"; do
    # Set the genome directory path
    GENOME="$folder/index/"
    
    # Set the sample directory path
    Sample_dir="rRNA_depleted/$folder/"
    
    # Loop through each file in the sample directory
    for file in "$Sample_dir"*trimmed_rRNA_depletedUnmapped.out.mate1; do
        # Check if the file does not have a .bam extension
        if [[ "$file" != *.bam ]]; then
            # Rename the file to add the .bam extension
            mv "$file" "$file.bam"
            file="$file.bam"
        fi
        
        # Get the base name of the file without the .bam extension
        BASENAME=$(basename "$file" .bam)
        echo "$BASENAME"
        
        # Set the output prefix for STAR
        OUTPUT_PREFIX="mapped/$folder/${BASENAME}_sol_sclero_mapped"
        
        # Run STAR with the specified parameters
        $HOME/programs/STAR-2.7.9a/source/STAR \
            --runThreadN $THREADS \
            --genomeDir $GENOME \
            --readFilesIn $file \
            --outFilterType BySJout \
            --outFilterMultimapNmax 200 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverLmax 0.6 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --limitOutSJcollapsed 5000000 \
            --outSAMattributes NH HI NM MD \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix $OUTPUT_PREFIX \
            --limitBAMsortRAM 1091046991 \
            --outReadsUnmapped Fastx
    done
done

# Generate mapping statistics
dirs=("schil" "shabr" "slyco" "spen" "spim")
for folder in "${dirs[@]}"; do
    for file in "$folder"/*Log.final.out; do
        echo "$file"
        
        # Extract the ID from the file name
        ID="${file:5:8}"
        
        # Extract mapping statistics from the log file
        MAPPING=$(sed '10q;d' "$file" | awk -F '|' '{print $2}')
        MULTIMAPPING=$(sed '25q;d' "$file" | awk -F '|' '{print $2}')
        MAPPED=$(sed '9q;d' "$file" | awk -F '|' '{print $2}')
        
        # Write the statistics to the output file
        echo -e "$ID\t$folder\t$MAPPED\t$MAPPING\t$MULTIMAPPING" >> mapping_stats.txt
    done
done
echo "finished"