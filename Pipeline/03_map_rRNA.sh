#!/bin/bash

# Map the reads onto the rRNAome
# Severin Einspanier

# define params
THREADS=64

dirs=("schil" "shabr" "slyco" "spen" "spim")


for folder in "${dirs[@]}"; do
    GENOME="$folder/index_complete_uncut_rRNA/"
    Sample_dir="cut_trimmed/$folder/"
    for file in "$Sample_dir"/*.fastq.gz; do
        #echo $file
        BASENAME=$(basename "$file" .fastq.gz)  
        #SAMPLE="$file"
        OUTPUT_PREFIX="rRNA_depleted/$folder/${BASENAME}_rRNA_depleted"
        #echo $SAMPLE
        #gzip -d $file
        SAMPLE="${file}"
        $HOME/programs/STAR-2.7.9a/source/STAR \
            --readFilesCommand zcat \
            --runThreadN $THREADS \
            --genomeDir $GENOME \
            --readFilesIn $SAMPLE \
            --outFilterType BySJout \
            --outFilterMultimapNmax 10 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.2 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outSAMattributes NH HI NM MD \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix $OUTPUT_PREFIX \
            --limitBAMsortRAM 1091046991 \
            --outReadsUnmapped Fastx
    done
done

echo "done MAPPING"
dirs=("schil" "shabr" "slyco" "spen" "spim")

# Extract mapping stats

for folder in "${dirs[@]}"; do
    for file in "$folder"/*Log.final.out; do
        echo "$file"
        ID="${file:5:8}"
        MAPPING=$(sed '10q;d' "$file" | awk -F '|' '{print $2}')
        READS_IN=$(sed '6q;d' "$file" | awk -F '|' '{print $2}')
        READS_mapped=$(sed '9q;d' "$file" | awk -F '|' '{print $2}')
        echo -e "$ID\t$folder\t$READS_IN\t$READS_mapped\t$MAPPING" >> mapping_stats_new.txt
    done
done
echo "finished"

