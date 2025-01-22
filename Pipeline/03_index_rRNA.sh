#!/bin/bash
# Generate Genome indices of rRNA references
# Severin Einspanier

THREADS=16
module load gcc12-env/12.3.0
cd  

dirs=("schil" "shabr" "slyco" "spen" "spim")

for folder in "${dirs[@]}"; do
    if [[ -d "$folder" ]]; then
    # Change directory (use pushd/popd for safer stack management)
        pushd "$folder" >/dev/null 2>&1  # Suppress potential output from pushd

        mkdir -p index_complete_rRNA/

        for file in complete_rRNA_uncut.fasta; do
            GENOME_fasta="$file"
            GENOME_dir="index_complete_uncut_rRNA/"
            /gxfs_home/cau/suaph281/programs/STAR-2.7.9a/source/STAR \
                --runThreadN "$THREADS" \
                --runMode genomeGenerate \
                --genomeDir "$GENOME_dir" \
                --genomeFastaFiles "$GENOME_fasta" \
                --genomeSAindexNbases 5
        done
        
        #for file in rRNA_sequences*.fa; do
        #    GENOME_fasta="$file"
        #    GENOME_dir="index/"
        #    /gxfs_home/cau/suaph281/programs/STAR-2.7.9a/source/STAR \
        #        --runThreadN "$THREADS" \
        #        --runMode genomeGenerate \
        #        --genomeDir "$GENOME_dir" \
        #        --genomeFastaFiles "$GENOME_fasta" \
        #        --genomeSAindexNbases 5
        #done
         # Pop back from the directory stack (assuming pushd was used)
        popd >/dev/null 2>&1  # Suppress potential output from popd
    fi
done




