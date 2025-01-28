#!/bin/bash
#SBATCH --job-name=run_proteome                  # Job name
#SBATCH --output=run_proteome.out                # Output file
#SBATCH --error=run_proteome.err                 # Error file
#SBATCH --nodes=1                                # Number of nodes
#SBATCH --tasks-per-node=2                       # Number of tasks (processes)
#SBATCH --cpus-per-task=16                       # Number of CPU cores per task
#SBATCH --mem=200G                               # Memory per node (adjust as needed)
#SBATCH --time=2:00:00                           # Maximum runtime (HH:MM:SS)

# QDR RNAseq Solanum species
# Perform Proteome-Filtering (to exclude spuriouse sequences) 
# multiple-tier appraoch: compare with uniprot, ITAG4-orthology and last Pannzer2 hits
# Severin Einspanier

# Load necessary modules
module load gcc12-env
module load miniconda3

# Activate the seqkit conda environment
conda activate seqkit

# Function to filter sequences that did not survive the Swissprot filter
filter_sequences() {
    local file1="$1"           # Raw transdecoder file
    local file2="$2"           # Swiss-prot-filtered file
    local output_file="$3"     # Output file for sequences that did not survive the Swissprot filter
    local log_file="$4"        # Log file

    # Echo the start of the function
    echo "Starting sequence filtering for $file1 and $file2..."

    # Make temporary ID files
    echo "Creating temporary ID files..."
    cat "$file1" | grep ">" | sed 's/>//' | cut -d " " -f 1 > "temp_ids_transdecoder.txt"
    cat "$file2" | grep ">" | sed 's/>//' | cut -d " " -f 1 > "temp_ids_swiss.txt"

    # Get the IDs that are not in the Swissprot file
    echo "Filtering sequences not in the Swissprot file..."
    grep -vf "temp_ids_swiss.txt" "temp_ids_transdecoder.txt" > "$output_file"

    # Log results
    echo "Logging results..."
    {
        echo "Original sequences count:"
        grep ">" "$file1" | wc -l
        echo "Swissprot Filtered (good) sequences count:"
        grep ">" "$file2" | wc -l
        echo "Filtered (bad) sequences count:"
        cat "$output_file" | wc -l
    } > "$log_file"

    # Clean up temporary files
    echo "Cleaning up temporary files..."
    rm "temp_ids_transdecoder.txt" "temp_ids_swiss.txt"

    # Echo the completion of the function
    echo "Sequence filtering completed. Results saved to $output_file and log saved to $log_file."
}

# Base directory
base_dir="PROTEOME"

# Species list
species_list=("spimp" "schil" "shabro" "slyco" "spen")

# Process each species
for species in "${species_list[@]}"; do
    echo "Processing species: $species"

    transdec_raw=($base_dir/transdecoder_longest/$species/*.fasta)
    swissprot_filtered=($base_dir/swiss_prot_filtered/$species/*_pep_clean.fasta)
    output_file="$base_dir/swiss_prot_filtered/$species/bad_swissprot_sequences_${species}.txt"
    log_file="$base_dir/swiss_prot_filtered/$species/${species}_log.txt"

    badsamples_sequences="$base_dir/orthofinder/${species}.fasta"

    # Filter sequences
    filter_sequences "$transdec_raw" "$swissprot_filtered" "$output_file" "$log_file"

    # Extract the sequences
    echo "Extracting sequences that survived the filter..."
    seqkit grep -f "$output_file" "$transdec_raw" -o "$badsamples_sequences"

    echo "Processing for species $species completed."
done

# Copy ITAG4.1 proteins file
cp ITAG4.1_proteins.fasta PROTEOME/orthofinder/ITAG4.1_proteins.fasta

# Activate orthofinder environment and run orthofinder
conda activate orthofinder_env
cd PROTEOME/orthofinder/
orthofinder -f PROTEOME/orthofinder/proteome/ -t 32

# Isolate orthologues
cd OrthoFinder/Results_Nov21/Orthologues/Orthologues_ITAG4.1_proteins

# Convert GeneExt~mRNA:SHch09g025890.1.p1 to GeneExt~mRNA_SHch09g025890.1.p1
files=("ITAG4.1_proteins__v__spen.tsv" "ITAG4.1_proteins__v__spimp.tsv" "ITAG4.1_proteins__v__schil.tsv" "ITAG4.1_proteins__v__slyco.tsv")

# Get the seuqence IDs of the ITAG orthologues

for file in "${files[@]}"; do
    # Make ID file of orthologues
    extract_ids="${file%.tsv}_ID.txt"
    extract_ids_newline="${extract_ids%.txt}_newline.txt"
    extract_ids_newline_chrtrs="${extract_ids%.txt}_newline_chrtrs.txt"

    cat "$file" | cut -f 3 > "$extract_ids"
    # Change ',' with new line
    wait
    sed 's/, /\n/g' "$extract_ids" > "$extract_ids_newline"
    # sed 's/, /\n/g' "$extract_ids" | sed 's/_/:/g' > "$extract_ids_newline_chrtrs"
    touch "$extract_ids_newline_chrtrs"
    #rm "$extract_ids"
    mv "$extract_ids_newline" PROTEOME/orthofinder/
    mv "$extract_ids_newline_chrtrs" /PROTEOME/orthofinder/
done

# Process Shabro
# needs to be done as shabro had issues with ':' in fasta header. those were replaced with '_' by orthofinder
cat ITAG4.1_proteins__v__shabro.tsv | cut -f 3 > ITAG4.1_proteins__v__shabro_ID.txt
sed 's/, /\n/g' ITAG4.1_proteins__v__shabro_ID.txt > ITAG4.1_proteins__v__shabro_ID_newline.txt
sed 's/, /\n/g' ITAG4.1_proteins__v__shabro_ID.txt | sed 's/_/:/g' > ITAG4.1_proteins__v__shabro_ID_newline_chrtrs.txt
rm ITAG4.1_proteins__v__shabro_ID.txt
mv "ITAG4.1_proteins__v__shabro_ID_newline.txt" PROTEOME/orthofinder/shabro/
mv "ITAG4.1_proteins__v__shabro_ID_newline_chrtrs.txt" PROTEOME/orthofinder/shabro/

# Move files to respective directories
cd PROTEOME/orthofinder/
mv *_schil_ID_newline*.txt schil/
mv *_slyco_ID_newline*.txt slyco/
# mv *_shabro_ID_newline*.txt shabro/
mv *_spen_ID_newline*.txt spen/
mv *_spimp_ID_newline*.txt spimp/

# Quick seqkit for all species
cd PROTEOME/

# Add PANNZER step
conda activate seqkit

for species in "${species_list[@]}"; do
    ITAG_good_sequence_id=($base_dir/orthofinder/${species}/*_ID_newline.txt)
    orthofinder_prot=($base_dir/orthofinder/${species}.fasta)
    ITAG_bad_sequences=($base_dir/orthofinder/${species}/${species}_bad_sequences.fasta)
    ITAG_bad_sequences_id=($base_dir/orthofinder/${species}/${species}_bad_sequences_ID.txt)
    go_term_file=($base_dir/pannzer2_filtered/GO/${species}/*filtered_IDs.txt)
    pannzer_recovered=($base_dir/pannzer2_filtered/${species}/${species}_pannzer_recovered.txt)

    # Filter sequences that dont pass ITAG OG filter
    seqkit grep -vf "$ITAG_good_sequence_id" "$orthofinder_prot" -o "$ITAG_bad_sequences"
    #rm "$ITAG_good_sequence_id"
    # Get IDs of bad sequences
    cat "$ITAG_bad_sequences" | grep ">" > "$ITAG_bad_sequences_id"

    # Filter GO terms
    grep -f "$go_term_file" "$ITAG_bad_sequences_id" | sed 's/>//' > "$pannzer_recovered"
    cat "$pannzer_recovered" | wc -l

    rm "$ITAG_bad_sequences_id"

    transdec_raw=($base_dir/transdecoder_longest/$species/*.fasta)
    swissprot_filtered=($base_dir/swiss_prot_filtered/$species/*_pep_clean.fasta)
    #orthologues=($base_dir/orthofinder/${species}/*_ID_newline.txt)
    orthologues_chrtrs=($base_dir/orthofinder/${species}/*_ID_newline_chrtrs.txt)

    echo "Processing species: $species"
    seqkit grep -f "$ITAG_good_sequence_id" "$transdec_raw" -o "orthofinder/$species/${species}_orthologues.fasta"
    seqkit grep -f "$orthologues_chrtrs" "$transdec_raw" -o "orthofinder/$species/${species}_orthologues_chrts.fasta"
    seqkit grep -nf "$pannzer_recovered" "$transdec_raw" -o "orthofinder/$species/${species}_pannzer_recovered.fasta"

    # Generate new curated proteome
    echo "Generating new curated proteome..."
    cat "$swissprot_filtered" "orthofinder/$species/${species}_orthologues.fasta" "orthofinder/$species/${species}_orthologues_chrts.fasta" "orthofinder/$species/${species}_pannzer_recovered.fasta" > "$base_dir/${species}_curated_proteome_OG_pannzer.fasta"
    cat "$swissprot_filtered" "orthofinder/$species/${species}_orthologues.fasta" "orthofinder/$species/${species}_orthologues_chrts.fasta" > "$base_dir/${species}_curated_proteome_OG.fasta"

    # Remove duplicates
    seqkit rmdup -s < "$base_dir/${species}_curated_proteome_OG.fasta" > "$base_dir/${species}_curated_proteome_OG_dedub.fasta"
    seqkit rmdup -s < "$base_dir/${species}_curated_proteome_OG_pannzer.fasta" > "$base_dir/${species}_curated_proteome_OG_pannzer_dedub.fasta"
    
done

echo -e "\n"
echo "All species processed."

echo -e "\e[32mStarting Length analysis\e[0m"

{
# Analysis for S. chilense
echo -e "\e[34mS.chilense\e[0m"
echo "Unfiltered:"
cat transdecoder_longest/schil/*.fasta | seqkit stats
echo "SwissProteome:"
cat swiss_prot_filtered/schil/*_pep_clean.fasta | seqkit stats
echo -e "\033[33mITAG filter:\033[0m"
cat schil_curated_proteome_OG_dedub.fasta | seqkit stats
echo -e "\033[33mITAG&PANNZER filter:\033[0m"
cat schil_curated_proteome_OG_pannzer_dedub.fasta | seqkit stats

# Analysis for S. lycopersicoides
echo -e "\e[34mS.lycopersicum\e[0m"
echo "Unfiltered:"
cat transdecoder_longest/slyco/*.fasta | seqkit stats
echo "SwissProteome:"
cat swiss_prot_filtered/slyco/*_pep_clean.fasta | seqkit stats
echo -e "\033[33mITAG filter:\033[0m"
cat slyco_curated_proteome_OG_dedub.fasta | seqkit stats
echo -e "\033[33mITAG&PANNZER filter:\033[0m"
cat slyco_curated_proteome_OG_pannzer_dedub.fasta | seqkit stats

# Analysis for S. habrochaites
echo -e "\e[34mS.habrochaites\e[0m"
echo "Unfiltered:"
cat transdecoder_longest/shabro/*.fasta | seqkit stats
echo "SwissProteome:"
cat swiss_prot_filtered/shabro/*_pep_clean.fasta | seqkit stats
echo -e "\033[33mITAG filter:\033[0m"
cat shabro_curated_proteome_OG_dedub.fasta | seqkit stats
echo -e "\033[33mITAG&PANNZER filter:\033[0m"
cat shabro_curated_proteome_OG_pannzer_dedub.fasta | seqkit stats

# Analysis for S. pennellii
echo -e "\e[34mS.pennellii\e[0m"
echo "Unfiltered:"
cat transdecoder_longest/spen/*.fasta | seqkit stats
echo "SwissProteome:"
cat swiss_prot_filtered/spen/*_pep_clean.fasta | seqkit stats
echo -e "\033[33mITAG filter:\033[0m"
cat spen_curated_proteome_OG_dedub.fasta | seqkit stats
echo -e "\033[33mITAG&PANNZER filter:\033[0m"
cat spen_curated_proteome_OG_pannzer_dedub.fasta | seqkit stats

# Analysis for S. pimpinellifolium
echo -e "\e[34mS.pimpinellifolium\e[0m"
echo "Unfiltered:"
cat transdecoder_longest/spimp/*.fasta | seqkit stats
echo "SwissProteome:"
cat swiss_prot_filtered/spimp/*_pep_clean.fasta | seqkit stats
echo -e "\033[33mITAG filter:\033[0m"
cat spimp_curated_proteome_OG_dedub.fasta | seqkit stats
echo -e "\033[33mITAG&PANNZER filter:\033[0m"
cat spimp_curated_proteome_OG_pannzer_dedub.fasta | seqkit stats

echo "All species processed."
} > LOG_TAIR_OG_FILTERING.txt