#!/bin/bash
# Perform BUSCO-analysis to derive phylogenetic tree
# Severin Einspanier

# Load necessary modules
module load gcc12-env
module load singularity

# Set proxy settings
export http_proxy=http://10.0.7.235:3128
export https_proxy=http://10.0.7.235:3128
export ftp_proxy=http://10.0.7.235:3128

# Bind paths for Singularity
export SINGULARITY_BIND="/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/:/mnt"

# Define species and corresponding input/output files
species=("pennellii" "pimpinellifolium" "habrochaites" "lycopersicoides" "chilense")
input_files=(
    "Filtered_SwissProt_pennellii_longest_orfs_pep_clean.fasta"
    "Filtered_SwissProt_pimpinellifolium_longest_orfs_pep_clean.fasta"
    "Filtered_SwissProt_habrochaites_longest_orfs_pep_clean.fasta"
    "Filtered_SwissProt_lycopersicoides_longest_orfs_pep_clean.fasta"
    "Filtered_SwissProt_chilense_longest_orfs_pep_clean.fasta"
)
output_dirs=(
    "data/BUSCO/spen"
    "data/BUSCO/spimp"
    "data/BUSCO/shabro"
    "data/BUSCO/slyco"
    "data/BUSCO/schil"
)
input_files_ed=(
    "Filtered_SwissProt_pennellii_longest_orfs_pep_clean_rm_astrex.fasta"
    "Filtered_SwissProt_pimpinellifolium_longest_orfs_pep_clean_rm_astrex.fasta"
    "Filtered_SwissProt_habrochaites_longest_orfs_pep_clean_rm_astrex.fasta"
    "Filtered_SwissProt_lycopersicoides_longest_orfs_pep_clean_rm_astrex.fasta"
    "Filtered_SwissProt_chilense_longest_orfs_pep_clean_rm_astrex.fasta"
)

# Loop through each species and run BUSCO
for i in "${!species[@]}"; do
    input_file="${input_files[$i]}"
    input_file_ed="${input_files_ed[$i]}"
    output_dir="${output_dirs[$i]}"
    edited_file="${input_file%.fasta}_rm_astrex.fasta"
    
    # Remove asterisks from the input file
    sed -e 's/\*$//' "$input_file" > "$edited_file"
    
    # Run BUSCO using Singularity
    singularity exec --cleanenv BUSCO.sif busco -i "$input_file_ed" -m protein -l solanales_odb10 -o "$output_dir" -f
done