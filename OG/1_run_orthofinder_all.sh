#!/bin/bash
#SBATCH --job-name=OF_WholeP                   # Job name
#SBATCH --output=OF_WholeP.out  # Output file
#SBATCH --error=OF_WholeP.err   # Error file
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --tasks-per-node=2                # Number of tasks (processes)
#SBATCH --cpus-per-task=16                # Number of CPU cores per task
#SBATCH --mem=128G                         # Memory per node (adjust as needed)
#SBATCH --time=3:00:00                   # Maximum runtime (HH:MM:SS)

# QDR RNAseq Solanum species
# Orthofinder Analysis on 5 Solanum species and pentapetlaeae
# Whole proteome analysis
# Severin Einspanier

module load gcc12-env
module load miniconda3

conda activate orthofinder_env

echo -e "\e[32mStarting Orthofinder Analysis\e[0m"
COREdir="whole_proteome/SwissItagPannzer_filtered/"

orthofinder -f "$COREdir" -t 32

echo "DONE"
