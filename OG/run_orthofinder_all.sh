#!/bin/bash
#SBATCH --job-name=OF_WholeP                   # Job name
#SBATCH --output=/gxfs_home/cau/suaph281/debug/LDT_rnaseq/OF_WholeP.out  # Output file
#SBATCH --error=/gxfs_home/cau/suaph281/debug/LDT_rnaseq/OF_WholeP.err   # Error file
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --tasks-per-node=2                # Number of tasks (processes)
#SBATCH --cpus-per-task=16                # Number of CPU cores per task
#SBATCH --mem=128G                         # Memory per node (adjust as needed)
#SBATCH --time=3:00:00                   # Maximum runtime (HH:MM:SS)


module load gcc12-env
module load miniconda3

conda activate orthofinder_env

echo -e "\e[32mStarting Orthofinder Analysis\e[0m"
COREdir="/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data/orthofinder/whole_proteome/SwissItagPannzer_filtered/"

orthofinder -f "$COREdir" -t 32

echo "DONE"
