#!/bin/bash
#SBATCH --job-name=GenEra_schil                  # Job name
#SBATCH --output=GenEra/GenEra_schil.out  # Output file
#SBATCH --error=GenEra/GenEra_schil.err   # Error file
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --tasks-per-node=2                # Number of tasks (processes)
#SBATCH --cpus-per-task=16                # Number of CPU cores per task
#SBATCH --mem=200G                        # Memory per node (adjust as needed)
#SBATCH --time=24:00:00                   # Maximum runtime (HH:MM:SS)

module load gcc12-env
module load singularity

# mktemp 

# Directory to check/create
DIR="schil/"

# Check if directory exists
if [ -d "$DIR" ]; then
    echo "Directory '$DIR' already exists."
else
    echo "Directory '$DIR' does not exist. Creating now..."
    mkdir -p "$DIR"
    if [ $? -eq 0 ]; then
        echo "Directory '$DIR' created successfully."
    else
        echo "Failed to create directory '$DIR'."
        exit 1
    fi
fi


export SINGULARITY_BIND="/data//:/mnt"

DIR="TAI/2024_11_21/schil"

singularity exec --cleanenv $HOME/programs/GenEra/genEra.sif genEra \
    -q /mnt/PROTEOME/schil_curated_proteome_OG_pannzer_dedub.fasta  \
    -t 4083  -b /mnt/TAI/db2/nr \
    -d /mnt/TAI/db2/taxdump/ \
    -x /mnt/"$DIR" \
    -o /mnt/"$DIR" \
    -n 32


echo -e "\e[32mDONE genEra\e[0m"
echo "S_chilense_genEra_pipeline_DONE"

jobinfo
