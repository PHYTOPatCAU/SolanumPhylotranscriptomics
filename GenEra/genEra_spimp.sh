#!/bin/bash
#SBATCH --job-name=GenEra_spimp                   # Job name
#SBATCH --output=GenEra_spimp.out  # Output file
#SBATCH --error=GenEra_spimp.err   # Error file
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --tasks-per-node=2                # Number of tasks (processes)
#SBATCH --cpus-per-task=16                # Number of CPU cores per task
#SBATCH --mem=200G                        # Memory per node (adjust as needed)
#SBATCH --time=24:00:00                   # Maximum runtime (HH:MM:SS)

module load gcc12-env
module load singularity


# mktemp 

# Directory to check/create
DIR="$WORK/RNAseq/RNAseq_work/data/TAI/2024_11_21/spimp/"

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





#cd $WORK/RNAseq/references/spim/spim_GeneExt_CDS.fasta.transdecoder_dir/

#cat longest_orfs.pep | sed 's/ .*//' > longest_orfs_noasterics.fasta

#cd $WORK/RNAseq/RNAseq_work/TAI


cd $WORK/RNAseq/RNAseq_work/data/TAI

export SINGULARITY_BIND="/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data//:/mnt"

DIR="TAI/2024_11_21/spimp"



singularity exec --cleanenv $HOME/programs/GenEra/genEra.sif genEra \
    -q /mnt/PROTEOME/spimp_curated_proteome_OG_pannzer_dedub.fasta \
    -t 4084 -b /mnt/TAI/db2/nr \
    -d /mnt/TAI/db2/taxdump/ \
    -x /mnt/"$DIR" \
    -o /mnt/"$DIR" \
    -n 32


echo -e "\e[32mDONE genEra\e[0m"
echo "S_pimpinellifolium_genEra_pipeline_DONE"

jobinfo