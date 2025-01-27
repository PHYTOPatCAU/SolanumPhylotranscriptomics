#!/bin/bash
#SBATCH --job-name=GenEra_spen                   # Job name
#SBATCH --output=/gxfs_home/cau/suaph281/debug/LDT_rnaseq/GenEra/GenEra_spen.out  # Output file
#SBATCH --error=/gxfs_home/cau/suaph281/debug/LDT_rnaseq/GenEra/GenEra_spen.err   # Error file
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --tasks-per-node=2                # Number of tasks (processes)
#SBATCH --cpus-per-task=16                # Number of CPU cores per task
#SBATCH --mem=200G                        # Memory per node (adjust as needed)
#SBATCH --time=14:00:00                   # Maximum runtime (HH:MM:SS)

module load gcc12-env
module load singularity

# mktemp 

# Directory to check/create
DIR="$WORK/RNAseq/RNAseq_work/data/TAI/2024_11_21/spen/"

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

# get raw proteome


#singularity exec --B $WORK -e transdecoder.sif  TransDecoder.LongOrfs -t  $WORK/RNAseq/references/spen/spen_GeneExt_CDS.fasta --output_dir  $WORK/RNAseq/references/spen/
# complete only and 200 bp
#cd $WORK/RNAseq/references/spen/spen_GeneExt_CDS.fasta.transdecoder_dir/

#cat longest_orfs.pep | sed 's/ .*//' > longest_orfs_noasterics.fasta

cd $WORK/RNAseq/RNAseq_work/data/TAI

export SINGULARITY_BIND="/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/data//:/mnt"

DIR="TAI/2024_11_21/spen/"

singularity run --cleanenv $HOME/programs/GenEra/genEra.sif genEra \
    -q /mnt/PROTEOME/spen_curated_proteome_OG_pannzer_dedub.fasta \
    -t 28526 -b /mnt/TAI/db/nr \
    -d /mnt/TAI/db/taxdump/ \
    -x /mnt/"$DIR" \
    -o /mnt/"$DIR" \
    -n 32


echo -e "\e[32mDONE genEra\e[0m"
echo "S_pennellii_genEra_pipeline_DONE"

jobinfo
