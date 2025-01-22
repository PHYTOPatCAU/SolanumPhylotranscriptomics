#!/bin/bash

module load gcc12-env/12.3.0
module load miniconda3/23.5.2

# Path to input directory containing paired-end fastq files
input_dir=""

output_dir=""


conda activate cutadapt 
for file in $input_dir/*.fastq.gz; do
    # Extract the sample name from the file name
    sample=$(basename "$file" .fastq.gz)
    cutadapt --quiet -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 $input_dir/${sample}.fastq.gz | \
    cutadapt --quiet -m 20 -O 3 --nextseq-trim=10 -a "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | \
    cutadapt --quiet -m 20 -O 3 -a "r1polyA=A{18}" - | \
    cutadapt --quiet -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o $output_dir/${sample}_trimmed.fastq.gz -
done

conda activate fastqc

# do fastqc
# do fastqc for each FASTQ file in the output directory
for file in $output_dir/*.fastq.gz; do
    fastqc "$file"
done

multiqc "$output_dir/."
echo "done"
