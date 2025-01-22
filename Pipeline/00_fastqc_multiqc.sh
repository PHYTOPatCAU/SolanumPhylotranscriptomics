#!/bin/bash

DIR="/gxfs_work/cau/suaph281/RNAseq/RNAseq_work/temp_trim_cuta"

for file in $DIR/cutadapt/*.fastq.gz; do
    fastqc "$file"
done

multiqc "$DIR/cutadapt/."

echo "done"