#!/bin/bash
# Pipeline for Genome-Annotation-Curation
# Using Reads and GeneExt to predict UTRs
# Severin Einspanier 


# Approach: edit GFF/GTF beforehand (remove scaffolds etc) and run AGAT gff2gtf on them
# Perform THIS script to merge bams and remove scaffold/sclero reads from all bams
# Then run GeneExt and FeatureCounts.

module load gcc12-env
module load miniconda3
module load singularity

# 1. Remove scaffolds from GFF and convert to GTF
# By using AGAT, possible issues from the GFFs are resolved as it checks for UTRs, exons and also level nomenclature.

cd references/

## 1.1 S. habrochaites 
awk '$1 ~ /^GWHBJTH0000000[1-9]$/ || $1 ~ /^GWHBJTH0000001[0-2]$/' shabr/GWHBJTH00000000.gff > shabr/GWHBJTH00000000_top12.gff
singularity run -B $WORK/RNAseq/ agat_1.0.0.sif agat_convert_sp_gff2gtf.pl \
    --gff shabr/GWHBJTH00000000_top12.gff  \
    -o shabr/GWHBJTH00000000_annot_top12.gtf
wait
echo "DONE SINGULARITY HABRO" 

## 1.2 S. lycopersicoides
cat slyco/SlydLA2951_v2.0_gene_models_all.gff3 | grep -vE "Contig" > slyco/SlydLA2951_v2.0_gene_models_all_top12.gff3 
singularity run -B $WORK/RNAseq/ agat_1.0.0.sif agat_convert_sp_gff2gtf.pl \
    --gff slyco/SlydLA2951_v2.0_gene_models_all_top12.gff3  \
    -o slyco/SlydLA2951_v2.0_gene_models_all_annot_top12.gtf 
wait
echo "DONE SINGULARITY LYCO"

## 1.3 S. pennellii
cat spen/spenn_v2.0_gene_models_annot.gff | grep -vE "Spenn-ch00" > spen/spenn_v2.0_gene_models_annot_top12.gff
singularity run -B $WORK/RNAseq/ agat_1.0.0.sif agat_convert_sp_gff2gtf.pl \
    --gff spen/spenn_v2.0_gene_models_annot_top12.gff \
    -o spen/spenn_v2.0_gene_models_annot_top12.gtf
wait
echo "DONE SINGULARITY SPEN"

## 1.4 S. pimpinellifolium 
cat spim/LA2093_v1.5.gff | grep -vE "SPIMPch00" | awk '($4 < $5){ print }' > spim/LA2093_top12.gff
singularity run -B $WORK/RNAseq/ agat_1.0.0.sif agat_convert_sp_gff2gtf.pl \
    --gff spim/LA2093_top12.gff \
    -o spim/LA2093_annot_top12.gtf
wait
echo "DONE SINGULARITY PIMP"

## 1.5 S. chilense
awk '$1 == "Scaffold_1__2_contigs__length_78070536" ||
     $1 == "Scaffold_2__1_contigs__length_61731525" ||
     $1 == "Scaffold_3__3_contigs__length_114446166" ||
     $1 == "Scaffold_4__2_contigs__length_60228994" ||
     $1 == "Scaffold_5__2_contigs__length_88504823" ||
     $1 == "Scaffold_6__2_contigs__length_83591330" ||
     $1 == "Scaffold_7__2_contigs__length_80977755" ||
     $1 == "Scaffold_8__2_contigs__length_76899493" ||
     $1 == "Scaffold_9__1_contigs__length_70652527" ||
     $1 == "Scaffold_10__2_contigs__length_69708562" ||
     $1 == "Scaffold_11__3_contigs__length_75158268" ||
     $1 == "Scaffold_12__3_contigs__length_71445296"' schil/S_chilense_Hirise.gff3 > schil/S_chilense_Hirise_top12.gff3
singularity run -B $WORK/RNAseq/ agat_1.0.0.sif agat_convert_sp_gff2gtf.pl \
    --gff schil/S_chilense_Hirise_top12.gff3 \
    -o schil/S_chilense_Hirise_annot_top12.gtf
wait
echo "DONE SINGULARITY SCHIL"

echo -e "\e[32mDONE SINGULARITY\e[0m"

# 2. Copy files to GeneExt folder
dirs=("schil" "shabr" "slyco" "spen" "spim")
for species in "${dirs[@]}"; do
    echo "$species"
    geneExtpath=GeneExt/"$species"
    if [ -d "$geneExtpath" ]; then
        echo "Directory $geneExtpath exists. Changing into it."
    else
        echo "Directory $geneExtpath does not exist. Creating it."
        mkdir -p "$geneExtpath"
    fi
    refpath=references/"$species"
    cp "$refpath"/*annot_top12.gtf "$geneExtpath"/
done
echo -e "\e[32mDONE COPY 2 GeneExt\e[0m"

# 3. Concatenate all .bam files into one to generate MASTER-bam with all reads for guided GTF-Extend 
cd mapped/clean/

conda activate samtools

for species in "${dirs[@]}"; do
    cd "$species"
    echo "$species"
    
    # Merge BAM files
    samtools merge -o all_samples.bam *_chr_rename.bam -f
    if [ $? -eq 0 ]; then
        echo "Successfully merged BAM files into all_samples.bam."
    else
        echo "Error: Failed to merge BAM files." >&2
        exit 1
    fi
    
    samtools index all_samples.bam
    
    # Index the merged BAM file
    if [ $? -eq 0 ]; then
        echo "Successfully indexed all_samples.bam."
    else
        echo "Error: Failed to index all_samples.bam." >&2
        exit 1
    fi
    
    # Filter and modify BAM headers
    samtools view -b all_samples.bam -L sclero_chromosomes.bed -U tmp_plant.bam -o sclero.bam
    if [ $? -eq 0 ]; then
        echo "Successfully filtered sclero.bam from all_samples.bam."
    else
        echo "Error: Failed to filter sclero.bam." >&2
        exit 1
    fi
    
    # Reheader
    samtools view -H all_samples.bam | sed '/SN:CP0178/d;/random/d;/chrUn/d' | samtools reheader - tmp_plant.bam > reads_sclero_depleted.bam
    if [ $? -eq 0 ]; then
        echo "Successfully reheadered BAM file."
    else
        echo "Error: Failed to reheader BAM file." >&2
        exit 1
    fi

    # Remove temp files
    rm tmp_plant.bam 
    rm sclero.bam
    
    # Check if the directory exists
    geneExtpath=GeneExt/"$species"
    if [ -d "$geneExtpath" ]; then
        echo "Directory $geneExtpath exists. Changing into it."
    else
        echo "Directory $geneExtpath does not exist. Creating it."
        mkdir -p "$geneExtpath"
    fi
    
    # Copy from Mapped to GeneExt dir
    cp reads_sclero_depleted.bam "$geneExtpath"/

    cd ..
done
wait
echo -e "\e[32mDONE BAM MERGING\e[0m"

# 4. Perform the actual GeneExt process
conda activate GeneExt
cd GeneExt/

# Many issues with non-indexed .bam files. --> adding samtools index in front. It would run through but just
# forget to continue after 'checking for exons'

for species in "${dirs[@]}"; do
    echo "$species"
    
    # Change into the directory
    cd "$species" || { echo "Failed to change directory to $species"; exit 1; }
    rm *_GeneExt.gtf # Just to be sure
    rm *.bai
    samtools index reads_sclero_depleted.bam # VERY MUCH NEEDED! OTHERWISE CRASH AFTER checking gene exons. 
    
    CURRENT_DATETIME=$(date +"%Y-%m-%d_%H-%M-%S")
    GTF=$(ls *top12.gtf)
    GTF_out="$species"_GeneExt.gtf
    # Use the date and time in a filename
    TMP="tmp_geneext${CURRENT_DATETIME}"
    # Run GeneExt
    python3 geneext.py -g "$GTF" -b reads_sclero_depleted.bam -o "$GTF_out" -t "$TMP" --peak_perc 10 --orphan -j 16 -v 1 -m 5000 --force
    
    if [ $? -eq 0 ]; then
        echo "GeneExt successfully completed for $species."
    else
        echo "Error: GeneExt failed for $species." >&2
        exit 1
    fi
    cd ..
done
echo -e "\e[32mDONE GENEeext\e[0m"
echo "GENE_Ext_pipeline_DONE"

