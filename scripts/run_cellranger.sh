#!/bin/bash

# Path to STAR index
GENOME_DIR="/home/raquelcr/ratopin_genome/HetGla_female_1.0"

# Path to sample sheet
SAMPLE_SHEET="/home/raquelcr/cellranger_output/sample_list_cellranger.tsv"

# Read the sample sheet line by line (skip header)
tail -n +2 "$SAMPLE_SHEET" | while IFS=$'\t' read -r Sample IDs; do
    echo "Processing sample: $Sample"
    # Create output directory for the sample
    mkdir -p "out_HetGla_female_1/$Sample"


    # Run cellranger
    cellranger count --id $Sample \
    --transcriptome /home/raquelcr/ratopin_genome/HetGla_female_1 \
    --fastqs /home/raquelcr/jcgorozco/scRNAseq/ \
   --sample $IDs \
   --create-bam false \
   --localcores 8 \
   --localmem 4 \
   --expect-cells 10000 \
   --nosecondary 
done