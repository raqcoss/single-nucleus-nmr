#!/bin/bash

# Define the directory containing your fastq files
INPUT_DIR="/home/raquelcr/jcgorozco/scRNAseq"
OUTPUT_DIR="/home/raquelcr/fastqc_output"

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through all fastq files in the input directory

for FILE in $INPUT_DIR/20230810scRNAAjoRat04*/*.fastq.gz; do
    # Run FastQC on each file
    fastqc $FILE -o $OUTPUT_DIR
done
for FILE in $INPUT_DIR/20230810scRNAAjoRat05*/*.fastq.gz; do
    # Run FastQC on each file
    fastqc $FILE -o $OUTPUT_DIR
done

echo "FastQC analyses are complete. Check the output directory for results."