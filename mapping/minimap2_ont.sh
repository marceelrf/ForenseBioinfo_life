#!/bin/bash

# Define the path to the FASTQ files and reference genome
fastq_directory="/path/to/fastq/files"
genome_reference="/path/to/genome/fasta"
output_directory="/path/to/output"
threads=16  # Set the number of threads to use (adjust as needed)

# Loop through each FASTQ file in the directory
for fastq_file in "$fastq_directory"/*.fastq; do
    # Extract sample name (the part before the first underscore or hyphen)
    sample_name=$(basename "$fastq_file" | cut -d'_' -f1 | cut -d'-' -f1)

    echo -e "Aligning sample: $sample_name\n"

    # Perform minimap2 alignment
    minimap2 -ax map-ont "$genome_reference" "$fastq_file" -R "@RG\tID:$sample_name\tSM:$sample_name" -t "$threads" > "$output_directory/$sample_name.sam"

    echo -e "Finished aligning sample: $sample_name\n\n"
done
