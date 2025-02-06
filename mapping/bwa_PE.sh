#!/bin/bash

# Define the path to the FASTQ files and genome index
fastq_directory="/path/to/fastq/files"
genome_index="/path/to/genome/fasta"
output_directory="/path/to/output"

# Index the genome (if it's not indexed already)
if [ ! -f "${genome_index}.bwt" ]; then
    echo "Indexing the genome..."
    bwa index "$genome_index"
fi

# Loop through all FASTQ pairs in the directory
for fastq_r1 in "$fastq_directory"/*_r1.fastq; do
    # Extract sample name (the part before the first underscore)
    sample_name=$(basename "$fastq_r1" | cut -d'_' -f1)

    # Get the corresponding R2 file
    fastq_r2="${fastq_r1/_r1.fastq/_r2.fastq}"

    echo -e "Aligning sample: $sample_name\n"

    # Perform BWA MEM alignment
    bwa mem -R "@RG\tID:$sample_name\tSM:$sample_name" "$genome_index" "$fastq_r1" "$fastq_r2" > "$output_directory/$sample_name.sam"

    echo -e "Finished aligning sample: $sample_name\n\n"
done
