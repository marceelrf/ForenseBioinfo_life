#!/bin/bash

# Define the path to the BAM files directory and the BED file
bam_directory="/path/to/bam/"
bed_file="/path/to/BEDs/file.bed"
output_directory="/path/to/output"
threads=20 #change here

# Loop through each BAM file in the directory
for bam_file in "$bam_directory"/*.bam; do
    # Extract the sample name (first part before the hyphen)
    sample_name=$(basename "$bam_file" | cut -d'-' -f1)
    
    echo -e "Começando a amostra: $sample_name\n"
    
    # Perform the samtools view command
    samtools view -@ $threads -b -L "$bed_file" "$bam_file" > "$output_directory/$sample_name.bam"
    
    # Add read group information with samtools addreplacerg
    samtools addreplacerg -@ $threads -r "ID:$sample_name" -r "SM:$sample_name" -o "$output_directory/tmp.bam" "$output_directory/$sample_name.bam"
    
    # Sort the BAM file
    samtools sort -@ $threads -o "$output_directory/sorted.bam" "$output_directory/tmp.bam"
    
    # Index the sorted BAM file
    samtools index -@ $threads "$output_directory/sorted.bam"
    
    # Rename the sorted BAM file to the sample name
    mv "$output_directory/sorted.bam" "$output_directory/$sample_name.bam"
    mv "$output_directory/sorted.bam.bai" "$output_directory/$sample_name.bam.bai"
    
    # Remove the temporary BAM file
    rm "$output_directory/tmp.bam"
    
    echo -e "Finalizando a amostra: $sample_name\n\n"
done
