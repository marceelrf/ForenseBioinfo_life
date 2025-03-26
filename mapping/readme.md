# Mapping

This repository provides a bioinformatics workflow for sequence alignment using BWA (for Illumina short reads) and Minimap2 (for Nanopore long reads), followed by processing steps with SAMtools. The included diagram illustrates the overall pipeline, starting from raw sequencing reads (.fastq files) and ending with sorted and indexed BAM files for downstream analysis.

## Workflow Overview
1. Input Data:
- Illumina paired-end reads (R1.fastq and R2.fastq);
- Nanopore single-end reads (.fastq);
- Reference genome (.fasta)

2. Alignment:
- BWA-MEM for short-read mapping;
- Minimap2 for long-read mapping;

3. Post-processing with SAMtools:
- Convert SAM to BAM format;
- Sort BAM files;
- Index BAM files

This workflow ensures a standardized pipeline for processing sequencing data, making it easier to analyze genomic data efficiently. Feel free to explore, modify, and contribute!

![mapping](https://github.com/marceelrf/ForenseBioinfo_life/blob/main/mapping/mermaid-diagram-2025-02-07-093838.png)
