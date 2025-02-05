# About Clair3

Clair3 is a deep-learning-based variant caller optimized for third-generation sequencing technologies such as Oxford Nanopore and PacBio. It combines neural networks with heuristic methods to accurately detect single nucleotide variants (SNVs) and small insertions/deletions (Indels). Clair3 is efficient and scalable, making it suitable for large-scale variant calling.  
GLnexus is a variant call aggregation and refinement tool for gVCF data. It merges variant calls from multiple samples into a unified set, optimizing the accuracy and consistency of the final variants. GLnexus is widely used in genomic projects involving multi-sample genotyping.

Both tools are commonly used together in bioinformatics pipelines for variant calling and aggregation. 🚀

## Variant Calling Pipeline with Clair3 and GLnexus
This repository contains scripts and a workflow for calling variants from sequencing data using Clair3 and GLnexus. The pipeline processes .BAM files to generate variant files in VCF format in two main steps.

### Workflow
1. Input: BAM file containing aligned reads.
2. Step 1:
  - The Clair3 variant calling model processes the BAM and generates .gVCF files.
  - GLnexus aggregates the .gVCF files to generate an initial list of VCF candidates.
3. Step 2:
  - The candidate variants undergo a second processing using GLnexus to generate the final VCF.
### Scripts
`run_mult_clair3.py`: Runs the first stage of the pipeline, generating .gVCF files.  
`run_mult_clair3_step2.py`: Runs the second step, refining the variant calls to produce the final VCF.

![Clair3 flowchat](https://github.com/marceelrf/ForenseBioinfo_life/blob/main/clair3/mermaid-diagram-2025-02-05-145122.png)
