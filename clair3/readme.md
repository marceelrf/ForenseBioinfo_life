# About Clair3

![Clair3 flowchat](https://github.com/marceelrf/ForenseBioinfo_life/blob/main/clair3/mermaid-diagram-2025-02-05-145122.png)
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
run_mult_clair3.py: Runs the first stage of the pipeline, generating .gVCF files.
run_mult_clair3_step2.py: Runs the second step, refining the variant calls to produce the final VCF.
