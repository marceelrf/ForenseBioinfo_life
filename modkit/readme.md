# About

This repository contains Bash scripts for automating DNA methylation analysis using modkit. The scripts process multiple samples efficiently by iterating over directories of BAM files and generating corresponding methylation-related outputs.
Pipeline Overview

The analysis includes four key steps:

    Pileup Analysis (run_pileup.sh)

        Computes methylation pileup from BAM files.

        Outputs BED files with methylation data.

    Haplotype-Specific Methylation (run_haplotype_pileup.sh)

        Generates haplotype-resolved methylation data.

        Saves results with haplotype-specific prefixes.

    Differential Methylation Regions (DMR) (run_dmrhap.sh)

        Compares haplotype-specific methylation levels.

        Identifies differentially methylated regions.

    Methylation Entropy Analysis (run_methyl_entropy.sh)

        Measures methylation entropy per gene and individual.
