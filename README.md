# PhotoSeq analysis code

This repository contains the code and compact analysis inputs used for PhotoSeq, from raw mouse FASTQ processing through downstream statistical analysis and figure generation.

## Repository layout

```text
.
├── 01_raw_data_processing/       # FASTQ to UMI-deduplicated gene counts (Linux)
└── 02_downstream_analysis/       # Statistical and bioinformatic analyses (Windows/R and Python)
```

The two stages were run on different operating systems. This is intentional and does not affect the analysis: each stage has its own environment record, and all repository paths are relative.

## Analysis overview

1. Download the GRCm39 primary assembly and GENCODE mouse M37 annotation.
2. Build the STAR index for the GRCm39 primary assembly and GENCODE mouse M37 annotation.
3. Trim the established motifs from the raw FASTQ files.
4. Extract the 12-nt UMI and identify the PhotoSeq ROI barcode.
5. Align reads uniquely to the mouse genome, assign reads to genes and deduplicate by UMI and gene.
6. Generate the technical-replicate, Spatial and Temporal gene-count matrices.
7. Prepare the downstream matrices and metadata, then run the numbered R scripts.
8. In script 12, identify the five-way DEG intersection, reconstruct the frozen STRING PPI network and generate the PPI-strength, MCC, Degree and MCODE panels.

Detailed commands for steps 1–6 are in [01_raw_data_processing/README.md](01_raw_data_processing/README.md). Instructions for steps 7 onward are in [02_downstream_analysis/README.md](02_downstream_analysis/README.md).

## Included data

The repository includes the count matrices, sample metadata and compact reference tables needed by the supplied downstream scripts. This includes the frozen STRING response used for the reported PPI network. Raw FASTQ files, genome/annotation files, STAR indices, BAM files and other generated analysis outputs are not included because they are large or reproducible from the documented workflow.

Generated figures and result tables are written to numbered output folders by the analysis scripts and are not versioned in this repository.

## Technical-replicate naming

The technical-replicate pair represents two measurements of the same biological sample:

| Count-matrix column | Meaning |
| --- | --- |
| `Breast_early_N1` | first measurement |
| `Breast_early_N1_techrep2` | second technical measurement of `Breast_early_N1` |

These names are used consistently in the upstream count-generation script, the supplied count matrix and the downstream quality-control script.

## Code availability statement

The custom scripts used for raw sequencing data processing and downstream analysis are provided in this repository. The repository also contains the compact input and reference tables required to reproduce the reported analyses; large raw sequencing files are made available through the data repository specified in the manuscript.

## Reproducibility note

Run each stage from the directory specified in its README. The code records the transformations and statistical thresholds used by the analyses.
