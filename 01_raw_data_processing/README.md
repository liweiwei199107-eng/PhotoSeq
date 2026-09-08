# PhotoSeq raw sequencing data processing

This directory contains the Linux raw sequencing data-processing workflow used for PhotoSeq.

PhotoSeq stands for **Photo-directed in situ Barcoding and Spatial Transcriptome Sequencing**. This workflow processes mouse PhotoSeq sequencing data from compressed FASTQ files to a UMI-deduplicated gene count matrix.

## 1. Reference data

PhotoSeq mouse reads were aligned to the GRCm39 genome using GENCODE mouse release M37.

```bash
mkdir -p reference/gencode_M37
cd reference/gencode_M37

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/GRCm39.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/gencode.vM37.annotation.gff3.gz

gzip -dk GRCm39.primary_assembly.genome.fa.gz
gzip -dk gencode.vM37.annotation.gff3.gz
cd ../..
```

## 2. Software environment

The upstream workflow was run on Linux with the following principal software versions:

- STAR v2.7.11b
- featureCounts v2.0.1 (distributed with Subread v2.0.1)
- UMI-tools v1.1.1

The analysis environment used Python 3.7.5.

```bash
conda create --name PhotoSeqEnv python=3.7.5 pandas matplotlib seaborn pytables biopython scikit-image
conda activate PhotoSeqEnv

conda install -c bioconda pysam umi_tools=1.1.1 subread=2.0.1 samtools=1.12
conda install -c biocore scikit-bio
```

Install STAR v2.7.11b separately. The official source archive can be installed as follows:

```bash
mkdir -p software
wget -O software/STAR-2.7.11b.tar.gz https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz
tar -xzf software/STAR-2.7.11b.tar.gz -C software
```

## 3. STAR genome index

The STAR index was generated with `--sjdbOverhang 249` using the GRCm39 primary assembly and GENCODE M37 annotation.

```bash
mkdir -p reference/STAR_GRCm39_M37

software/STAR-2.7.11b/source/STAR \
  --runThreadN 12 \
  --runMode genomeGenerate \
  --genomeDir reference/STAR_GRCm39_M37 \
  --genomeFastaFiles reference/gencode_M37/GRCm39.primary_assembly.genome.fa \
  --sjdbGTFfile reference/gencode_M37/gencode.vM37.annotation.gff3 \
  --sjdbOverhang 249
```

## 4. Directory structure

The workflow uses the following project layout:

```text
PhotoSeq_raw_processing/
├── README.md
├── environment.yml
├── raw_fastq/
│   └── <sample>_R1.fastq.gz
├── inFiles/
│   └── <motif-trimmed sample>_R1.fastq.gz
├── outFiles/
│   └── <all intermediate files and deduplicated BAM files>
├── outFilesR2/
│   └── <deduplicated BAM files for technical-replicate analysis>
├── outFilesS/
│   └── <deduplicated BAM files for Spatial analysis>
├── outFilesT/
│   └── <deduplicated BAM files for Temporal analysis>
├── ET/
├── S/
├── T/
├── reference/
│   ├── gencode_M37/
│   └── STAR_GRCm39_M37/
├── software/
│   └── STAR-2.7.11b/
└── scripts/
    ├── trim_photoseq_fastq_at_motifs.py
    ├── extract_photoseq_barcodes.py
    ├── align_photoseq_reads_to_mouse.py
    ├── assign_genes_and_deduplicate_photoseq_reads.py
    ├── generate_photoseq_technical_replicate_counts.py
    ├── generate_photoseq_spatial_counts.py
    └── generate_photoseq_temporal_counts.py
```

The scripts use the FASTQ naming pattern and sample identifiers described below.

## 5. PhotoSeq barcodes

Three 6-nt barcode-identifying segments were used:

| Barcode | Sequence |
| --- | --- |
| barcode1 | `GTTAGG` |
| barcode2 | `AGGGTA` |
| barcode3 | `TATGGA` |

In the UMI-tools extraction regular expression, `cell_1` is the reserved group that transfers the PhotoSeq ROI barcode into the read identifier. The group name does not imply that each barcode represents a single cell.

In the Spatial dataset, barcode1 identifies the tumour ROI, barcode2 identifies the normal-adjacent ROI and barcode3 identifies the normal ROI. In the Temporal dataset, the barcode-to-sample mappings are encoded explicitly in `generate_photoseq_temporal_counts.py`.

## 6. Run order

Run motif trimming in the directory containing the raw FASTQ files, then place the resulting FASTQ files in `inFiles/`. Run all subsequent scripts from the project root.

All samples are processed together through barcode extraction, alignment, gene assignment and UMI deduplication. These steps use the common `outFiles/` directory. Only after deduplication are the required `*_Dedup.bam` files copied into the three analysis-specific directories used to generate count matrices.

```bash
conda activate PhotoSeqEnv

cd raw_fastq
python ../scripts/trim_photoseq_fastq_at_motifs.py
cd ..

# Place the motif-trimmed FASTQ files in inFiles/, then continue:
python scripts/extract_photoseq_barcodes.py
python scripts/align_photoseq_reads_to_mouse.py
python scripts/assign_genes_and_deduplicate_photoseq_reads.py

# Copy deduplicated BAM files into the analysis-specific input directories.
# Copy rather than move files because some samples are used in more than one analysis.
mkdir -p outFilesR2 outFilesS outFilesT ET S T
cp outFiles/{EL1rm,E1Trm}_Dedup.bam outFilesR2/
cp outFiles/{TAN1rm,TAN2rm,TAN3rm,MAN1rm,MAN2rm,MAN3rm}_Dedup.bam outFilesS/
cp outFiles/{EL1rm,EL2rm,EL3rm,TAN1rm,TAN2rm,TAN3rm,MAN1rm,MAN2rm,MAN3rm}_Dedup.bam outFilesT/

# Generate the three count matrices from their corresponding directories:
python scripts/generate_photoseq_technical_replicate_counts.py
python scripts/generate_photoseq_spatial_counts.py
python scripts/generate_photoseq_temporal_counts.py
```

The count-generation inputs are grouped as follows:

| Directory | Analysis | Deduplicated BAM files |
| --- | --- | --- |
| `outFilesR2/` | Technical replicates | `EL1rm`, `E1Trm` |
| `outFilesS/` | Spatial | `TAN1rm`, `TAN2rm`, `TAN3rm`, `MAN1rm`, `MAN2rm`, `MAN3rm` |
| `outFilesT/` | Temporal | `EL1rm`, `EL2rm`, `EL3rm`, `TAN1rm`, `TAN2rm`, `TAN3rm`, `MAN1rm`, `MAN2rm`, `MAN3rm` |

Each name in the table denotes the corresponding `<name>_Dedup.bam` file. The technical-replicate, Spatial and Temporal count-generation scripts read `outFilesR2/`, `outFilesS/` and `outFilesT/`, respectively. The upstream processing scripts do not write directly to these analysis-specific directories.

The technical-replicate output names are assigned as follows:

| BAM/barcode input | Output column |
| --- | --- |
| `EL1rm/TATGGA` | `Breast_early_N1` |
| `E1Trm/TATGGA` | `Breast_early_N1_techrep2` |

Both columns describe the same biological sample; the second is its technical replicate.

Expected stages:

1. trim reads at the established PhotoSeq motifs;
2. extract the UMI and one of the three PhotoSeq ROI barcodes;
3. align the processed mouse reads to the GRCm39/M37 STAR index;
4. assign aligned reads to genes with featureCounts;
5. deduplicate reads by UMI and gene using UMI-tools;
6. copy the deduplicated BAM files from `outFiles/` into the applicable analysis-specific directories;
7. generate gene-by-sample/barcode count matrices.

The resulting matrix is a **UMI-deduplicated gene count matrix** (also termed a **UMI-collapsed gene count matrix**).

`extract_photoseq_barcodes.py` recreates the common `outFiles/` directory at the start of a run, replacing intermediate outputs from a preceding run.

## 7. Analysis parameters

The barcode-extraction pattern discards the initial 27 nt, extracts a 12-nt UMI and identifies one of the three 6-nt PhotoSeq barcodes. Reads are aligned to the mouse reference with STAR using unique-only mapping. Gene assignment is performed with featureCounts, followed by per-gene UMI deduplication with UMI-tools. The gene-assignment workflow excludes `ENSMUSG00000136525.1` and `ENSMUSG00000119584.1`, as encoded in `assign_genes_and_deduplicate_photoseq_reads.py`.

## 8. Validation

The Python scripts can be syntax-checked with:

```bash
python -m py_compile scripts/*.py
```

## References

- [GENCODE mouse release M37 (GRCm39)](https://www.gencodegenes.org/mouse/release_M37.html)
- [STAR releases](https://github.com/alexdobin/STAR/releases)
- [UMI-tools documentation](https://umi-tools.readthedocs.io/en/latest/)
- [Subread and featureCounts](https://subread.sourceforge.net/)
