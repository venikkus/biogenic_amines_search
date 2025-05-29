# Biogenic amines search in bacterial genomes

This repository contains a **Snakemake**-based pipeline for searching for biogenic amines in bacterial genomes, which performs genome assembly and annotation of prokaryotic genomes.

## Research overview

Biogenic amines (histamine, tyramine) can cause poisoning and allergies.

The consumption of foods containing high concentrations of biogenic amines has been associated with health hazards [1].

## Research objectives

The pipeline will allow testing the genetic data of bacterial strains for the presence of genes responsible for the synthesis of biogenic amines. The pipeline screens bacterial genome assemblies or raw reads to detect the following genes:


## Relevance of the study

Biogenic amines can cause poisoning and allergies. The consumption of foods containing high concentrations of biogenic amines has been associated with health hazards [1]. 

## Content

- [Pipeline Overview](#pipeline-overview)
- [Requirements](#requirements)
- [Install](#install)
- [Running instructions](#running-instructions)
- [Examples](#examples)
- [Contact](#contact)

## Pipeline Overview

The pipeline consists of the following steps:

1. **Read Downloading (optional)**  
   - Download reads from the SRA database (`prefetch` + `fasterq-dump`).
   - Copy reads from local storage if `SOURCE=offline`.

2. **Read Trimming**  
   - Clean paired-end reads with `Trimmomatic` to remove low-quality bases.

3. **Genome Assembly**  
   - Use `SPAdes` to assemble reads into contigs.

4. **Assembly Downloading (optional)**  
   - Download existing assemblies from online sources or copy from local storage.

5. **Genome Annotation**  
   - Annotate the assembly with `Prokka` to identify genes and other features.

6. **Assembly Quality Assessment**  
   - Evaluate assembly quality with `QUAST`.

## System requirements

- Linux (tested on Ubuntu 22.04)

## Install

1. Clone this repository to your local machine:
```bash
git clone https://github.com/venikkus/biogenic_amines_search.git
cd biogenic_amines_search
```

2. Use the provided environment.yaml file to create the conda environment.  If you don't have conda, use these installation instuctions:

```bash
conda env create -f environment.yaml
conda activate amines_search
```


Running instructions
The pipeline is configured via config.yaml, where you can specify:

>`source` (online / offline)

>`input_type` (sra / assembly)

>`id` (NCBI accession or sample ID)

Run the full pipeline
```bash
snakemake -p --cores [NUM_CORES] --config id=[ID] input_type=sra source=[SOURCE]

```

- [NUM_CORES] — number of CPU cores to use (e.g., 4, 16, etc.)
- [ID] - SRA ID or assembly ID.
- [SOURCE] -- offline (from folder) or online (from )


## Features
Output structure:

1. Input Reads
input_reads/{id}_1.fastq, input_reads/{id}_2.fastq

2. Clean Reads
clean_data/{id}_1P.fastq, clean_data/{id}_2P.fastq

3. Assembly
assembly/{id}/contigs.fasta

4. Annotation
prokka_output/{id}/
Prokka-generated GFF, GBK, and other annotation files.

5. Assembly Quality
quast_output/{id}/
QUAST reports, including report.html and statistics.

## References

prokka 1.14.6
snakemake 9.5.1
sra-toolkit 2.11.3
SPAdes genome assembler v4.0.0
QUAST v5.3.0

## Contact

Please report any problems directly to the GitHub [issue tracker](https://github.com/venikkus/biogenic_amines_search/issues).

Also, you can send your feedback to authors:
- [niksamusik@gmail.com](mailto:niksamusik@gmail.com)