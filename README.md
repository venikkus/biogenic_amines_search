# Biogenic amines search in bacterial genomes

This repository contains a **Snakemake**-based pipeline for searching for biogenic amines in bacterial genomes, which performs genome assembly and annotation of prokaryotic genomes.

## Content

- [Research overview](#research-overview)
- [Research objectives](#research-objectives)
- [Features](#key-features)
- [Pipeline Overview](#pipeline-overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Run the pipeline](#run-the-pipeline)
- [Configuration](#configuration)
- [Output structure](#output-structure)
- [Contact](#contact)

## Research overview

Biogenic amines (histamine, tyramine) can cause poisoning and allergies.

The consumption of foods containing high concentrations of biogenic amines has been associated with health hazards [1].

## Research objectives

The pipeline will allow testing the genetic data of bacterial strains for the presence of genes responsible for the synthesis of biogenic amines. The pipeline screens bacterial genome assemblies or raw reads to detect the following genes:

| Genome ID | Organism | Gene | Source | Assembly Link |
|-----------------|---------------------------------------------------|-------|--------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------|
| GCF_000350465.1 | *Enterococcus durans* IPLA 655                    | hdcA  | [Ladero et al. (2013)](https://pubmed.ncbi.nlm.nih.gov/23682153/)                                                                                                                       | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000350465.1/)            |
| GCF_000009685.1 | *Clostridium perfringens* str. 13                 | hdcA  | [Landete et al. (2008)](https://www.tandfonline.com/doi/full/10.1080/10408390701639041)                                                                | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000009685.1/)            |
| GCF_901482635.1 | *Staphylococcus capitis*                          | hdcA  | [Landete et al. (2008)](https://www.tandfonline.com/doi/full/10.1080/10408390701639041)                                                                | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_901482635.1/)            |
| GCF_000159435.1 | *Limosilactobacillus vaginalis* DSM 5837          | hdcA  | [Landete et al. (2015)](https://www.sciencedirect.com/science/article/pii/S0168160515301136)                                                           | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000159435.1/)            |
| GCF_000317165.1 | *Ligilactobacillus saerimneri* 30a                | hdcA  | [Landete et al. (2015)](https://www.sciencedirect.com/science/article/pii/S0168160515301136)                                                           | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000317165.1/)            |
| GCF_001677035.1 | *Lentilactobacillus parabuchneri*                 | hdc   | [Schirone et al. (2017)](https://ifst.onlinelibrary.wiley.com/doi/full/10.1111/ijfs.13826)                                                            | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001677035.1/)            |
| GCF_001687155.1 | *Lentilactobacillus parabuchneri*                 | hdc   | [Schirone et al. (2017)](https://ifst.onlinelibrary.wiley.com/doi/full/10.1111/ijfs.13826)                                                            | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001687155.1/)            |
| GCF_001687145.1 | *Lentilactobacillus parabuchneri*                 | hdc   | [Schirone et al. (2017)](https://ifst.onlinelibrary.wiley.com/doi/full/10.1111/ijfs.13826)                                                            | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001687145.1/)            |
| GCF_900258435.1 | *Carnobacterium divergens*                        | tdc   | [Marcobal et al. (2012)](https://www.tandfonline.com/doi/abs/10.1080/10408398.2010.500545)                                                            | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_900258435.1/)            |
| GCF_003966405.1 | *Enterococcus faecalis*                           | tdc   | [Marcobal et al. (2012)](https://www.tandfonline.com/doi/abs/10.1080/10408398.2010.500545)                                                            | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003966405.1/)            |
| GCF_008365415.1 | *Enterococcus faecium*                            | tdc   | [Marcobal et al. (2012)](https://www.tandfonline.com/doi/abs/10.1080/10408398.2010.500545)                                                            | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_008365415.1/)            |
| GCF_000026105.1 | *Pseudomonas entomophila* L48                     | tdc   | [Marcobal et al. (2012)](https://www.tandfonline.com/doi/abs/10.1080/10408398.2010.500545)                                                            | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000026105.1/)            |
| GCF_000007565.2 | *Pseudomonas putida* KT2440                       | tdc   | [Marcobal et al. (2012)](https://www.tandfonline.com/doi/abs/10.1080/10408398.2010.500545)                                                            | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000007565.2/)            |
| GCF_000014285.2 | *Granulibacter bethesdensis* CGDNIH1              | tdc   | [Marcobal et al. (2012)](https://www.tandfonline.com/doi/abs/10.1080/10408398.2010.500545)                                                            | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000014285.2/)            |
| GCF_029581645.1 | *Latilactobacillus curvatus*                      | tdcA  | [La Gioia et al. (2012)](https://journals.asm.org/doi/full/10.1128/aem.01928-10)                                                                                                                        | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_029581645.1/)            |
| GCF_000014465.1 | *Levilactobacillus brevis* ATCC 367               | agdi, tdc | [Lucas et al. (2007)](https://www.sciencedirect.com/science/article/pii/S074000200800138X)<br>[Lucas et al. (2014)](https://www.sciencedirect.com/science/article/pii/S0168160514000373) | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000014465.1/)            |
| GCF_000161935.1 | *Limosilactobacillus coleohominis* 101-4-CHN      | tdc   | [Bonnin-Jusserand et al. (2012)](https://link.springer.com/article/10.1186/1471-2180-12-199)                                                           | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000161935.1/)            |
| GCF_000180015.1 | *Limosilactobacillus oris* PB013-T2-3             | tdc   | [Bonnin-Jusserand et al. (2012)](https://link.springer.com/article/10.1186/1471-2180-12-199)                                                           | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000180015.1/)            |
| GCF_000008065.1 | *Lactobacillus johnsonii* NCC 533                 | odc   | [Marcobal et al. (2004)](https://academic.oup.com/femsle/article-abstract/239/2/213/622022)<br>[Lucas et al. (2014)](https://www.sciencedirect.com/science/article/pii/S0168160514000373) | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000008065.1/)            |
| GCF_000011985.1 | *Lactobacillus acidophilus* NCFM                  | odc   | [Lucas et al. (2014)](https://www.sciencedirect.com/science/article/pii/S0168160514000373)                                                             | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000011985.1/)            |
| GCF_000014425.1 | *Lactobacillus gasseri* ATCC 33323                | odc   | [Lucas et al. (2014)](https://www.sciencedirect.com/science/article/pii/S0168160514000373)                                                             | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000014425.1/)            |
| GCF_000014505.1 | *Pediococcus pentosaceus* ATCC 25745              | agdi  | [Lucas et al. (2014)](https://www.sciencedirect.com/science/article/pii/S0168160514000373)                                                             | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000014505.1/)            |
| GCF_000008285.1 | *Listeria monocytogenes* serotype 4b str. F2365   | agdi  | [Lucas et al. (2014)](https://www.sciencedirect.com/science/article/pii/S0168160514000373)                                                             | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000008285.1/)            |
| GCF_000026065.1 | *Latilactobacillus sakei* subsp. *sakei* 23K      | agdi  | [Lucas et al. (2014)](https://www.sciencedirect.com/science/article/pii/S0168160514000373)                                                             | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000026065.1/)            |

## Key Features

**Dual input support:**
- SRA/ENA reads (online download)  
- Local FASTQ files  
- Public assemblies (NCBI)  
- Local assemblies  

**Comprehensive analysis:**
- Quality control (Trimmomatic)  
- Genome assembly (SPAdes)  
- Annotation (Prokka)  
- Targeted gene search  
- Quality metrics (QUAST)  

**Interactive reporting:**
- Sample-level HTML reports  
- Presence/absence matrix  
- Heatmap visualization  


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

## Installation

### Conda (recommended)
```bash
git clone https://github.com/venikkus/biogenic_amines_search.git
cd biogenic_amines_search
conda env create -f environment.yaml
conda activate amines_search
```


## Quick Start

### Edit the configuration file:

- [NUM_CORES] — number of CPU cores to use (e.g., 4, 16, etc.)
- [ID] — SRA ID or assembly ID.
- [SOURCE] — offline (from folder) or online (from NCBI)

```bash
# config.yaml
assembly_samples:
  id: ["GCF_000350465.1", "GCF_000009685.1"]  # NCBI Assembly IDs
  source: "online"  # or "offline" for local files

sra_samples:
  id: ["SRR123456", "SRR789012"]  # SRA accessions
  source: "online"  # or "offline"
```

### Run the pipeline:

```bash
snakemake --cores 8 --use-conda
```


## Configuration

### Input Options
|Parameter	|Description	|Example |Values|
|---|---|---|--|
|`assembly_samples`	|Pre-assembled genomes	|NCBI accession| IDs|
|`sra_samples`	|Raw sequencing reads|	SRA/ENA accessions
|`source`	|Data |location	|"online"/"offline"


## Output structure

```
results/
├── input_reads/               # SRA
├── assembly/               # Final assemblies
├── quast_output/           # Quality metrics
├── prokka_output/          # Annotation files
├── decarboxylases/         # Detected genes
└── reports/
    ├── sample_*.html       # Individual reports
    └── full_report.html    # Summary dashboard
```

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