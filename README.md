# Biogenic amines search in bacterial genomes

This repository contains a **Snakemake**-based pipeline for searching for biogenic amines in bacterial genomes, which performs genome assembly and annotation of prokaryotic genomes.

## Content

- [Research objectives](#research-objectives)
- [Features](#key-features)
- [Pipeline Overview](#pipeline-overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Run the pipeline](#run-the-pipeline)
- [Configuration](#configuration)
- [Output structure](#output-structure)
- [Contact](#contact)

## Research objectives

The pipeline will allow testing the genetic data of bacterial strains for the presence of genes responsible for the synthesis of biogenic amines. We used these genomes:

| Genome ID | Organism | Gene (mentioned in article) | Source | Assembly Link |
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


The pipeline screens bacterial genome assemblies or raw reads to detect the following genes:

| Name               | Decarboxylase Name | Gene(s) | Classification | Chemical Structure | Precursor Amino Acid | References |
|--------------------|--------------------|---------|----------------|--------------------|----------------------|------------|
| Methylamine        | Glycine decarboxylase | gcv, gdc | Monoamine | Aliphatic | Glycine | [1](https://www.cell.com/trends/plant-science/fulltext/S1360-1385(01)01892-1), [2](https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1365-294X.1997.00279.x), [3](https://www.sciencedirect.com/topics/medicine-and-dentistry/glycine-dehydrogenase-decarboxylating) |
| Ethylamine         | Alanine decarboxylase | aldc | Monoamine | Aliphatic | Alanine | [1](https://www.sciencedirect.com/science/article/abs/pii/S0031942200941736), [2](https://link.springer.com/article/10.1186/s12896-021-00674-x) |
| Phenylethylamine   | Tyrosine decarboxylase | tdc/DDC, aadc | Monoamine | Aromatic | Phenylalanine | [1](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2015.00259/full), [2](https://www.tandfonline.com/doi/abs/10.1080/10408398.2010.500545) |
| Tyramine           | Tyrosine decarboxylase | tdc, tdcA, tdcB, TDC-1, DDC, TYDC, tyrDC, mfnA, adc | Monoamine | Aromatic | Tyrosine | [1](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2015.00259/full), [2](https://www.tandfonline.com/doi/abs/10.1080/10408398.2010.500545) |
| Octopamine         | Tyrosine decarboxylase | tdc | Monoamine | Aromatic | Tyrosine (Tyramine intermediate) | [1](https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/octopamine), [2](https://www.annualreviews.org/content/journals/10.1146/annurev.ento.50.071803.130404) |
| Dopamine           | Tyrosine decarboxylase, L-tryptophan decarboxylase | TYDC, DDC, TDC | Monoamine | Aromatic | Tyrosine | [1](https://www.annualreviews.org/content/journals/10.1146/annurev.ento.50.071803.130404) |
| Norepinephrine     | Dopamine beta-monooxygenase | - | Monoamine | Aromatic | Tyrosine | - |
| Adrenalin          | Phenylethanolamine N-methyltransferase | - | Monoamine | Aromatic | Tyrosine | - |
| Histamine          | Histidine decarboxylase | hdc | Monoamine | Heterocyclic | Histidine | [1](https://journals.asm.org/doi/full/10.1128/aem.01496-07) |
| Tryptamine         | Tryptophan decarboxylase | tdc | Monoamine | Heterocyclic | Tryptophan | [1](https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1365-313X.1997.11061167.x) |
| Serotonin          | Aromatic-L-amino-acid/L-tryptophan decarboxylase, Tryptamine 5-hydroxylase | DDC, TDC, CYP71P1 | Monoamine | Heterocyclic | Hydroxytryptophane | - |
| Putrescine         | Ornithine decarboxylase | odc, speC, speF | Diamine | Aliphatic | Ornithine | [1](https://www.jbc.org/article/S0021-9258(20)72392-6/fulltext), [2](https://journals.asm.org/doi/abs/10.1128/jb.175.5.1221-1234.1993), [3](https://www.sciencedirect.com/science/article/pii/S0362028X22064948) |
|                    | Arginine decarboxylase | adc, adiA, speA | Polyamine | Aliphatic | Arginine | [1](https://www.jbc.org/article/S0021-9258(20)72392-6/fulltext), [2](https://academic.oup.com/mbe/article/15/10/1312/999792), [3](https://journals.asm.org/doi/abs/10.1128/jb.175.5.1221-1234.1993) |
| Agmatine           | Arginine decarboxylase | speA, adiA | Polyamine | Aliphatic | Arginine | [1](https://nyaspubs.onlinelibrary.wiley.com/doi/full/10.1196/annals.1304.004), [2](https://pmc.ncbi.nlm.nih.gov/articles/PMC9241758/), [3](https://pubmed.ncbi.nlm.nih.gov/39546898/) |
| Cadaverine         | Lysine decarboxylase | ldcC, cadA | Diamine | Aliphatic | Lysine | [1](https://pubs.acs.org/doi/pdf/10.1021/bi00701a005), [2](https://www.microbiologyresearch.org/content/journal/micro/10.1099/00221287-144-3-751), [3](https://journals.asm.org/doi/abs/10.1128/jb.175.5.1221-1234.1993) |
|                    | Cadaverine antiporter | cadB | - | - | - | [1](https://academic.oup.com/jimb/article/41/4/701/5995075) |
| Spermidine         | Spermidine synthase | speE, speD | Polyamine | Aliphatic | Arginine/Ornithine | [1](https://pmc.ncbi.nlm.nih.gov/articles/PMC9241758/) |
| Spermine           | Spermine synthase | speE, speD | Polyamine | Aliphatic | Arginine/Ornithine | [1](https://pmc.ncbi.nlm.nih.gov/articles/PMC11107197/) |

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
- Sample-level HTML reports and summary report
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

Use the provided `environment.yaml` file to create the conda environment. If you don't have conda, use these [installation instuctions](https://www.anaconda.com/docs/getting-started/miniconda/install#linux):
```bash
git clone https://github.com/venikkus/biogenic_amines_search.git
cd biogenic_amines_search
conda env create -f environment.yml
conda activate amines_search
```

## Quick Start

### Edit the configuration file:

Example:
```bash
# For ONLINE data (NCBI downloads)
assembly_samples:
  id: 
    - GCF_000350465.1  # Enterococcus durans IPLA 655
    - GCF_000009685.1  # Clostridium perfringens str.13
  source: "online"

sra_samples:
  id:
    - SRR8594964  # Lactobacillus acidophilus
  source: "online"

# For OFFLINE data (local files)
assembly_samples:
  id:
    - offline_GCF_901482635.1  # Local assembly
  source: "offline"
  local_path: "/path/to/assemblies"  # Directory containing *.fasta files

sra_samples:
  id:
    - offline_SRR8594964  # Local reads
  source: "offline" 
  local_path: "/path/to/reads"  # Directory containing *_1.fastq and *_2.fastq

genes_of_interest: "genes.txt"  # Optional gene list
```

### Run the pipeline:

```bash
snakemake -p --cores [NUM_CORES]
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