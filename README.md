# 16S rRNA Amplicon Analysis of New York African Burial Ground Soil Samples
# By Carter Clinton, Ph.D.

## Overview

This repository contains the complete bioinformatic pipeline for 16S rRNA
amplicon sequencing analysis of burial soils from the New York African Burial
Ground (NYABG) and nearby urban control sites. The pipeline is end-to-end,
fully parameterized, and designed for reproducibility.

Key capabilities:
- QIIME 2 wrappers for demultiplexing, denoising, taxonomy, diversity, and differential abundance
- Region-specific SILVA classifier training and taxonomy export
- iTOL-ready phylogenetic tree generation with color annotations
- Human-associated and pathogenic taxa curation (Python)
- Ancient DNA (aDNA) damage profiling for authentication

## Publication

Clinton, C.K. and Jackson, F.L.C. (2025). Persistent human-associated
microbial signatures in burial soils from the 17th and 18th century New York
African Burial Ground. *ISME Communications*, 5(1).
DOI: [10.1093/ismeco/ycaf181](https://doi.org/10.1093/ismeco/ycaf181)

## Getting Started

### Prerequisites

- [QIIME 2](https://qiime2.org/) (tested with 2022.2 and amplicon-2024.10)
- [R](https://www.r-project.org/) (>=4.0) with packages: ALDEx2, phyloseq, tidyverse
- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/)

### Quick Setup

```bash
# 1. Clone this repository
git clone https://github.com/The-ASHES-Laboratory/ABG_16S.git
cd ABG_16S

# 2. Create QIIME 2 environment
mamba create -n qiime2-2024.5 -c qiime2 -c conda-forge qiime2 biom-format
conda activate qiime2-2024.5

# 3. (Optional) Create aDNA authentication environment
mamba create -n adna -c bioconda -c conda-forge bwa samtools damageprofiler
conda activate adna

# 4. (Optional) Python tools for pathogen curation
pip install pandas
```

## Pipeline Stages

The analysis is organized as 14 numbered scripts, run sequentially:

| Script | Description |
|--------|-------------|
| `01_demultiplex.sh` | Import and summarize demultiplexed paired-end reads |
| `02_denoise.sh` | DADA2 paired-end denoising with parameterized trimming |
| `03_seq_depth_and_QC.sh` | Feature table and rep-seq QC summaries for depth selection |
| `04_rarefaction.sh` | Alpha-rarefaction curves (observed features, Shannon, Faith's PD) |
| `05_taxonomy.sh` | SILVA 138 Naive Bayes classification (optional region-specific training) |
| `06_taxa_barplots.sh` | Taxa bar plots and rank-collapsed frequency tables |
| `07_build_phylogeny.sh` | Phylogenetic tree via MAFFT alignment + FastTree + midpoint rooting |
| `08_core_metrics_diversity.sh` | Alpha/beta diversity metrics, PERMANOVA, ANOSIM, pairwise tests |
| `09_permdisp.sh` | PERMDISP homogeneity of dispersion testing |
| `10_ancom.sh` | ANCOM compositional differential abundance analysis |
| `11_aldex2.R` | ALDEx2 differential abundance with CLR effect sizes (R) |
| `12_generate_itol_tree.sh` | iTOL-ready phylogenetic trees with color-strip annotations |
| `13_human_pathogen_curation.py` | Curate human-associated and pathogenic taxa from results (Python) |
| `14_adna_damage_profiler.sh` | Ancient DNA damage profiling for authentication (DamageProfiler) |

Every script accepts `--help` for full usage details, parameters, and examples.

### Example Workflow

```bash
# Activate QIIME 2
conda activate qiime2-2024.5

# Import and summarize reads
bash 01_demultiplex.sh -M manifest.csv -o 01_demux/

# Denoise with DADA2
bash 02_denoise.sh -i 01_demux/demux.qza -o 02_denoise/ --trunc-len-f 240 --trunc-len-r 220

# QC summary to choose rarefaction depth
bash 03_seq_depth_and_QC.sh -t 02_denoise/dada2_table.qza -s 02_denoise/dada2_rep-seqs.qza -o 03_qc/

# Rarefaction curves
bash 04_rarefaction.sh -i 02_denoise/dada2_table.qza -o 04_rare/ --max-depth 1000

# Taxonomy classification
bash 05_taxonomy.sh -t 02_denoise/dada2_table.qza -s 02_denoise/dada2_rep-seqs.qza \
  -c silva-138-classifier.qza -m metadata.tsv -o 05_taxonomy/

# Continue through scripts 06-14 as needed...
```

## Repository Structure

```
ABG_16S/
├── 01_demultiplex.sh              # Read import and demux summary
├── 02_denoise.sh                  # DADA2 denoising
├── 03_seq_depth_and_QC.sh         # Depth assessment and QC
├── 04_rarefaction.sh              # Alpha-rarefaction curves
├── 05_taxonomy.sh                 # SILVA taxonomy classification
├── 06_taxa_barplots.sh            # Taxa bar plots and collapsed tables
├── 07_build_phylogeny.sh          # Phylogenetic tree construction
├── 08_core_metrics_diversity.sh   # Core diversity metrics and tests
├── 09_permdisp.sh                 # PERMDISP dispersion testing
├── 10_ancom.sh                    # ANCOM differential abundance
├── 11_aldex2.R                    # ALDEx2 differential abundance (R)
├── 12_generate_itol_tree.sh       # iTOL tree preparation
├── 13_human_pathogen_curation.py  # Pathogen/human taxa curation (Python)
├── 14_adna_damage_profiler.sh     # aDNA damage profiling
├── LICENSE                        # MIT License
├── .gitignore                     # Standard Python/Conda ignores
└── README.md                      # This file
```

## Data Requirements (not included)

This repository does **not** include raw sequencing data. To reproduce the
analysis, you will need:

1. **16S rRNA amplicon reads** — 81 burial soil samples + 6 urban control
   samples, sequenced on Illumina MiSeq (V3-V4 region, paired-end)
2. **Sample metadata** — TSV with sample IDs, burial/control classification,
   and any grouping variables
3. **SILVA 138 classifier** — Pre-trained or train your own with script 05

Raw sequence data availability is described in the publication.

## Notes

- All scripts are fully parameterized via CLI arguments with sensible defaults.
  No hardcoded paths or system-specific configuration is required.

- Scripts capture the QIIME 2 environment (`qiime info`) in each output
  directory for reproducibility tracking.

- Default DADA2 parameters (truncation lengths 240/220, consensus chimera
  removal) reflect the settings used in the manuscript. Adjust for your
  sequencing run quality.

## Citation

If you use this repository, please cite:

> Clinton, C.K. and Jackson, F.L.C. (2025). Persistent human-associated
> microbial signatures in burial soils from the 17th and 18th century New York
> African Burial Ground. *ISME Communications*, 5(1).
> https://doi.org/10.1093/ismeco/ycaf181
