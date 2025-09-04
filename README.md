# 16S rRNA Analysis of New York African Burial Ground Soil Samples

### Description
Microbiome Workflow (QIIME 2) with Pathogen Curation & aDNA Authentication for 16S rRNA reads from New York African Burial Ground soil samples.
End-to-end, reproducible tooling for 16S/amplicon microbiome analyses, plus utilities for human-associated/pathogen curation and ancient DNA (aDNA) authentication.
All scripts are cleaned, parameterized, and ready for third-party reuse.

### Highlights
- QIIME 2 wrappers for core metrics, ordinations, heatmaps, and bar plots
-  region-specific classifier training + taxonomy export utilities
- iTOL-ready targeted trees for top genera
- Human-associated/pathogen curation (Python, TSV-in/TSV-out)
- aDNA verification

### Purpose
This repository contains scripts, configuration files, and workflow documentation for the analysis of 16S rRNA amplicon sequencing data from the New York African Burial Ground (NYABG) burial soils and nearby urban controls. The pipeline is designed to be reproducible and mirrors the bioinformatic workflow used in our manuscript "Human-Associated Microbial Persistence and Ecological Succession in Burial Soils from the New York African Burial Ground".

### Overview
The analysis was conducted using QIIME2 (2022.2 and amplicon-2024.10), with additional statistical and visualization steps in R. Key stages include:
Demultiplexing and quality control
Denoising and chimera removal (DADA2)
Taxonomic assignment (SILVA 138)
Diversity analyses (alpha and beta diversity)
Differential abundance testing (ANCOM, ALDEx2)
Visualization and statistical validation (phyloseq in R)
Phylogenetic Reconstruciton (iTOL)

### Repository Structure
.
├── scripts/                               
│   ├── core_metrics_diversity.sh
│   ├── pcoa_feature_biplot.sh
│   ├── longitudinal_diversity.sh
│   ├── taxa_barplots_basic.sh
│   ├── taxa_barplots_plus.sh
│   ├── abundance_heatmap.sh
│   ├── silva_v4_train_and_classify.sh
│   ├── species_taxonomy_export.sh
│   ├── prepare_itol_top_genera_tree.sh
│   ├── sample_classifier_rf.sh
│   ├── adna_damageprofiler_pipeline.sh
│   └── human_pathogen_curation.py
├── README.md                               #  Project documentation



### Requirements

QIIME2 (tested with 2022.2 and amplicon-2024.10)
R (≥4.5.0) with packages: phyloseq, ALDEx2, tidyverse
Conda or Docker for reproducible environments

### Usage
Clone the repository:
git clone https://github.com/your-username/NYABG-16S-QIIME2.git
cd NYABG-16S-QIIME2

### Quickstart
### Set up the environment:
#### QIIME 2 stack (example release; adjust if needed)
mamba create -n qiime2-2024.5 -c qiime2 -c conda-forge qiime2 biom-format
mamba activate qiime2-2024.5

#### aDNA tooling
mamba create -n adna -c bioconda -c conda-forge bwa samtools damageprofiler
mamba activate adna

#### Python tools for curation (optional)
python -m venv .venv && source .venv/bin/activate
pip install pandas

##### Citation
If you use this repository, please cite:
[Clinton, C.K. and Jackson, F.L.C.]. (2025). Human-Associated Microbial Persistence and Ecological Succession in Burial Soils from the New York African Burial Ground. ISME Communications. [DOI not yet available]
