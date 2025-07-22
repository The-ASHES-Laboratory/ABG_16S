# 16S rRNA Analysis of New York African Burial Ground Soil Samples
Description

Scripts and workflows for 16S rRNA amplicon analysis of New York African Burial Ground soil samples using QIIME2. Includes demultiplexing, taxonomy assignment, diversity analyses, and differential abundance testing with ALDEx2 and ANCOM, ensuring reproducibility of the study.

Purpose

This repository contains scripts, configuration files, and workflow documentation for the analysis of 16S rRNA amplicon sequencing data from the New York African Burial Ground (NYABG) burial soils and nearby urban controls. The pipeline is designed to be reproducible and mirrors the bioinformatic workflow used in our manuscript "Human-Associated Microbial Persistence and Ecological Succession in Burial Soils from the New York African Burial Ground".

Overview

The analysis was conducted using QIIME2 (2022.2 and amplicon-2024.10), with additional statistical and visualization steps in R. Key stages include:
Demultiplexing and quality control
Denoising and chimera removal (DADA2)
Taxonomic assignment (SILVA 138)
Diversity analyses (alpha and beta diversity)
Differential abundance testing (ANCOM, ALDEx2)
Visualization and statistical validation (phyloseq in R)
Phylogenetic Reconstruciton (iTOL)

Repository Structure

├── data/             # Example input data formats or metadata templates
├── scripts/          # QIIME2 and R scripts for analysis
├── results/          # Example outputs (optional)
├── environment.yml   # Conda environment for QIIME2
├── README.md         # Project documentation


Requirements

QIIME2 (tested with 2022.2 and amplicon-2024.10)
R (≥4.5.0) with packages: phyloseq, ALDEx2, tidyverse
Conda or Docker for reproducible environments

Usage

Clone the repository:
git clone https://github.com/your-username/NYABG-16S-QIIME2.git
cd NYABG-16S-QIIME2


Set up the environment:
conda env create -f environment.yml
conda activate nyabg-qiime2
Run the pipeline using the scripts in the scripts/ directory. Each step is documented with input/output details.

Citation
If you use this repository, please cite:

[Clinton, C.K. and Jackson, F.L.C.]. (2025). Human-Associated Microbial Persistence and Ecological Succession in Burial Soils from the New York African Burial Ground. ISME Journal. [DOI not yet available]
