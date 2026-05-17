# Changelog

All notable changes to this project will be documented in this file. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to semantic versioning.

## [Unreleased]

### Planned

- Companion notebooks rendering the manuscript figures from the QIIME 2 outputs.
- Snakemake / Nextflow wrappers for end-to-end pipeline execution.
- Support for QIIME 2 amplicon-2025 releases as they stabilize.

## [1.1] - 2026-02

### Changed

- Repository reorganized to a folder-based layout: pipeline scripts moved into `scripts/`.
- License renamed `LICENSE` → `LICENSE.md`.
- README updated to use the new script paths.

### Added

- `CHANGELOG.md`, `CITATION.cff`, `CONTRIBUTING.md`.
- `docs/architecture.md`, `docs/extending.md`, `docs/parameters.md`.

## [1.0] - 2025-09

### Added

- **QIIME 2 wrappers (scripts 01-10, 12, 14):**
  - Read import and demultiplex summary (`01_demultiplex.sh`)
  - DADA2 paired-end denoising with parameterized truncation (`02_denoise.sh`)
  - Feature table and rep-seq QC summaries (`03_seq_depth_and_QC.sh`)
  - Alpha-rarefaction curves over observed features, Shannon, Faith's PD (`04_rarefaction.sh`)
  - SILVA 138 Naive Bayes classification with optional region-specific training (`05_taxonomy.sh`)
  - Taxa bar plots and rank-collapsed frequency tables (`06_taxa_barplots.sh`)
  - Phylogenetic tree via MAFFT + FastTree + midpoint rooting (`07_build_phylogeny.sh`)
  - Core diversity metrics, PERMANOVA, ANOSIM, pairwise tests (`08_core_metrics_diversity.sh`)
  - PERMDISP homogeneity of dispersion testing (`09_permdisp.sh`)
  - ANCOM compositional differential abundance (`10_ancom.sh`)
  - iTOL-ready phylogenetic trees with color-strip annotations (`12_generate_itol_tree.sh`)
  - Ancient DNA damage profiling via DamageProfiler (`14_adna_damage_profiler.sh`)
- **R differential abundance (script 11):**
  - ALDEx2 differential abundance with CLR effect sizes (`11_aldex2.R`)
- **Python curation (script 13):**
  - Human-associated and pathogenic taxa curation (`13_human_pathogen_curation.py`)
- Reproducibility scaffolding: every wrapper captures the active QIIME 2 environment (`qiime info`) in its output directory.
