# Architecture

This document describes the analytical structure of the `ABG_16S` pipeline.

## Goal

Characterize the bacterial community composition of 17th- and 18th-century burial soils from the New York African Burial Ground (NYABG) and nearby urban control sites, and assess whether human-associated taxa persist in the burial soils above the background expected for urban soils.

## Data model

The pipeline operates on:

- **81 burial soil samples + 6 urban control samples** sequenced on Illumina MiSeq (V3-V4 region, paired-end).
- **Sample metadata** in QIIME 2 TSV format (sample ID, burial vs control, additional grouping variables).
- **SILVA 138 classifier** (either pre-trained or trained in-pipeline by `05_taxonomy.sh`).

All intermediate artifacts are QIIME 2 `.qza` (or `.qzv` for visualizations); the pipeline reads and writes these natively and does not require manual format conversion between steps.

## Stages

```
01_demultiplex.sh
      │  imports demultiplexed paired-end reads, generates demux summary
      ▼
02_denoise.sh
      │  DADA2 paired-end with --trunc-len-f / --trunc-len-r and consensus chimera removal
      ▼
03_seq_depth_and_QC.sh ─────────────┐
      │  table and rep-seq QC,      │   (depth selection
      │  feeds rarefaction depth    │    feeds 04 and 08)
      ▼                             │
04_rarefaction.sh ──────────────────┤
      │  observed features, Shannon,│
      │  Faith's PD alpha curves    │
      ▼                             │
05_taxonomy.sh ─────────────────────┤
      │  SILVA 138 Naive Bayes       │
      │  classification             │
      ▼                             │
06_taxa_barplots.sh                 │
      │  taxa bar plots,            │
      │  rank-collapsed tables      │
      ▼                             │
07_build_phylogeny.sh ──────────────┤
      │  MAFFT + FastTree +         │
      │  midpoint rooting           │
      ▼                             ▼
08_core_metrics_diversity.sh ◄──────┘
      │  alpha + beta metrics, PERMANOVA, ANOSIM, pairwise tests
      ▼
09_permdisp.sh
      │  homogeneity of dispersion testing
      ▼
10_ancom.sh + 11_aldex2.R
      │  parallel differential abundance methods
      ▼
12_generate_itol_tree.sh
      │  iTOL-ready trees with color-strip annotations
      ▼
13_human_pathogen_curation.py
      │  human-associated and pathogenic taxa curation
      ▼
14_adna_damage_profiler.sh
      │  DamageProfiler authentication of ancient DNA signal
      ▼
   (manuscript figures and tables)
```

## Reproducibility model

Every wrapper script:

1. Captures the active QIIME 2 environment (`qiime info`) into the output directory.
2. Accepts all parameters via CLI flags (`--help` lists them).
3. Writes a log file recording the parameter values used.
4. Produces QIIME 2 `.qza` / `.qzv` artifacts that can be inspected via [QIIME 2 View](https://view.qiime2.org/) without re-running the upstream steps.

This means a reviewer can re-run a single stage (for example, recompute differential abundance with a tighter prevalence filter) without re-running denoising, taxonomy, or phylogeny.

## Authentication

Burial soils can be contaminated by modern handling during excavation and storage. The pipeline addresses this with two complementary checks:

- **Damage profiling (`14_adna_damage_profiler.sh`):** Computes C-to-T deamination patterns at read termini for taxa identified as human-associated. Authentic ancient DNA shows characteristic 5'-C-to-T and 3'-G-to-A patterns; modern contamination does not.
- **Urban control comparison (PERMANOVA in `08_core_metrics_diversity.sh`):** Burial samples are tested against the 6 urban control samples to confirm that the human-associated signal is not a background feature of urban soils.

## Outputs

Each script writes outputs to a numbered directory matching its script number (`01_demux/`, `02_denoise/`, ...). These directories are not tracked in the repository; the pipeline is designed to be re-run on the user's hardware against the published raw data.
