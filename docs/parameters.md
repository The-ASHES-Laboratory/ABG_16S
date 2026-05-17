# Key parameters

This document collects the parameters that most influence pipeline output. Defaults reproduce the manuscript results (Clinton & Jackson, 2025, *ISME Communications* 5(1)).

## Denoising (`02_denoise.sh`)

| Parameter | Default | Rationale |
|---|---|---|
| `--trunc-len-f` | 240 | Forward-read DADA2 truncation; chosen to keep median Q-score above 30. |
| `--trunc-len-r` | 220 | Reverse-read DADA2 truncation; reverse reads degrade earlier than forward on V3-V4 MiSeq runs. |
| `--chimera-method` | consensus | DADA2 default; balances sensitivity and specificity for chimera detection. |
| `--p-n-threads` | 4 | Reasonable for a workstation; raise for cluster runs. |

Adjust truncation lengths based on the rep-seq quality plot from `03_seq_depth_and_QC.sh`. Re-running denoising is by far the most expensive step; verify your truncation choice before committing to a full pipeline run.

## Sequencing depth and rarefaction (`03_seq_depth_and_QC.sh`, `04_rarefaction.sh`)

| Parameter | Default | Rationale |
|---|---|---|
| Rarefaction `--max-depth` | 1,000 reads/sample | Captures the burial samples (median depth ~6,000 post-denoising); raise if your samples have higher depth. |
| Rarefaction iterations | 10 | QIIME 2 default; tighter than 1 (which is the visual default) without becoming computationally expensive. |

The rarefaction depth chosen for downstream analysis (`08_core_metrics_diversity.sh`) should be informed by the curves from `04_rarefaction.sh`. The manuscript used 1,000 as a depth that retains 100% of samples while approaching the alpha-diversity plateau.

## Taxonomy classification (`05_taxonomy.sh`)

| Parameter | Default | Rationale |
|---|---|---|
| Classifier | SILVA 138 Naive Bayes | Standard reference for 16S; manuscript-tracked version. |
| Region trim | V3-V4 (515F/806R primers can be enabled by `--train-region`) | Region-specific classification improves precision for short amplicons. |

Region-specific classifier training is optional but recommended when the amplicon region is shorter than full 16S; pre-trained classifiers for common regions are available from QIIME 2 data resources.

## Diversity analysis (`08_core_metrics_diversity.sh`)

| Parameter | Default | Rationale |
|---|---|---|
| `--sampling-depth` | 1,000 | Set by the rarefaction analysis. |
| PERMANOVA permutations | 999 | Standard convention; tighter than 99. |
| `--m-metadata-column` | `burial_vs_control` | The primary axis of comparison; pairwise tests are run for any additional metadata columns specified. |

## Differential abundance (`10_ancom.sh`, `11_aldex2.R`)

| Parameter | Default | Rationale |
|---|---|---|
| ANCOM `--m-metadata-column` | `burial_vs_control` | Primary comparison. |
| ALDEx2 effect-size threshold | |effect| > 1 | Standard ALDEx2 convention for "biologically meaningful". |
| Cross-method consensus | Both methods must flag | Reduces method-specific false positives. |

## Phylogeny (`07_build_phylogeny.sh`)

| Parameter | Default | Rationale |
|---|---|---|
| Aligner | MAFFT | QIIME 2 default; fast and accurate for 16S. |
| Tree method | FastTree | Fast approximate-maximum-likelihood; sufficient for UniFrac. Switch to IQ-TREE if a publication-grade tree is required (see `docs/extending.md`). |
| Rooting | Midpoint | Standard for unrooted trees from FastTree. |

## Authentication (`14_adna_damage_profiler.sh`)

| Parameter | Default | Rationale |
|---|---|---|
| `--reference` | per-taxon Genus-level reference | Authentic ancient signal is expected at the taxon level; per-sample alignment is unnecessary. |
| `--minreadlen` | 30 | DamageProfiler default; shorter reads have unreliable damage patterns. |

## Reproducibility flags common to every script

- `--seed` — RNG seed (default 42); guarantees deterministic outputs where stochasticity is involved.
- `--threads` — parallelism, where supported.
- `qiime info` capture into the output directory is automatic.
