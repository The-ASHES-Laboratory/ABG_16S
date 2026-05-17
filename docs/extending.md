# Extending the framework

This document gives entry points for the most common extensions.

## Adapting to a different 16S variable region

The pipeline defaults to V3-V4 (the region used in the manuscript). To adapt to V4-only (515F/806R), V4-V5, or full-length 16S:

1. Adjust DADA2 truncation lengths in `02_denoise.sh` (`--trunc-len-f`, `--trunc-len-r`) based on the rep-seq quality of your run.
2. If using a different region, retrain the SILVA classifier in `05_taxonomy.sh`:
   ```bash
   bash scripts/05_taxonomy.sh --train-region 515f-806r \
     --silva-seqs silva-138-seqs.qza --silva-tax silva-138-tax.qza
   ```
3. Re-run downstream stages (06 onward); the rep-seq QC in `03_seq_depth_and_QC.sh` will surface any issues with the new truncation choice before you commit to a full pipeline run.

## Adding a new differential abundance method

The pipeline currently runs ANCOM (`10_ancom.sh`) and ALDEx2 (`11_aldex2.R`) in parallel. To add MaAsLin2, Corncob, or LinDA:

1. Add the new method as a numbered script following the same convention (e.g., `11a_maaslin2.R`).
2. Write outputs to a TSV with the same column schema as the existing methods: `feature_id`, `effect_size`, `pvalue`, `qvalue`.
3. Document any new R/conda dependencies in the script header.
4. Surface the new method in the manuscript Supplementary Methods.

## Adding an additional reference panel for classification

The current pipeline uses SILVA 138. To add or compare against Greengenes2 or RDP:

1. Add a `--classifier` flag to `05_taxonomy.sh` and gate the SILVA-specific training on its value.
2. Run `05_taxonomy.sh` once per panel, writing to per-panel output directories.
3. The downstream pipeline is panel-agnostic; only the taxonomy artifact path needs to change.

## Replacing FastTree with IQ-TREE

`07_build_phylogeny.sh` uses MAFFT + FastTree for speed. For a more rigorous phylogeny:

1. Replace the FastTree call with an IQ-TREE invocation (`iqtree -s aligned.fasta -m TEST -bb 1000`).
2. Keep the midpoint-rooting step downstream — IQ-TREE outputs are unrooted by default.
3. Note in your manuscript that the tree topology underlying UniFrac was inferred with IQ-TREE; reviewers will want this documented.

## Reproducibility expectations

Every extension should preserve the existing reproducibility guarantees:

- QIIME 2 environment captured in each output directory.
- All parameters exposed via CLI flags.
- Random seeds (DADA2, FastTree, classifier training) settable and logged.
- Output artifacts produced as QIIME 2 `.qza` / `.qzv` where applicable so they remain inspectable without re-running upstream stages.

If a proposed extension would break reproducibility of the v1.0 ISME Communications figures, surface this in the pull request description so reviewers can decide whether to gate it behind a version bump.
