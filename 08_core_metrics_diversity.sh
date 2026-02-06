#!/usr/bin/env bash
# QIIME 2 Core Diversity Metrics and Group Significance Tests
# By Carter Clinton, Ph.D.
# core_metrics_diversity.sh — Run QIIME 2 core-metrics-phylogenetic + group significance tests
# Maps original file: Core_Metrics_SILVA.txt
# Analysis Stage: Phylogeny & Diversity
# Language: Bash

set -euo pipefail

# ---------------------------
# Logging / errors
# ---------------------------
log() { printf "[%s] %s\n" "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

# ---------------------------
# Help
# ---------------------------
usage() {
  cat <<'USAGE'
Run QIIME 2 `core-metrics-phylogenetic` with a rooted tree and perform alpha/beta group significance tests.

Required:
  -t, --table-qza FILE           Feature table (.qza; FeatureTable[Frequency])
  -p, --phylogeny-qza FILE       Rooted phylogeny (.qza; Phylogeny[Rooted])
  -m, --metadata FILE            Sample metadata TSV (QIIME 2 format)
  -o, --output-dir DIR           Output directory
  -d, --depth INT                Sampling depth for rarefaction (e.g., 1000)

Optional:
  -l, --label STR                Basename/label for outputs (default: div)
  -q, --qiime-cmd PATH           'qiime' executable (default: qiime)
  -j, --threads INT              Threads for underlying steps when supported (default: 0 → auto)
  --beta-columns LIST            Comma-separated metadata columns to test for beta-group-significance (default: none)
  --alpha-columns LIST           Comma-separated metadata columns just for documentation (alpha tests ignore column and use full metadata)
  --permutations INT             Number of permutations for significance tests (default: 999)
  --pairwise                     Enable pairwise post-hoc tests (default: off)
  --beta-method {permanova,anosim,permdisp}  Method for beta tests (default: permanova)
  --also-permdisp                Additionally run PERMDISP alongside the chosen beta method (default: off)
  --skip-alpha-tests             Skip alpha-group-significance (default: run)
  --skip-beta-tests              Skip beta-group-significance (default: run)
  --no-emperor                   Do not attach metadata to emperor plots in core-metrics (default: attach)
  -h, --help                     Show help

Outputs (under <OUTPUT_DIR>):
  core_metrics/                  (core-metrics-phylogenetic artifacts: distances, PCoAs, emperors, alpha vectors, rarefied table)
  alpha_significance/            (QZVs per alpha metric)
  beta_significance/<col>/       (QZVs per distance × method for each metadata column)
  qiime-info.txt                 (environment capture)

Examples:
  core_metrics_diversity.sh \
    -t <TABLE_QZA> -p <ROOTED_TREE_QZA> -m <METADATA_TSV> \
    -o <OUTPUT_DIR> -d 1000 -l runA \
    --beta-columns Group,Sex --permutations 999 --pairwise --also-permdisp
USAGE
}

# ---------------------------
# Defaults
# ---------------------------
QIIME_CMD="${QIIME_CMD:-qiime}"
TABLE_QZA=""
TREE_QZA=""
METADATA=""
OUTPUT_DIR=""
LABEL="div"
DEPTH=""
THREADS=0
BETA_COLUMNS=""
ALPHA_COLUMNS=""
PERMS=999
PAIRWISE=0
BETA_METHOD="permanova"
ALSO_PERMDISP=0
SKIP_ALPHA=0
SKIP_BETA=0
ATTACH_META=1

# ---------------------------
# Parse CLI
# ---------------------------
OPTS_SHORT="t:p:m:o:d:l:q:j:h"
OPTS_LONG="table-qza:,phylogeny-qza:,metadata:,output-dir:,depth:,label:,qiime-cmd:,threads:,help,beta-columns:,alpha-columns:,permutations:,pairwise,beta-method:,also-permdisp,skip-alpha-tests,skip-beta-tests,no-emperor"
PARSED=$(getopt --options="$OPTS_SHORT" --longoptions="$OPTS_LONG" --name "$0" -- "$@") || { usage; exit 2; }
eval set -- "$PARSED"
while true; do
  case "$1" in
    -t|--table-qza)     TABLE_QZA="$2"; shift 2 ;;
    -p|--phylogeny-qza) TREE_QZA="$2"; shift 2 ;;
    -m|--metadata)      METADATA="$2"; shift 2 ;;
    -o|--output-dir)    OUTPUT_DIR="$2"; shift 2 ;;
    -d|--depth)         DEPTH="$2"; shift 2 ;;
    -l|--label)         LABEL="$2"; shift 2 ;;
    -q|--qiime-cmd)     QIIME_CMD="$2"; shift 2 ;;
    -j|--threads)       THREADS="$2"; shift 2 ;;
    --beta-columns)     BETA_COLUMNS="$2"; shift 2 ;;
    --alpha-columns)    ALPHA_COLUMNS="$2"; shift 2 ;;
    --permutations)     PERMS="$2"; shift 2 ;;
    --pairwise)         PAIRWISE=1; shift 1 ;;
    --beta-method)      BETA_METHOD="$2"; shift 2 ;;
    --also-permdisp)    ALSO_PERMDISP=1; shift 1 ;;
    --skip-alpha-tests) SKIP_ALPHA=1; shift 1 ;;
    --skip-beta-tests)  SKIP_BETA=1; shift 1 ;;
    --no-emperor)       ATTACH_META=0; shift 1 ;;
    -h|--help)          usage; exit 0 ;;
    --) shift; break ;;
    *) die "Unknown option: $1" ;;
  esac
done

# ---------------------------
# Validate
# ---------------------------
[ -n "$TABLE_QZA" ]  || { usage; die "Missing --table-qza"; }
[ -n "$TREE_QZA" ]   || { usage; die "Missing --phylogeny-qza"; }
[ -n "$METADATA" ]   || { usage; die "Missing --metadata"; }
[ -n "$OUTPUT_DIR" ] || { usage; die "Missing --output-dir"; }
[ -n "$DEPTH" ]      || { usage; die "Missing --depth"; }
[ -f "$TABLE_QZA" ]  || die "Feature table not found: $TABLE_QZA"
[ -f "$TREE_QZA" ]   || die "Rooted phylogeny not found: $TREE_QZA"
[ -f "$METADATA" ]   || die "Metadata TSV not found: $METADATA"
[[ "$DEPTH" =~ ^[0-9]+$ ]] || die "--depth must be an integer"
[[ "$PERMS" =~ ^[0-9]+$ ]] || die "--permutations must be an integer"
[[ "$BETA_METHOD" =~ ^(permanova|anosim|permdisp)$ ]] || die "--beta-method must be permanova, anosim, or permdisp"
mkdir -p "$OUTPUT_DIR"/{core_metrics,alpha_significance,beta_significance}
command -v "$QIIME_CMD" >/dev/null 2>&1 || die "Cannot find '$QIIME_CMD' in PATH. Activate QIIME 2."

# ---------------------------
# Capture environment
# ---------------------------
"$QIIME_CMD" info > "$OUTPUT_DIR/qiime-info.txt" || true

# ---------------------------
# core-metrics-phylogenetic
# ---------------------------
log "Running core-metrics-phylogenetic at depth=$DEPTH"
EMPEROR_ARGS=()
[ "$ATTACH_META" -eq 1 ] && EMPEROR_ARGS=( --m-metadata-file "$METADATA" )

CORE_OUT="$OUTPUT_DIR/core_metrics"
CMD=( "$QIIME_CMD" diversity core-metrics-phylogenetic
  --i-phylogeny "$TREE_QZA"
  --i-table "$TABLE_QZA"
  --p-sampling-depth "$DEPTH"
  --o-rarefied-table "$CORE_OUT/${LABEL}_rarefied_table.qza"
  --o-faith-pd-vector "$CORE_OUT/${LABEL}_faith_pd_vector.qza"
  --o-evenness-vector "$CORE_OUT/${LABEL}_pielou_evenness_vector.qza"
  --o-observed-features-vector "$CORE_OUT/${LABEL}_observed_features_vector.qza"
  --o-shannon-vector "$CORE_OUT/${LABEL}_shannon_vector.qza"
  --o-unweighted-unifrac-distance-matrix "$CORE_OUT/${LABEL}_unweighted_unifrac_dm.qza"
  --o-weighted-unifrac-distance-matrix "$CORE_OUT/${LABEL}_weighted_unifrac_dm.qza"
  --o-jaccard-distance-matrix "$CORE_OUT/${LABEL}_jaccard_dm.qza"
  --o-bray-curtis-distance-matrix "$CORE_OUT/${LABEL}_bray_curtis_dm.qza"
  --o-unweighted-unifrac-pcoa-results "$CORE_OUT/${LABEL}_unweighted_unifrac_pcoa.qza"
  --o-weighted-unifrac-pcoa-results "$CORE_OUT/${LABEL}_weighted_unifrac_pcoa.qza"
  --o-jaccard-pcoa-results "$CORE_OUT/${LABEL}_jaccard_pcoa.qza"
  --o-bray-curtis-pcoa-results "$CORE_OUT/${LABEL}_bray_curtis_pcoa.qza"
  --o-unweighted-unifrac-emperor "$CORE_OUT/${LABEL}_unweighted_unifrac_emperor.qzv"
  --o-weighted-unifrac-emperor "$CORE_OUT/${LABEL}_weighted_unifrac_emperor.qzv"
  --o-jaccard-emperor "$CORE_OUT/${LABEL}_jaccard_emperor.qzv"
  --o-bray-curtis-emperor "$CORE_OUT/${LABEL}_bray_curtis_emperor.qzv"
  "${EMPEROR_ARGS[@]}"
)
# threads param is not universally supported here; relying on plugin defaults
log "Running: ${CMD[*]}"
"${CMD[@]}"

# ---------------------------
# Alpha group significance
# ---------------------------
if [ "$SKIP_ALPHA" -eq 0 ]; then
  log "Running alpha-group-significance on four alpha metrics"
  declare -A ALPHA_MAP=(
    ["faith_pd"]="$CORE_OUT/${LABEL}_faith_pd_vector.qza"
    ["shannon"]="$CORE_OUT/${LABEL}_shannon_vector.qza"
    ["evenness"]="$CORE_OUT/${LABEL}_pielou_evenness_vector.qza"
    ["observed_features"]="$CORE_OUT/${LABEL}_observed_features_vector.qza"
  )
  for metric in "${!ALPHA_MAP[@]}"; do
    vec="${ALPHA_MAP[$metric]}"
    [ -f "$vec" ] || { log "WARNING: Missing alpha vector for $metric at $vec"; continue; }
    out_qzv="$OUTPUT_DIR/alpha_significance/${LABEL}_alpha_${metric}_group_significance.qzv"
    "$QIIME_CMD" diversity alpha-group-significance \
      --i-alpha-diversity "$vec" \
      --m-metadata-file "$METADATA" \
      --o-visualization "$out_qzv"
  done
else
  log "Skipping alpha-group-significance (per --skip-alpha-tests)."
fi

# ---------------------------
# Beta group significance
# ---------------------------
if [ "$SKIP_BETA" -eq 0 ] && [ -n "$BETA_COLUMNS" ]; then
  log "Running beta-group-significance (method=$BETA_METHOD; perms=$PERMS; pairwise=$PAIRWISE)"
  declare -A DM_MAP=(
    ["unweighted_unifrac"]="$CORE_OUT/${LABEL}_unweighted_unifrac_dm.qza"
    ["weighted_unifrac"]="$CORE_OUT/${LABEL}_weighted_unifrac_dm.qza"
    ["jaccard"]="$CORE_OUT/${LABEL}_jaccard_dm.qza"
    ["bray_curtis"]="$CORE_OUT/${LABEL}_bray_curtis_dm.qza"
  )
  IFS=',' read -r -a COLS <<< "$BETA_COLUMNS"
  for col in "${COLS[@]}"; do
    COL_TRIM="$(echo "$col" | sed 's/^ *//;s/ *$//')"
    SAFE_COL="$(echo "$COL_TRIM" | tr ' ' '_' )"
    OUT_DIR="$OUTPUT_DIR/beta_significance/${SAFE_COL}"
    mkdir -p "$OUT_DIR"

    for metric in "${!DM_MAP[@]}"; do
      dm="${DM_MAP[$metric]}"
      [ -f "$dm" ] || { log "WARNING: Missing distance matrix for $metric at $dm"; continue; }
      base="$OUT_DIR/${LABEL}_beta_${metric}_${SAFE_COL}_${BETA_METHOD}"
      args=( --i-distance-matrix "$dm"
             --m-metadata-file "$METADATA"
             --m-metadata-column "$COL_TRIM"
             --p-permutations "$PERMS"
             --o-visualization "${base}.qzv"
             --p-method "$BETA_METHOD" )
      [ "$PAIRWISE" -eq 1 ] && args+=( --p-pairwise )
      "$QIIME_CMD" diversity beta-group-significance "${args[@]}"

      if [ "$ALSO_PERMDISP" -eq 1 ] && [ "$BETA_METHOD" != "permdisp" ]; then
        base_pd="$OUT_DIR/${LABEL}_beta_${metric}_${SAFE_COL}_permdisp"
        "$QIIME_CMD" diversity beta-group-significance \
          --i-distance-matrix "$dm" \
          --m-metadata-file "$METADATA" \
          --m-metadata-column "$COL_TRIM" \
          --p-permutations "$PERMS" \
          --o-visualization "${base_pd}.qzv" \
          --p-method permdisp
      fi
    done
  done
elif [ -z "$BETA_COLUMNS" ] && [ "$SKIP_BETA" -eq 0 ]; then
  log "No --beta-columns provided; skipping beta-group-significance."
else
  log "Skipping beta-group-significance (per --skip-beta-tests)."
fi

log "Done. See outputs in $OUTPUT_DIR."
