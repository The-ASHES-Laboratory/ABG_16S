#!/usr/bin/env bash
# QIIME 2 ANCOM Differential Abundance Analysis
# By Carter Clinton, Ph.D.
# ancom.sh — Run QIIME 2 ANCOM at one or more taxonomy ranks with optional species-only variant
# Maps original file: ANCOM_SILVA.txt
# Analysis Stage: Differential Abundance
# Language: Bash

set -euo pipefail

log() { printf "[%s] %s\n" "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

usage() {
  cat <<'USAGE'
Run QIIME 2 ANCOM (composition ancom) to test differential abundance across metadata groups.
Optionally collapse the feature table to one or more taxonomy ranks and/or create a species-only variant.

Required:
  -t, --table-qza FILE           Input feature table (.qza; FeatureTable[Frequency])
  -x, --taxonomy-qza FILE        Taxonomy (.qza; FeatureData[Taxonomy]) for collapsing/filtering
  -m, --metadata FILE            Sample metadata TSV (QIIME 2 format)
  -c, --columns LIST             Comma-separated metadata columns to test (e.g., Group,Sex)
  -o, --output-dir DIR           Output directory

Optional:
  -l, --label STR                Basename/label for outputs (default: ancom)
  --levels LIST                  Comma-separated taxonomy levels to collapse (1=Kingdom..7=Species). If unset, runs on ASV table.
  --species-only                 Also run a species-only variant (filters to taxonomy strings containing 's__')
  --pseudocount INT              Pseudocount added before ANCOM (default: 1)
  --min-samples INT              (optional) Remove features present in <N samples before ANCOM (default: off)
  -q, --qiime-cmd PATH           'qiime' executable (default: qiime)
  -h, --help                     Show this help

Outputs (under <OUTPUT_DIR>/ancom/):
  base/
    <label>_ASV_composition.qza                 (if levels unset)
    <label>_ASV_<COLUMN>_ancom.qzv
  L<r>/
    <label>_L<r>_collapsed.qza
    <label>_L<r>_composition.qza
    <label>_L<r>_<COLUMN>_ancom.qzv             (per tested column)
  species_only/ (if --species-only)
    <label>_speciesonly_table.qza
    <label>_speciesonly_composition.qza
    <label>_speciesonly_<COLUMN>_ancom.qzv
  qiime-info.txt

Examples:
  # Genus-level and phylum-level ANCOM for two metadata columns, plus species-only variant
  ancom_silva.sh \
    -t <TABLE_QZA> -x <TAXONOMY_QZA> -m <METADATA_TSV> \
    -c Group,Sex -o <OUTPUT_DIR> -l runA \
    --levels 2,6 --species-only

  # Run ANCOM on the uncollapsed ASV table only
  ancom_silva.sh -t <TABLE_QZA> -x <TAXONOMY_QZA> -m <METADATA_TSV> \
    -c Group -o <OUTPUT_DIR> -l asv_only
USAGE
}

# Defaults
QIIME_CMD="${QIIME_CMD:-qiime}"
TABLE_QZA=""
TAXONOMY_QZA=""
METADATA=""
COLUMNS=""
OUTPUT_DIR=""
LABEL="ancom"
LEVELS=""
SPECIES_ONLY=0
PSEUDOCOUNT=1
MIN_SAMPLES=""

# Parse CLI
OPTS_SHORT="t:x:m:c:o:l:q:h"
OPTS_LONG="table-qza:,taxonomy-qza:,metadata:,columns:,output-dir:,label:,qiime-cmd:,help,levels:,species-only,pseudocount:,min-samples:"
PARSED=$(getopt --options="$OPTS_SHORT" --longoptions="$OPTS_LONG" --name "$0" -- "$@") || { usage; exit 2; }
eval set -- "$PARSED"
while true; do
  case "$1" in
    -t|--table-qza)    TABLE_QZA="$2"; shift 2 ;;
    -x|--taxonomy-qza) TAXONOMY_QZA="$2"; shift 2 ;;
    -m|--metadata)     METADATA="$2"; shift 2 ;;
    -c|--columns)      COLUMNS="$2"; shift 2 ;;
    -o|--output-dir)   OUTPUT_DIR="$2"; shift 2 ;;
    -l|--label)        LABEL="$2"; shift 2 ;;
    -q|--qiime-cmd)    QIIME_CMD="$2"; shift 2 ;;
    --levels)          LEVELS="$2"; shift 2 ;;
    --species-only)    SPECIES_ONLY=1; shift 1 ;;
    --pseudocount)     PSEUDOCOUNT="$2"; shift 2 ;;
    --min-samples)     MIN_SAMPLES="$2"; shift 2 ;;
    -h|--help)         usage; exit 0 ;;
    --) shift; break ;;
    *) die "Unknown option: $1" ;;
  esac
done

# Validate
[ -n "$TABLE_QZA" ]  || { usage; die "Missing --table-qza"; }
[ -n "$TAXONOMY_QZA" ] || { usage; die "Missing --taxonomy-qza"; }
[ -n "$METADATA" ]   || { usage; die "Missing --metadata"; }
[ -n "$COLUMNS" ]    || { usage; die "Missing --columns"; }
[ -n "$OUTPUT_DIR" ] || { usage; die "Missing --output-dir"; }
[ -f "$TABLE_QZA" ]  || die "Feature table not found: $TABLE_QZA"
[ -f "$TAXONOMY_QZA" ] || die "Taxonomy artifact not found: $TAXONOMY_QZA"
[ -f "$METADATA" ]   || die "Metadata TSV not found: $METADATA"
[[ "$PSEUDOCOUNT" =~ ^[0-9]+$ ]] || die "--pseudocount must be an integer"
[ -n "$MIN_SAMPLES" ] && [[ ! "$MIN_SAMPLES" =~ ^[0-9]+$ ]] && die "--min-samples must be an integer if provided"
mkdir -p "$OUTPUT_DIR/ancom"/{base} || true
command -v "$QIIME_CMD" >/dev/null 2>&1 || die "Cannot find '$QIIME_CMD' in PATH. Activate QIIME 2."

# Capture environment
"$QIIME_CMD" info > "$OUTPUT_DIR/qiime-info.txt" || true

# Helper: optionally prevalence-filter features
prevalence_filter() {
  local in_table="$1" out_table="$2" min_samples="$3"
  if [ -n "$min_samples" ]; then
    log "Prevalence-filter: keeping features present in >= $min_samples samples"
    "$QIIME_CMD" feature-table filter-features \
      --i-table "$in_table" \
      --p-min-samples "$min_samples" \
      --o-filtered-table "$out_table"
  else
    cp "$in_table" "$out_table"
  fi
}

# Helper: add pseudocount → composition table
to_composition() {
  local in_table="$1" out_comp="$2" pc="$3"
  if [ "$pc" -ne 1 ]; then
    log "NOTE: QIIME 2 add-pseudocount uses 1 by default; custom pseudocount=$pc requires prior scaling."
    log "      Proceeding with standard add-pseudocount (1)."
  fi
  "$QIIME_CMD" composition add-pseudocount \
    --i-table "$in_table" \
    --o-composition-table "$out_comp"
}

# Run ANCOM for a given composition table and one metadata column
run_ancom() {
  local comp_qza="$1" column="$2" out_qzv="$3"
  "$QIIME_CMD" composition ancom \
    --i-table "$comp_qza" \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$column" \
    --o-visualization "$out_qzv"
}

# Parse columns list
IFS=',' read -r -a COLS <<< "$COLUMNS"

# 0) ASV-level (uncollapsed) ANCOM (if no levels provided OR user still wants ASV)
if [ -z "$LEVELS" ]; then
  log "Preparing ASV-level ANCOM"
  PRE_ASV="$OUTPUT_DIR/ancom/base/${LABEL}_ASV_prefilter.qza"
  prevalence_filter "$TABLE_QZA" "$PRE_ASV" "$MIN_SAMPLES"
  COMP_ASV="$OUTPUT_DIR/ancom/base/${LABEL}_ASV_composition.qza"
  to_composition "$PRE_ASV" "$COMP_ASV" "$PSEUDOCOUNT"
  for col in "${COLS[@]}"; do
    COL_TRIM="$(echo "$col" | sed 's/^ *//;s/ *$//')"
    SAFE_COL="$(echo "$COL_TRIM" | tr ' ' '_' )"
    OUT_QZV="$OUTPUT_DIR/ancom/base/${LABEL}_ASV_${SAFE_COL}_ancom.qzv"
    log "ANCOM (ASV) for column '$COL_TRIM' → $OUT_QZV"
    run_ancom "$COMP_ASV" "$COL_TRIM" "$OUT_QZV"
  done
fi

# 1) Collapsed ranks, if provided
if [ -n "$LEVELS" ]; then
  mkdir -p "$OUTPUT_DIR/ancom"
  IFS=',' read -r -a LVLS <<< "$LEVELS"
  for L in "${LVLS[@]}"; do
    [[ "$L" =~ ^[1-7]$ ]] || { log "Skipping invalid level '$L' (1..7)"; continue; }
    DIR_L="$OUTPUT_DIR/ancom/L${L}"
    mkdir -p "$DIR_L"
    COLLAPSED="$DIR_L/${LABEL}_L${L}_collapsed.qza"
    PRE_L="$DIR_L/${LABEL}_L${L}_prefilter.qza"
    COMP_L="$DIR_L/${LABEL}_L${L}_composition.qza"

    log "Collapsing to taxonomy level $L → $COLLAPSED"
    "$QIIME_CMD" taxa collapse \
      --i-table "$TABLE_QZA" \
      --i-taxonomy "$TAXONOMY_QZA" \
      --p-level "$L" \
      --o-collapsed-table "$COLLAPSED"

    prevalence_filter "$COLLAPSED" "$PRE_L" "$MIN_SAMPLES"
    to_composition "$PRE_L" "$COMP_L" "$PSEUDOCOUNT"

    for col in "${COLS[@]}"; do
      COL_TRIM="$(echo "$col" | sed 's/^ *//;s/ *$//')"
      SAFE_COL="$(echo "$COL_TRIM" | tr ' ' '_' )"
      OUT_QZV="$DIR_L/${LABEL}_L${L}_${SAFE_COL}_ancom.qzv"
      log "ANCOM (L$L) for column '$COL_TRIM' → $OUT_QZV"
      run_ancom "$COMP_L" "$COL_TRIM" "$OUT_QZV"
    done
  done
fi

# 2) Species-only variant (taxonomy strings containing "s__")
if [ "$SPECIES_ONLY" -eq 1 ]; then
  DIR_S="$OUTPUT_DIR/ancom/species_only"
  mkdir -p "$DIR_S"
  TABLE_S="$DIR_S/${LABEL}_speciesonly_table.qza"
  PRE_S="$DIR_S/${LABEL}_speciesonly_prefilter.qza"
  COMP_S="$DIR_S/${LABEL}_speciesonly_composition.qza"

  log "Filtering to species-only features (taxonomy contains 's__') → $TABLE_S"
  "$QIIME_CMD" taxa filter-table \
    --i-table "$TABLE_QZA" \
    --i-taxonomy "$TAXONOMY_QZA" \
    --p-include "s__" \
    --o-filtered-table "$TABLE_S"

  prevalence_filter "$TABLE_S" "$PRE_S" "$MIN_SAMPLES"
  to_composition "$PRE_S" "$COMP_S" "$PSEUDOCOUNT"

  for col in "${COLS[@]}"; do
    COL_TRIM="$(echo "$col" | sed 's/^ *//;s/ *$//')"
    SAFE_COL="$(echo "$COL_TRIM" | tr ' ' '_' )"
    OUT_QZV="$DIR_S/${LABEL}_speciesonly_${SAFE_COL}_ancom.qzv"
    log "ANCOM (species-only) for column '$COL_TRIM' → $OUT_QZV"
    run_ancom "$COMP_S" "$COL_TRIM" "$OUT_QZV"
  done
fi

log "Done. ANCOM outputs under: $OUTPUT_DIR/ancom"
