#!/usr/bin/env bash
# QIIME 2 Taxa Bar Plots and Rank-Collapsed Tables
# By Carter Clinton, Ph.D.
# taxa_barplots.sh — Generate SILVA-based taxa bar plots and (optionally) collapsed tables
# Maps original file: Taxa_Bar_Plots_SILVA.txt
# Analysis Stage: Taxonomy & Composition
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
Generate QIIME 2 taxa bar plots from a feature table + taxonomy + metadata.
Optionally produce rank-collapsed tables (e.g., genus) and summaries, plus relative frequency tables.

Required:
  -t, --table-qza FILE        Feature table (.qza; FeatureTable[Frequency])
  -x, --taxonomy-qza FILE     Taxonomy (.qza; FeatureData[Taxonomy])
  -m, --metadata FILE         Sample metadata TSV (QIIME 2 format)
  -o, --output-dir DIR        Output directory

Optional:
  -l, --label STR             Basename for outputs (default: taxbars)
  -q, --qiime-cmd PATH        'qiime' executable (default: qiime)
  -j, --threads INT           Threads for any threaded ops (reserved; not used by barplot) (default: 0)
  --collapse-levels LIST      Comma-separated taxonomy levels to collapse (1=K..7=S) (default: 6)
  --make-relative             Also compute relative-frequency tables for collapsed outputs (default: off)
  --no-summaries              Skip 'feature-table summarize' for collapsed tables (default: summaries on)
  -h, --help                  Show this help

Outputs (in <OUTPUT_DIR>):
  barplots/<label>_taxa-bar-plots.qzv
  collapse/<label>_L<level>_table.qza                 (for each requested level)
  collapse/<label>_L<level>_table.qzv                 (unless --no-summaries)
  collapse/<label>_L<level>_relfreq.qza               (if --make-relative)
  qiime-info.txt

Examples:
  taxa_barplots_silva.sh \
    -t <TABLE_QZA> -x <TAXONOMY_QZA> -m <METADATA_TSV> \
    -o <OUTPUT_DIR> -l runA

  # Also generate genus + species collapsed tables, with relative frequency
  taxa_barplots_silva.sh \
    -t <TABLE_QZA> -x <TAXONOMY_QZA> -m <METADATA_TSV> \
    -o <OUTPUT_DIR> -l runB \
    --collapse-levels 6,7 --make-relative
USAGE
}

# ---------------------------
# Defaults
# ---------------------------
QIIME_CMD="${QIIME_CMD:-qiime}"
TABLE_QZA=""
TAXONOMY_QZA=""
METADATA=""
OUTPUT_DIR=""
LABEL="taxbars"
THREADS=0
COLLAPSE_LEVELS="6"
MAKE_REL=0
DO_SUMMARIES=1

# ---------------------------
# Parse CLI
# ---------------------------
OPTS_SHORT="t:x:m:o:l:q:j:h"
OPTS_LONG="table-qza:,taxonomy-qza:,metadata:,output-dir:,label:,qiime-cmd:,threads:,help,collapse-levels:,make-relative,no-summaries"
PARSED=$(getopt --options="$OPTS_SHORT" --longoptions="$OPTS_LONG" --name "$0" -- "$@") || { usage; exit 2; }
eval set -- "$PARSED"
while true; do
  case "$1" in
    -t|--table-qza)    TABLE_QZA="$2"; shift 2 ;;
    -x|--taxonomy-qza) TAXONOMY_QZA="$2"; shift 2 ;;
    -m|--metadata)     METADATA="$2"; shift 2 ;;
    -o|--output-dir)   OUTPUT_DIR="$2"; shift 2 ;;
    -l|--label)        LABEL="$2"; shift 2 ;;
    -q|--qiime-cmd)    QIIME_CMD="$2"; shift 2 ;;
    -j|--threads)      THREADS="$2"; shift 2 ;;
    --collapse-levels) COLLAPSE_LEVELS="$2"; shift 2 ;;
    --make-relative)   MAKE_REL=1; shift 1 ;;
    --no-summaries)    DO_SUMMARIES=0; shift 1 ;;
    -h|--help)         usage; exit 0 ;;
    --) shift; break ;;
    *) die "Unknown option: $1" ;;
  esac
done

# ---------------------------
# Validate
# ---------------------------
[ -n "$TABLE_QZA" ]    || { usage; die "Missing --table-qza"; }
[ -n "$TAXONOMY_QZA" ] || { usage; die "Missing --taxonomy-qza"; }
[ -n "$METADATA" ]     || { usage; die "Missing --metadata"; }
[ -n "$OUTPUT_DIR" ]   || { usage; die "Missing --output-dir"; }
[ -f "$TABLE_QZA" ]    || die "Feature table not found: $TABLE_QZA"
[ -f "$TAXONOMY_QZA" ] || die "Taxonomy artifact not found: $TAXONOMY_QZA"
[ -f "$METADATA" ]     || die "Metadata TSV not found: $METADATA"
mkdir -p "$OUTPUT_DIR"/{barplots,collapse}
command -v "$QIIME_CMD" >/dev/null 2>&1 || die "Cannot find '$QIIME_CMD' in PATH. Activate QIIME 2."

# ---------------------------
# Capture environment
# ---------------------------
"$QIIME_CMD" info > "$OUTPUT_DIR/qiime-info.txt" || true

# ---------------------------
# Taxa bar plots
# ---------------------------
BARS_QZV="$OUTPUT_DIR/barplots/${LABEL}_taxa-bar-plots.qzv"
log "Generating taxa bar plots → $BARS_QZV"
"$QIIME_CMD" taxa barplot \
  --i-table "$TABLE_QZA" \
  --i-taxonomy "$TAXONOMY_QZA" \
  --m-metadata-file "$METADATA" \
  --o-visualization "$BARS_QZV"

# ---------------------------
# Optional: collapse tables by rank
# ---------------------------
IFS=',' read -r -a LVLS <<< "$COLLAPSE_LEVELS"
for L in "${LVLS[@]}"; do
  [[ "$L" =~ ^[1-7]$ ]] || { log "Skipping invalid collapse level '$L' (must be 1..7)"; continue; }
  OUT_QZA="$OUTPUT_DIR/collapse/${LABEL}_L${L}_table.qza"
  OUT_QZV="$OUTPUT_DIR/collapse/${LABEL}_L${L}_table.qzv"
  log "Collapsing table to taxonomy level $L → $OUT_QZA"
  "$QIIME_CMD" taxa collapse \
    --i-table "$TABLE_QZA" \
    --i-taxonomy "$TAXONOMY_QZA" \
    --p-level "$L" \
    --o-collapsed-table "$OUT_QZA"

  if [ "$DO_SUMMARIES" -eq 1 ]; then
    log "Summarizing collapsed table (L$L) → $OUT_QZV"
    "$QIIME_CMD" feature-table summarize \
      --i-table "$OUT_QZA" \
      --o-visualization "$OUT_QZV" \
      ${METADATA:+--m-sample-metadata-file "$METADATA"}
  fi

  if [ "$MAKE_REL" -eq 1 ]; then
    REL_QZA="$OUTPUT_DIR/collapse/${LABEL}_L${L}_relfreq.qza"
    log "Computing relative frequencies for collapsed table (L$L) → $REL_QZA"
    "$QIIME_CMD" feature-table relative-frequency \
      --i-table "$OUT_QZA" \
      --o-relative-frequency-table "$REL_QZA"
  fi
done

log "Done. Outputs in $OUTPUT_DIR:"
log "  barplots/${LABEL}_taxa-bar-plots.qzv"
for L in "${LVLS[@]}"; do
  [[ "$L" =~ ^[1-7]$ ]] || continue
  log "  collapse/${LABEL}_L${L}_table.qza"
  [ "$DO_SUMMARIES" -eq 1 ] && log "  collapse/${LABEL}_L${L}_table.qzv"
  [ "$MAKE_REL" -eq 1 ]     && log "  collapse/${LABEL}_L${L}_relfreq.qza"
done
