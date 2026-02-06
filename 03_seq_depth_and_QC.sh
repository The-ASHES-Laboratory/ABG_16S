#!/usr/bin/env bash
# QIIME 2 Sequencing Depth Assessment and Quality Control Summary
# By Carter Clinton, Ph.D.
# seq_depth_and_QC.sh — Summarize feature table, rep seqs, and denoising stats after DADA2
# Analysis Stage: QC & Depth Assessment
# Language: Bash
#
# Description:
#   Generates interactive .qzv visualizations for the feature table (sample/feature
#   counts, per-sample depth distribution) and representative sequences. Optionally
#   visualizes DADA2 denoising statistics. Use the outputs to choose a rarefaction
#   depth and verify that denoising retained sufficient reads.

set -euo pipefail

log() { printf "[%s] %s\n" "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

usage() {
  cat <<'USAGE'
Summarize DADA2 outputs: feature table depth, representative sequences, and denoising statistics.
Produces .qzv visualizations for interactive review in QIIME 2 View (https://view.qiime2.org).

Required:
  -t, --table-qza FILE          Feature table (.qza; FeatureTable[Frequency])
  -s, --rep-seqs-qza FILE       Representative sequences (.qza; FeatureData[Sequence])
  -o, --output-dir DIR           Output directory

Optional:
  -d, --denoising-stats FILE     Denoising stats (.qza; SampleData[DADA2Stats]) for loss summary
  -m, --metadata FILE            Sample metadata TSV (adds per-group summaries to table viz)
  -l, --label STR                Basename for outputs (default: qc)
  -q, --qiime-cmd PATH           qiime executable (default: qiime)
  -h, --help                     Show this help

Outputs (written to <OUTPUT_DIR>):
  <label>_table_summary.qzv          Per-sample/feature depth summary
  <label>_rep-seqs_summary.qzv       Tabulated representative sequences
  <label>_denoising-stats.qzv        Denoising statistics (if --denoising-stats provided)
  qiime-info.txt                     QIIME 2 environment capture

Examples:
  # Basic QC after DADA2
  03_seq_depth_and_QC.sh \
    -t dada2_table.qza \
    -s dada2_rep-seqs.qza \
    -o qc_output/ \
    -d dada2_denoising-stats.qza

  # With metadata for per-group summaries
  03_seq_depth_and_QC.sh \
    -t dada2_table.qza \
    -s dada2_rep-seqs.qza \
    -o qc_output/ \
    -m sample-metadata.tsv \
    -l burial_qc
USAGE
}

# ---------------------------
# Defaults
# ---------------------------
QIIME_CMD="${QIIME_CMD:-qiime}"
TABLE_QZA=""
REP_SEQS_QZA=""
DENOISING_STATS=""
METADATA=""
OUTPUT_DIR=""
LABEL="qc"

# ---------------------------
# Parse CLI
# ---------------------------
OPTS_SHORT="t:s:o:d:m:l:q:h"
OPTS_LONG="table-qza:,rep-seqs-qza:,output-dir:,denoising-stats:,metadata:,label:,qiime-cmd:,help"
PARSED=$(getopt --options="$OPTS_SHORT" --longoptions="$OPTS_LONG" --name "$0" -- "$@") || { usage; exit 2; }
eval set -- "$PARSED"
while true; do
  case "$1" in
    -t|--table-qza)        TABLE_QZA="$2"; shift 2 ;;
    -s|--rep-seqs-qza)     REP_SEQS_QZA="$2"; shift 2 ;;
    -o|--output-dir)       OUTPUT_DIR="$2"; shift 2 ;;
    -d|--denoising-stats)  DENOISING_STATS="$2"; shift 2 ;;
    -m|--metadata)         METADATA="$2"; shift 2 ;;
    -l|--label)            LABEL="$2"; shift 2 ;;
    -q|--qiime-cmd)        QIIME_CMD="$2"; shift 2 ;;
    -h|--help)             usage; exit 0 ;;
    --) shift; break ;;
    *) die "Unknown option: $1" ;;
  esac
done

# ---------------------------
# Validate
# ---------------------------
[ -n "$TABLE_QZA" ]    || { usage; die "Missing --table-qza"; }
[ -n "$REP_SEQS_QZA" ] || { usage; die "Missing --rep-seqs-qza"; }
[ -n "$OUTPUT_DIR" ]   || { usage; die "Missing --output-dir"; }
[ -f "$TABLE_QZA" ]    || die "Feature table not found: $TABLE_QZA"
[ -f "$REP_SEQS_QZA" ] || die "Rep seqs not found: $REP_SEQS_QZA"
[ -n "$DENOISING_STATS" ] && [ ! -f "$DENOISING_STATS" ] && die "Denoising stats not found: $DENOISING_STATS"
[ -n "$METADATA" ] && [ ! -f "$METADATA" ] && die "Metadata not found: $METADATA"
mkdir -p "$OUTPUT_DIR"
command -v "$QIIME_CMD" >/dev/null 2>&1 || die "Cannot find '$QIIME_CMD' in PATH. Activate your QIIME 2 env."

# ---------------------------
# Capture environment
# ---------------------------
"$QIIME_CMD" info > "$OUTPUT_DIR/qiime-info.txt" || true

# ---------------------------
# Summarize feature table
# ---------------------------
TABLE_QZV="$OUTPUT_DIR/${LABEL}_table_summary.qzv"
log "Summarizing feature table → $TABLE_QZV"
TABLE_CMD=( "$QIIME_CMD" feature-table summarize
  --i-table "$TABLE_QZA"
  --o-visualization "$TABLE_QZV" )
[ -n "$METADATA" ] && TABLE_CMD+=( --m-sample-metadata-file "$METADATA" )
"${TABLE_CMD[@]}"

# ---------------------------
# Tabulate representative sequences
# ---------------------------
SEQS_QZV="$OUTPUT_DIR/${LABEL}_rep-seqs_summary.qzv"
log "Tabulating representative sequences → $SEQS_QZV"
"$QIIME_CMD" feature-table tabulate-seqs \
  --i-data "$REP_SEQS_QZA" \
  --o-visualization "$SEQS_QZV"

# ---------------------------
# Visualize denoising stats (optional)
# ---------------------------
if [ -n "$DENOISING_STATS" ]; then
  STATS_QZV="$OUTPUT_DIR/${LABEL}_denoising-stats.qzv"
  log "Visualizing denoising statistics → $STATS_QZV"
  "$QIIME_CMD" metadata tabulate \
    --m-input-file "$DENOISING_STATS" \
    --o-visualization "$STATS_QZV"
fi

# ---------------------------
# Summary
# ---------------------------
log "Done. QC outputs in $OUTPUT_DIR:"
log "  ${LABEL}_table_summary.qzv      — Open in https://view.qiime2.org to check per-sample depth"
log "  ${LABEL}_rep-seqs_summary.qzv   — Review ASV sequences and lengths"
[ -n "$DENOISING_STATS" ] && log "  ${LABEL}_denoising-stats.qzv  — Check read retention through DADA2 steps"
log ""
log "Next step: Review the table summary to choose a rarefaction depth for 04_rarefaction.sh"
