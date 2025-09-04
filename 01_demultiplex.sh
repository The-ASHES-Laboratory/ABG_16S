##### 01_Demultiplex
##### By Carter Clinton



#!/usr/bin/env bash
# demultiplex_summarize.sh — Summarize demultiplexed reads (QIIME 2 demux summarize)
# Analysis Stage: Import
# Language: Bash

set -euo pipefail

# ---------------------------
# Logging helpers
# ---------------------------
log() { printf "[%s] %s\n" "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

# ---------------------------
# Usage
# ---------------------------
usage() {
  cat <<'USAGE'
Summarize demultiplexed reads in QIIME 2 (.qza) and optionally import from a manifest.

Required:
  -o, --output-dir             Output directory for artifacts/reports

One of (A) or (B):
  (A) Use an existing demultiplexed artifact:
      -i, --input-demux        Path to demultiplexed artifact (demux.qza)

  (B) Import from a manifest (optional convenience):
      -M, --manifest           Path to QIIME 2 manifest CSV/TSV
      -t, --type               QIIME 2 semantic type (default: SampleData[PairedEndSequencesWithQuality])
      -f, --format             QIIME 2 input format (default: PairedEndFastqManifestPhred33V2)

Optional:
  -l, --label                  Basename/label for outputs (default: demux)
  -q, --qiime-cmd              qiime executable (default: qiime in PATH)
  -h, --help                   Show this help

Environment:
  Ensure a QIIME 2 environment is active (e.g., conda activate qiime2-<version>).

Outputs:
  <output-dir>/<label>.qza           (if importing)
  <output-dir>/<label>-summary.qzv   (visualization)

Example:
  demultiplex_summarize.sh -i <DEMUX_QZA> -o <OUTPUT_DIR> -l runA
  demultiplex_summarize.sh -M <MANIFEST.csv> -o <OUTPUT_DIR> -l runB
USAGE
}

# ---------------------------
# Defaults
# ---------------------------
QIIME_CMD="${QIIME_CMD:-qiime}"
LABEL="demux"
Q2_TYPE="SampleData[PairedEndSequencesWithQuality]"
Q2_FORMAT="PairedEndFastqManifestPhred33V2"

INPUT_DEMUX=""
MANIFEST=""
OUTPUT_DIR=""

# ---------------------------
# Parse CLI
# ---------------------------
short="i:o:l:q:M:t:f:h"
long="input-demux:,output-dir:,label:,qiime-cmd:,manifest:,type:,format:,help"
PARSED=$(getopt --options="$short" --longoptions="$long" --name "$0" -- "$@") || { usage; exit 2; }
eval set -- "$PARSED"
while true; do
  case "$1" in
    -i|--input-demux) INPUT_DEMUX="$2"; shift 2 ;;
    -o|--output-dir)  OUTPUT_DIR="$2"; shift 2 ;;
    -l|--label)       LABEL="$2"; shift 2 ;;
    -q|--qiime-cmd)   QIIME_CMD="$2"; shift 2 ;;
    -M|--manifest)    MANIFEST="$2"; shift 2 ;;
    -t|--type)        Q2_TYPE="$2"; shift 2 ;;
    -f|--format)      Q2_FORMAT="$2"; shift 2 ;;
    -h|--help)        usage; exit 0 ;;
    --) shift; break ;;
    *) die "Unknown option: $1" ;;
  esac
done

# ---------------------------
# Validate inputs
# ---------------------------
[ -n "$OUTPUT_DIR" ] || { usage; die "Missing --output-dir"; }
mkdir -p "$OUTPUT_DIR"

command -v "$QIIME_CMD" >/dev/null 2>&1 || die "Cannot find '$QIIME_CMD' in PATH. Activate your QIIME 2 env first."

if [ -n "$INPUT_DEMUX" ] && [ -n "$MANIFEST" ]; then
  die "Provide either --input-demux OR --manifest, not both."
fi
if [ -z "$INPUT_DEMUX" ] && [ -z "$MANIFEST" ]; then
  die "Provide one of --input-demux or --manifest."
fi

# ---------------------------
# If manifest provided, import to demux.qza
# ---------------------------
DEMUX_QZA="$INPUT_DEMUX"
if [ -n "$MANIFEST" ]; then
  [ -f "$MANIFEST" ] || die "Manifest not found: $MANIFEST"
  DEMUX_QZA="$OUTPUT_DIR/${LABEL}.qza"
  log "Importing reads from manifest → $DEMUX_QZA"
  "$QIIME_CMD" tools import \
    --type "$Q2_TYPE" \
    --input-path "$MANIFEST" \
    --output-path "$DEMUX_QZA" \
    --input-format "$Q2_FORMAT"
else
  [ -f "$DEMUX_QZA" ] || die "Input demux artifact not found: $DEMUX_QZA"
fi

# ---------------------------
# Summarize demultiplexed reads
# ---------------------------
SUMMARY_QZV="$OUTPUT_DIR/${LABEL}-summary.qzv"
log "Summarizing demultiplexed reads → $SUMMARY_QZV"
"$QIIME_CMD" demux summarize \
  --i-data "$DEMUX_QZA" \
  --o-visualization "$SUMMARY_QZV"

log "Done."
log "Artifacts:"
log "  DEMUX_QZA: $DEMUX_QZA"
log "  SUMMARY_QZV: $SUMMARY_QZV"


