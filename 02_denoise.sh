##### 02_denoise
##### By Carter Clinton



#!/usr/bin/env bash
# dada2_denoise_paired.sh — QIIME 2 DADA2 denoising for paired-end reads
# Analysis Stage: QC & Denoising
# Language: Bash
# Description: Run DADA2 (paired-end) with parameterized trimming/truncation and produce core artifacts and summaries.

set -euo pipefail

log() { printf "[%s] %s\n" "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

usage() {
  cat <<'USAGE'
Run QIIME 2 DADA2 denoise-paired on demultiplexed reads (.qza) with parameterized trimming/truncation.
Produces: feature table, representative sequences, denoising stats, and optional .qzv summaries.

Required:
  -i, --input-demux        Path to demultiplexed sequences artifact (.qza; SampleData[PairedEndSequencesWithQuality])
  -o, --output-dir         Output directory for artifacts/visualizations

Optional (I/O naming):
  -l, --label              Basename label for outputs (default: dada2)

Optional (DADA2 params; defaults reflect the manuscript):
  --trim-left-f N          Trim N bases from 5' of forward reads (default: 0)
  --trim-left-r N          Trim N bases from 5' of reverse reads (default: 0)
  --trunc-len-f N          Truncate forward reads at position N (default: 240)
  --trunc-len-r N          Truncate reverse reads at position N (default: 220)
  --trunc-q N              Truncate reads at first quality score <= N (unset by default)
  --max-ee-f X             Maximum expected errors for forward reads (unset by default)
  --max-ee-r X             Maximum expected errors for reverse reads (unset by default)
  --chimera-method M       Chimera removal method: consensus|pooled|none (default: consensus)
  --min-overlap N          Minimum overlap for merging (unset by default; DADA2 default applies)
  --n-reads-learn N        Subset of reads to learn error rates (default: 1000000 per DADA2)
  --p-hashed-seed S        Hashed seed for reproducible partitioning (optional; integer)
  --threads N              Threads for DADA2 (default: 0 → auto)

Optional (tooling):
  -q, --qiime-cmd PATH     qiime executable (default: qiime)
  --no-summaries           Skip generating table/rep-seqs .qzv summaries

Environment:
  Activate a QIIME 2 environment first (e.g., conda activate qiime2-<version>).

Outputs (written to <output-dir>):
  <label>_table.qza
  <label>_rep-seqs.qza
  <label>_denoising-stats.qza
  <label>_table.qzv                 (unless --no-summaries)
  <label>_rep-seqs.qzv              (unless --no-summaries)
  qiime-info.txt                    (QIIME env capture)

Examples:
  dada2_denoise_paired.sh -i <DEMUX_QZA> -o <OUTPUT_DIR> -l runA --trunc-len-f 240 --trunc-len-r 220
  dada2_denoise_paired.sh -i <DEMUX_QZA> -o <OUTPUT_DIR> --threads 8 --chimera-method consensus
USAGE
}

# Defaults
QIIME_CMD="${QIIME_CMD:-qiime}"
LABEL="dada2"
INPUT_DEMUX=""
OUTPUT_DIR=""
TRIM_LEFT_F=0
TRIM_LEFT_R=0
TRUNC_LEN_F=240
TRUNC_LEN_R=220
TRUNC_Q=""
MAX_EE_F=""
MAX_EE_R=""
CHIMERA_METHOD="consensus"
MIN_OVERLAP=""
N_READS_LEARN=""
HASHED_SEED=""
THREADS=0
DO_SUMMARIES=1

# Parse args
OPTS_SHORT="i:o:l:q:h"
OPTS_LONG="input-demux:,output-dir:,label:,qiime-cmd:,help,trim-left-f:,trim-left-r:,trunc-len-f:,trunc-len-r:,trunc-q:,max-ee-f:,max-ee-r:,chimera-method:,min-overlap:,n-reads-learn:,p-hashed-seed:,threads:,no-summaries"
PARSED=$(getopt --options="$OPTS_SHORT" --longoptions="$OPTS_LONG" --name "$0" -- "$@") || { usage; exit 2; }
eval set -- "$PARSED"
while true; do
  case "$1" in
    -i|--input-demux) INPUT_DEMUX="$2"; shift 2 ;;
    -o|--output-dir)  OUTPUT_DIR="$2"; shift 2 ;;
    -l|--label)       LABEL="$2"; shift 2 ;;
    -q|--qiime-cmd)   QIIME_CMD="$2"; shift 2 ;;
    --trim-left-f)    TRIM_LEFT_F="$2"; shift 2 ;;
    --trim-left-r)    TRIM_LEFT_R="$2"; shift 2 ;;
    --trunc-len-f)    TRUNC_LEN_F="$2"; shift 2 ;;
    --trunc-len-r)    TRUNC_LEN_R="$2"; shift 2 ;;
    --trunc-q)        TRUNC_Q="$2"; shift 2 ;;
    --max-ee-f)       MAX_EE_F="$2"; shift 2 ;;
    --max-ee-r)       MAX_EE_R="$2"; shift 2 ;;
    --chimera-method) CHIMERA_METHOD="$2"; shift 2 ;;
    --min-overlap)    MIN_OVERLAP="$2"; shift 2 ;;
    --n-reads-learn)  N_READS_LEARN="$2"; shift 2 ;;
    --p-hashed-seed)  HASHED_SEED="$2"; shift 2 ;;
    --threads)        THREADS="$2"; shift 2 ;;
    --no-summaries)   DO_SUMMARIES=0; shift 1 ;;
    -h|--help)        usage; exit 0 ;;
    --) shift; break ;;
    *) die "Unknown option: $1" ;;
  esac
done

# Validate
[ -n "$INPUT_DEMUX" ] || { usage; die "Missing --input-demux"; }
[ -n "$OUTPUT_DIR" ]  || { usage; die "Missing --output-dir"; }
[ -f "$INPUT_DEMUX" ] || die "Input artifact not found: $INPUT_DEMUX"
mkdir -p "$OUTPUT_DIR"
command -v "$QIIME_CMD" >/dev/null 2>&1 || die "Cannot find '$QIIME_CMD' in PATH. Activate QIIME 2."

# Capture environment
"$QIIME_CMD" info > "$OUTPUT_DIR/qiime-info.txt" || true

# Build DADA2 command
CMD=( "$QIIME_CMD" dada2 denoise-paired
  --i-demultiplexed-seqs "$INPUT_DEMUX"
  --p-trim-left-f "$TRIM_LEFT_F"
  --p-trim-left-r "$TRIM_LEFT_R"
  --p-trunc-len-f "$TRUNC_LEN_F"
  --p-trunc-len-r "$TRUNC_LEN_R"
  --p-n-threads "$THREADS"
  --p-chimera-method "$CHIMERA_METHOD"
  --o-table "$OUTPUT_DIR/${LABEL}_table.qza"
  --o-representative-sequences "$OUTPUT_DIR/${LABEL}_rep-seqs.qza"
  --o-denoising-stats "$OUTPUT_DIR/${LABEL}_denoising-stats.qza"
)

# Optional params
[ -n "$TRUNC_Q" ]       && CMD+=( --p-trunc-q "$TRUNC_Q" )
[ -n "$MAX_EE_F" ]      && CMD+=( --p-max-ee-f "$MAX_EE_F" )
[ -n "$MAX_EE_R" ]      && CMD+=( --p-max-ee-r "$MAX_EE_R" )
[ -n "$MIN_OVERLAP" ]   && CMD+=( --p-min-overlap "$MIN_OVERLAP" )
[ -n "$N_READS_LEARN" ] && CMD+=( --p-n-reads-learn "$N_READS_LEARN" )
[ -n "$HASHED_SEED" ]   && CMD+=( --p-hashed-seed "$HASHED_SEED" )

# Run
log "Running: ${CMD[*]}"
"${CMD[@]}"

# Summaries (optional)
if [ "$DO_SUMMARIES" -eq 1 ]; then
  log "Summarizing feature table → ${LABEL}_table.qzv"
  "$QIIME_CMD" feature-table summarize \
    --i-table "$OUTPUT_DIR/${LABEL}_table.qza" \
    --o-visualization "$OUTPUT_DIR/${LABEL}_table.qzv"

  log "Tabulating representative sequences → ${LABEL}_rep-seqs.qzv"
  "$QIIME_CMD" feature-table tabulate-seqs \
    --i-data "$OUTPUT_DIR/${LABEL}_rep-seqs.qza" \
    --o-visualization "$OUTPUT_DIR/${LABEL}_rep-seqs.qzv"
fi

log "Done. Outputs:"
log "  $OUTPUT_DIR/${LABEL}_table.qza"
log "  $OUTPUT_DIR/${LABEL}_rep-seqs.qza"
log "  $OUTPUT_DIR/${LABEL}_denoising-stats.qza"
[ "$DO_SUMMARIES" -eq 1 ] && {
  log "  $OUTPUT_DIR/${LABEL}_table.qzv"
  log "  $OUTPUT_DIR/${LABEL}_rep-seqs.qzv"
}



