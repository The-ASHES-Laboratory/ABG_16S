##### 04_Rarefaction
##### By Carter Clinton



#!/usr/bin/env bash
# rarefaction.sh — Generate alpha-rarefaction curves (QIIME 2)
# Analysis Stage: QC & Denoising
# Language: Bash
#
# Description:
#   Wraps `qiime diversity alpha-rarefaction` to create alpha-rarefaction curves
#   from an unrarefied FeatureTable[Frequency]. Optionally include a rooted
#   phylogeny to compute Faith’s PD. All paths/parameters are CLI-configurable.

set -euo pipefail

log() { printf "[%s] %s\n" "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

usage() {
  cat <<'USAGE'
Generate alpha-rarefaction curves from an unrarefied feature table.

Required:
  -i, --input-table FILE     Input feature table (.qza; FeatureTable[Frequency])
  -o, --output-dir  DIR      Output directory

Optional:
  -t, --tree FILE            Rooted phylogeny (.qza; Phylogeny[Rooted]) to include Faith's PD
  -m, --metadata FILE        Sample metadata TSV to overlay in the visualization
  -l, --label STR            Basename/label for outputs (default: rare)
  --max-depth INT            Maximum rarefaction depth (default: 1000)
  --min-depth INT            Minimum rarefaction depth (default: 1)
  --steps INT                Number of depth steps between min and max (default: 10)
  --iterations INT           Iterations per depth (default: 10)
  -q, --qiime-cmd PATH       qiime executable (default: qiime)
  --no-metadata              Do not attach metadata even if provided
  -h, --help                 Show this help

Outputs (written to <OUTPUT_DIR>):
  <label>_alpha_rarefaction.qzv
  qiime-info.txt

Examples:
  alpha_rarefaction.sh \
    -i <TABLE_QZA> \
    -o <OUTPUT_DIR> \
    --max-depth 1000 \
    -m <METADATA_TSV>

  # Include phylogeny to add Faith's PD curves
  alpha_rarefaction.sh \
    -i <TABLE_QZA> -o <OUTPUT_DIR> -t <ROOTED_TREE_QZA> --max-depth 1500
USAGE
}

# Defaults
QIIME_CMD="${QIIME_CMD:-qiime}"
INPUT_TABLE=""
ROOTED_TREE=""
OUTPUT_DIR=""
METADATA=""
LABEL="rare"
MAX_DEPTH=1000
MIN_DEPTH=1
STEPS=10
ITERATIONS=10
ATTACH_META=1

# Parse args
OPTS_SHORT="i:o:t:m:l:q:h"
OPTS_LONG="input-table:,output-dir:,tree:,metadata:,label:,qiime-cmd:,help,max-depth:,min-depth:,steps:,iterations:,no-metadata"
PARSED=$(getopt --options="$OPTS_SHORT" --longoptions="$OPTS_LONG" --name "$0" -- "$@") || { usage; exit 2; }
eval set -- "$PARSED"
while true; do
  case "$1" in
    -i|--input-table) INPUT_TABLE="$2"; shift 2 ;;
    -o|--output-dir)  OUTPUT_DIR="$2"; shift 2 ;;
    -t|--tree)        ROOTED_TREE="$2"; shift 2 ;;
    -m|--metadata)    METADATA="$2"; shift 2 ;;
    -l|--label)       LABEL="$2"; shift 2 ;;
    -q|--qiime-cmd)   QIIME_CMD="$2"; shift 2 ;;
    --max-depth)      MAX_DEPTH="$2"; shift 2 ;;
    --min-depth)      MIN_DEPTH="$2"; shift 2 ;;
    --steps)          STEPS="$2"; shift 2 ;;
    --iterations)     ITERATIONS="$2"; shift 2 ;;
    --no-metadata)    ATTACH_META=0; shift 1 ;;
    -h|--help)        usage; exit 0 ;;
    --) shift; break ;;
    *) die "Unknown option: $1" ;;
  esac
done

# Validate
[ -n "$INPUT_TABLE" ] || { usage; die "Missing --input-table"; }
[ -n "$OUTPUT_DIR" ]  || { usage; die "Missing --output-dir"; }
[ -f "$INPUT_TABLE" ] || die "Input table not found: $INPUT_TABLE"
[[ "$MAX_DEPTH" =~ ^[0-9]+$ ]] || die "--max-depth must be integer"
[[ "$MIN_DEPTH" =~ ^[0-9]+$ ]] || die "--min-depth must be integer"
[[ "$STEPS" =~ ^[0-9]+$ ]]     || die "--steps must be integer"
[[ "$ITERATIONS" =~ ^[0-9]+$ ]]|| die "--iterations must be integer"
[ "$MIN_DEPTH" -le "$MAX_DEPTH" ] || die "--min-depth must be <= --max-depth"
[ -n "$ROOTED_TREE" ] && [ ! -f "$ROOTED_TREE" ] && die "Tree not found: $ROOTED_TREE"
[ -n "$METADATA" ] && [ ! -f "$METADATA" ] && die "Metadata not found: $METADATA"
mkdir -p "$OUTPUT_DIR"
command -v "$QIIME_CMD" >/dev/null 2>&1 || die "Cannot find '$QIIME_CMD' in PATH. Activate QIIME 2."

# Capture environment
"$QIIME_CMD" info > "$OUTPUT_DIR/qiime-info.txt" || true

# Build command
OUT_QZV="$OUTPUT_DIR/${LABEL}_alpha_rarefaction.qzv"
CMD=( "$QIIME_CMD" diversity alpha-rarefaction
      --i-table "$INPUT_TABLE"
      --p-max-depth "$MAX_DEPTH"
      --p-min-depth "$MIN_DEPTH"
      --p-steps "$STEPS"
      --p-iterations "$ITERATIONS"
      --o-visualization "$OUT_QZV" )

[ -n "$ROOTED_TREE" ] && CMD+=( --i-phylogeny "$ROOTED_TREE" )
if [ "$ATTACH_META" -eq 1 ] && [ -n "$METADATA" ]; then
  CMD+=( --m-metadata-file "$METADATA" )
fi

log "Running: ${CMD[*]}"
"${CMD[@]}"

log "Done. Outputs:"
log "  $OUT_QZV"

