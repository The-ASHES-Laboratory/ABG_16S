#!/usr/bin/env bash
# QIIME 2 Phylogenetic Tree Construction (MAFFT + FastTree)
# By Carter Clinton, Ph.D.
# build_phylogeny.sh — Build SILVA-based phylogeny via MAFFT → mask → FastTree → midpoint-root (QIIME 2)
# Maps original file: Phylogeentic_Reconstruction_SILVA.txt
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
Build a phylogeny for downstream diversity metrics using the QIIME 2 pipeline:
MAFFT alignment → alignment mask → FastTree → midpoint-root.

Required:
  -s, --rep-seqs-qza FILE        Representative sequences (.qza; FeatureData[Sequence])
  -o, --output-dir DIR           Output directory

Optional:
  -l, --label STR                Basename/label for outputs (default: tree)
  -q, --qiime-cmd PATH           'qiime' executable (default: qiime)
  -j, --threads INT              Threads for MAFFT alignment (default: 0 → auto)
  --mask-max-gap FLOAT           Max gap frequency for masking (default: 0.2)
  --mask-min-cons FLOAT          Min conservation for masking (default: 0.4)
  --export-newick                Also export rooted tree as Newick (.nwk)
  --keep-intermediates           Keep intermediate .qza files (default: keep; set to disable with --no-keep)
  --no-keep                      Remove intermediates (aligned/masked/unrooted) after finish
  -h, --help                     Show help

Outputs (written to <OUTPUT_DIR>/phylogeny/):
  <label>_aligned-rep-seqs.qza
  <label>_masked-alignment.qza
  <label>_unrooted-tree.qza
  <label>_rooted-tree.qza
  <label>_rooted-tree.nwk                (if --export-newick)
  qiime-info.txt

Examples:
  build_phylogeny_silva.sh \
    -s <REP_SEQS_QZA> \
    -o <OUTPUT_DIR> \
    -l runA -j 8 --mask-max-gap 0.2 --mask-min-cons 0.4 --export-newick
USAGE
}

# ---------------------------
# Defaults
# ---------------------------
QIIME_CMD="${QIIME_CMD:-qiime}"
REP_SEQS_QZA=""
OUTPUT_DIR=""
LABEL="tree"
THREADS=0
MASK_MAX_GAP=0.2
MASK_MIN_CONS=0.4
EXPORT_NEWICK=0
KEEP_INTERM=1

# ---------------------------
# Parse CLI
# ---------------------------
OPTS_SHORT="s:o:l:q:j:h"
OPTS_LONG="rep-seqs-qza:,output-dir:,label:,qiime-cmd:,threads:,help,mask-max-gap:,mask-min-cons:,export-newick,keep-intermediates,no-keep"
PARSED=$(getopt --options="$OPTS_SHORT" --longoptions="$OPTS_LONG" --name "$0" -- "$@") || { usage; exit 2; }
eval set -- "$PARSED"
while true; do
  case "$1" in
    -s|--rep-seqs-qza) REP_SEQS_QZA="$2"; shift 2 ;;
    -o|--output-dir)   OUTPUT_DIR="$2"; shift 2 ;;
    -l|--label)        LABEL="$2"; shift 2 ;;
    -q|--qiime-cmd)    QIIME_CMD="$2"; shift 2 ;;
    -j|--threads)      THREADS="$2"; shift 2 ;;
    --mask-max-gap)    MASK_MAX_GAP="$2"; shift 2 ;;
    --mask-min-cons)   MASK_MIN_CONS="$2"; shift 2 ;;
    --export-newick)   EXPORT_NEWICK=1; shift 1 ;;
    --keep-intermediates) KEEP_INTERM=1; shift 1 ;;
    --no-keep)         KEEP_INTERM=0; shift 1 ;;
    -h|--help)         usage; exit 0 ;;
    --) shift; break ;;
    *) die "Unknown option: $1" ;;
  esac
done

# ---------------------------
# Validate
# ---------------------------
[ -n "$REP_SEQS_QZA" ] || { usage; die "Missing --rep-seqs-qza"; }
[ -n "$OUTPUT_DIR" ]   || { usage; die "Missing --output-dir"; }
[ -f "$REP_SEQS_QZA" ] || die "Rep seqs artifact not found: $REP_SEQS_QZA"
mkdir -p "$OUTPUT_DIR/phylogeny"
command -v "$QIIME_CMD" >/dev/null 2>&1 || die "Cannot find '$QIIME_CMD' in PATH. Activate QIIME 2."

# ---------------------------
# Capture environment
# ---------------------------
"$QIIME_CMD" info > "$OUTPUT_DIR/qiime-info.txt" || true

# ---------------------------
# Run align-to-tree pipeline
# ---------------------------
ALIGNED="$OUTPUT_DIR/phylogeny/${LABEL}_aligned-rep-seqs.qza"
MASKED="$OUTPUT_DIR/phylogeny/${LABEL}_masked-alignment.qza"
UNROOTED="$OUTPUT_DIR/phylogeny/${LABEL}_unrooted-tree.qza"
ROOTED="$OUTPUT_DIR/phylogeny/${LABEL}_rooted-tree.qza"

CMD=( "$QIIME_CMD" phylogeny align-to-tree-mafft-fasttree
      --i-sequences "$REP_SEQS_QZA"
      --o-alignment "$ALIGNED"
      --o-masked-alignment "$MASKED"
      --o-tree "$UNROOTED"
      --o-rooted-tree "$ROOTED"
      --p-mask-max-gap-frequency "$MASK_MAX_GAP"
      --p-mask-min-conservation "$MASK_MIN_CONS" )

# MAFFT threads; 0 means "auto" in many QIIME releases
[ "$THREADS" -ge 0 ] && CMD+=( --p-n-threads "$THREADS" )

log "Running: ${CMD[*]}"
"${CMD[@]}"

# ---------------------------
# (Optional) Export Newick
# ---------------------------
if [ "$EXPORT_NEWICK" -eq 1 ]; then
  TMPDIR="$(mktemp -d "$OUTPUT_DIR/phylogeny/export.XXXXXX")"
  log "Exporting rooted tree to Newick (.nwk)"
  "$QIIME_CMD" tools export \
    --input-path "$ROOTED" \
    --output-path "$TMPDIR"
  if [ -f "$TMPDIR/tree.nwk" ]; then
    mv "$TMPDIR/tree.nwk" "$OUTPUT_DIR/phylogeny/${LABEL}_rooted-tree.nwk"
  else
    log "WARNING: Expected tree.nwk not found in export directory."
  fi
  rm -rf "$TMPDIR" || true
fi

# ---------------------------
# Cleanup (optional)
# ---------------------------
if [ "$KEEP_INTERM" -eq 0 ]; then
  log "Removing intermediate artifacts (--no-keep enabled)"
  rm -f "$ALIGNED" "$MASKED" "$UNROOTED" || true
fi

log "Done. Outputs in $OUTPUT_DIR/phylogeny:"
log "  $(basename "$ALIGNED")"
log "  $(basename "$MASKED")"
log "  $(basename "$UNROOTED")"
log "  $(basename "$ROOTED")"
[ "$EXPORT_NEWICK" -eq 1 ] && log "  ${LABEL}_rooted-tree.nwk"
