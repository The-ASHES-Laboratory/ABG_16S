##### 05_Taxonomy
##### By Carter Clinton



#!/usr/bin/env bash
# taxonomy.sh  —  Train (optional) and run SILVA Naïve Bayes taxonomy classification; make bar plots
# Maps original file: taxonomy_with_silva.txt
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
Train (optionally) a Naïve Bayes classifier for the 515F/806R region on SILVA (or another reference)
and classify ASVs with QIIME 2 'feature-classifier classify-sklearn'. Optionally generate taxa bar plots.

Required:
  -t, --table-qza FILE           Feature table (.qza; FeatureTable[Frequency])
  -s, --rep-seqs-qza FILE        Representative sequences (.qza; FeatureData[Sequence])
  -o, --output-dir DIR           Output directory

One of the following classifier sources:
  (A) Use a pre-trained classifier:
      -c, --classifier-qza FILE  Pre-trained classifier (.qza)

  (B) Train a classifier from reference:
      --ref-seqs-qza FILE        Reference sequences (.qza; FeatureData[Sequence])
      --ref-tax-qza FILE         Reference taxonomy (.qza; FeatureData[Taxonomy])
      --f-primer STR             Forward primer sequence (default: GTGCCAGCMGCCGCGGTAA  # 515F)
      --r-primer STR             Reverse primer sequence (default: GGACTACHVGGGTWTCTAAT # 806R)
      --min-length INT           Min amplicon length after extraction (optional)
      --max-length INT           Max amplicon length after extraction (optional)
      --train-label STR          Basename for trained classifier (default: silva_515f_806r)

Optional:
  -m, --metadata FILE            Sample metadata TSV (for taxa bar plots)
  -l, --label STR                Basename for outputs (default: tax)
  -q, --qiime-cmd PATH           qiime executable (default: qiime)
  -j, --threads INT              Threads for classify-sklearn (default: 0 → auto)
  --confidence FLOAT             Naïve Bayes confidence threshold (default: 0.7)
  --skip-barplots                Skip taxa bar plots (default: generate if metadata is provided)
  -h, --help                     Show this help

Outputs (in <OUTPUT_DIR>):
  taxonomy/<label>_taxonomy.qza
  taxonomy/<label>_taxonomy.qzv
  barplots/<label>_taxa-bar-plots.qzv          (if metadata provided and barplots not skipped)
  trained/<train-label>_classifier.qza         (if training performed)
  qiime-info.txt

Examples:
  # A) Use a pre-trained SILVA 138 classifier
  taxonomy_classify_silva.sh \
    -t <TABLE_QZA> \
    -s <REP_SEQS_QZA> \
    -c <SILVA_CLASSIFIER_QZA> \
    -m <METADATA_TSV> \
    -o <OUTPUT_DIR> -l runA --confidence 0.7 -j 8

  # B) Train a region-specific classifier from reference seqs + taxonomy
  taxonomy_classify_silva.sh \
    -t <TABLE_QZA> -s <REP_SEQS_QZA> -o <OUTPUT_DIR> -l runA \
    --ref-seqs-qza <REFERENCE_SEQS_QZA> \
    --ref-tax-qza  <REFERENCE_TAX_QZA> \
    --f-primer GTGCCAGCMGCCGCGGTAA \
    --r-primer GGACTACHVGGGTWTCTAAT \
    --min-length 100 --max-length 500
USAGE
}

# ---------------------------
# Defaults
# ---------------------------
QIIME_CMD="${QIIME_CMD:-qiime}"
TABLE_QZA=""
REP_SEQS_QZA=""
OUTPUT_DIR=""
CLASSIFIER_QZA=""

REF_SEQS_QZA=""
REF_TAX_QZA=""
F_PRIMER="GTGCCAGCMGCCGCGGTAA"
R_PRIMER="GGACTACHVGGGTWTCTAAT"
MIN_LEN=""
MAX_LEN=""
TRAIN_LABEL="silva_515f_806r"

METADATA=""
LABEL="tax"
THREADS=0
CONFIDENCE=0.7
SKIP_BARPLOTS=0

# ---------------------------
# Parse CLI
# ---------------------------
OPTS_SHORT="t:s:o:c:m:l:q:j:h"
OPTS_LONG="table-qza:,rep-seqs-qza:,output-dir:,classifier-qza:,metadata:,label:,qiime-cmd:,threads:,help,ref-seqs-qza:,ref-tax-qza:,f-primer:,r-primer:,min-length:,max-length:,train-label:,confidence:,skip-barplots"
PARSED=$(getopt --options="$OPTS_SHORT" --longoptions="$OPTS_LONG" --name "$0" -- "$@") || { usage; exit 2; }
eval set -- "$PARSED"
while true; do
  case "$1" in
    -t|--table-qza)       TABLE_QZA="$2"; shift 2 ;;
    -s|--rep-seqs-qza)    REP_SEQS_QZA="$2"; shift 2 ;;
    -o|--output-dir)      OUTPUT_DIR="$2"; shift 2 ;;
    -c|--classifier-qza)  CLASSIFIER_QZA="$2"; shift 2 ;;
    --ref-seqs-qza)       REF_SEQS_QZA="$2"; shift 2 ;;
    --ref-tax-qza)        REF_TAX_QZA="$2"; shift 2 ;;
    --f-primer)           F_PRIMER="$2"; shift 2 ;;
    --r-primer)           R_PRIMER="$2"; shift 2 ;;
    --min-length)         MIN_LEN="$2"; shift 2 ;;
    --max-length)         MAX_LEN="$2"; shift 2 ;;
    --train-label)        TRAIN_LABEL="$2"; shift 2 ;;
    -m|--metadata)        METADATA="$2"; shift 2 ;;
    -l|--label)           LABEL="$2"; shift 2 ;;
    -q|--qiime-cmd)       QIIME_CMD="$2"; shift 2 ;;
    -j|--threads)         THREADS="$2"; shift 2 ;;
    --confidence)         CONFIDENCE="$2"; shift 2 ;;
    --skip-barplots)      SKIP_BARPLOTS=1; shift 1 ;;
    -h|--help)            usage; exit 0 ;;
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
mkdir -p "$OUTPUT_DIR"/{taxonomy,barplots,trained}
command -v "$QIIME_CMD" >/dev/null 2>&1 || die "Cannot find '$QIIME_CMD' in PATH. Activate your QIIME 2 env."

# Classifier source validation
if [ -n "$CLASSIFIER_QZA" ]; then
  [ -f "$CLASSIFIER_QZA" ] || die "Classifier not found: $CLASSIFIER_QZA"
else
  # Require training inputs
  [ -n "$REF_SEQS_QZA" ] && [ -n "$REF_TAX_QZA" ] || die "Provide --classifier-qza OR both --ref-seqs-qza and --ref-tax-qza to train."
  [ -f "$REF_SEQS_QZA" ] || die "Reference sequences not found: $REF_SEQS_QZA"
  [ -f "$REF_TAX_QZA" ]  || die "Reference taxonomy not found: $REF_TAX_QZA"
fi

# ---------------------------
# Capture environment
# ---------------------------
"$QIIME_CMD" info > "$OUTPUT_DIR/qiime-info.txt" || true

# ---------------------------
# (Optional) Train classifier
# ---------------------------
if [ -z "$CLASSIFIER_QZA" ]; then
  log "Training region-specific classifier with primers F:'$F_PRIMER' R:'$R_PRIMER'"
  TRAIN_EXTRACT="$OUTPUT_DIR/trained/${TRAIN_LABEL}_refseqs_extracted.qza"
  TRAIN_CLASSIFIER="$OUTPUT_DIR/trained/${TRAIN_LABEL}_classifier.qza"

  EXTRACT_CMD=( "$QIIME_CMD" feature-classifier extract-reads
    --i-sequences "$REF_SEQS_QZA"
    --p-f-primer "$F_PRIMER"
    --p-r-primer "$R_PRIMER"
    --o-reads "$TRAIN_EXTRACT" )
  [ -n "$MIN_LEN" ] && EXTRACT_CMD+=( --p-min-length "$MIN_LEN" )
  [ -n "$MAX_LEN" ] && EXTRACT_CMD+=( --p-max-length "$MAX_LEN" )
  log "Running: ${EXTRACT_CMD[*]}"
  "${EXTRACT_CMD[@]}"

  log "Fitting Naïve Bayes classifier"
  "$QIIME_CMD" feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$TRAIN_EXTRACT" \
    --i-reference-taxonomy "$REF_TAX_QZA" \
    --o-classifier "$TRAIN_CLASSIFIER"

  CLASSIFIER_QZA="$TRAIN_CLASSIFIER"
fi

# ---------------------------
# Classify ASVs
# ---------------------------
TAX_QZA="$OUTPUT_DIR/taxonomy/${LABEL}_taxonomy.qza"
TAX_QZV="$OUTPUT_DIR/taxonomy/${LABEL}_taxonomy.qzv"
CLASSIFY_CMD=( "$QIIME_CMD" feature-classifier classify-sklearn
  --i-classifier "$CLASSIFIER_QZA"
  --i-reads "$REP_SEQS_QZA"
  --p-confidence "$CONFIDENCE"
  --o-classification "$TAX_QZA" )
[ "$THREADS" -gt 0 ] && CLASSIFY_CMD+=( --p-n-jobs "$THREADS" )
log "Running: ${CLASSIFY_CMD[*]}"
"${CLASSIFY_CMD[@]}"

log "Tabulating taxonomy → $TAX_QZV"
"$QIIME_CMD" metadata tabulate \
  --m-input-file "$TAX_QZA" \
  --o-visualization "$TAX_QZV"

# ---------------------------
# (Optional) Bar plots
# ---------------------------
if [ "$SKIP_BARPLOTS" -eq 0 ] && [ -n "$METADATA" ]; then
  [ -f "$METADATA" ] || die "Metadata not found: $METADATA"
  BARS_QZV="$OUTPUT_DIR/barplots/${LABEL}_taxa-bar-plots.qzv"
  log "Generating taxa bar plots → $BARS_QZV"
  "$QIIME_CMD" taxa barplot \
    --i-table "$TABLE_QZA" \
    --i-taxonomy "$TAX_QZA" \
    --m-metadata-file "$METADATA" \
    --o-visualization "$BARS_QZV"
else
  log "Skipping taxa bar plots (no metadata or --skip-barplots)."
fi

log "Done. Artifacts:"
log "  TAXONOMY_QZA: $TAX_QZA"
log "  TAXONOMY_QZV: $TAX_QZV"
[ -f "${OUTPUT_DIR}/barplots/${LABEL}_taxa-bar-plots.qzv" ] && log "  BARPLOTS_QZV: ${OUTPUT_DIR}/barplots/${LABEL}_taxa-bar-plots.qzv"
[ -d "$OUTPUT_DIR/trained" ] && ls "$OUTPUT_DIR/trained"/*.qza >/dev/null 2>&1 && log "  TRAINED_CLASSIFIER: $(ls "$OUTPUT_DIR/trained"/*_classifier.qza)"



