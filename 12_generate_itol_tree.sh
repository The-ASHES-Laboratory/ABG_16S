##### 12_Generate_iTOL_tree
##### By Carter Clinton



#!/usr/bin/env bash
# generate_itol_tree.sh — Select top-N genera, prune ASVs & tree, and emit iTOL annotation files
# Maps original file: ITOL_trees_SILVA.txt
# Analysis Stage: Species-level & Targeted Trees
# Language: Bash

set -euo pipefail

# ---------------------------
# Logging / errors
# ---------------------------
log() { printf "[%s] %s\n" "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

usage() {
  cat <<'USAGE'
Select the top-N genera by total abundance, collect their member ASVs, prune the feature table,
representative sequences, and rooted tree accordingly, and export a Newick tree with iTOL
annotation files (color strip per genus).

Required:
  -t, --table-qza FILE         Feature table (.qza; FeatureTable[Frequency])
  -x, --taxonomy-qza FILE      Taxonomy (.qza; FeatureData[Taxonomy])
  -s, --rep-seqs-qza FILE      Representative sequences (.qza; FeatureData[Sequence])
  -r, --rooted-tree-qza FILE   Rooted phylogeny (.qza; Phylogeny[Rooted])
  -o, --output-dir DIR         Output directory

Optional:
  -l, --label STR              Basename label for outputs (default: itol)
  -n, --top-n INT              Number of genera to keep (default: 10)
  --level INT                  Taxonomy level for "genus" (default: 6; QIIME levels 1..7 = K..S)
  -q, --qiime-cmd PATH         'qiime' executable (default: qiime)
  -b, --biom-cmd PATH          'biom' executable (default: biom)  # used to convert BIOM → TSV
  --palette LIST               Comma-separated hex colors used cyclically for genera (default: preset 12-color palette)
  --keep-intermediates         Keep intermediary exports (default: off; intermediates removed)
  -h, --help                   Show this help

Outputs (under <OUTPUT_DIR>/itol/):
  tables/
    <label>_topG<GENUSLEVEL>_collapsed.qza          (collapsed L<level> table)
    <label>_top< N >_genera.txt                     (list of top-N genera by abundance)
    <label>_keep_feature_ids.txt                    (ASV IDs to keep)
    <label>_filtered_table.qza                      (table with only ASVs from top-N genera)
  seqs/
    <label>_filtered_rep-seqs.qza                   (rep seqs pruned to kept ASVs)
  tree/
    <label>_filtered_rooted-tree.qza                (phylogeny pruned to kept ASVs)
    <label>_filtered_rooted-tree.nwk                (exported Newick)
  itol/
    dataset_colorstrip.txt                          (iTOL color-strip dataset for leaves by genus)
    legend.txt                                      (legend mapping genus → color)
  taxonomy/
    taxonomy_export.tsv                             (exported taxonomy table; FeatureID→Taxon)
  qiime-info.txt

Example:
  prepare_itol_top_genera_tree.sh \
    -t <TABLE_QZA> -x <TAXONOMY_QZA> -s <REP_SEQS_QZA> -r <ROOTED_TREE_QZA> \
    -o <OUTPUT_DIR> -l runA -n 10 --level 6
USAGE
}

# ---------------------------
# Defaults
# ---------------------------
QIIME_CMD="${QIIME_CMD:-qiime}"
BIOM_CMD="${BIOM_CMD:-biom}"
TABLE_QZA=""
TAX_QZA=""
REP_SEQS_QZA=""
ROOTED_TREE_QZA=""
OUTPUT_DIR=""
LABEL="itol"
TOP_N=10
LEVEL=6
KEEP_INTERM=0
PALETTE_DEFAULT="#1f77b4,#ff7f0e,#2ca02c,#d62728,#9467bd,#8c564b,#e377c2,#7f7f7f,#bcbd22,#17becf,#a55194,#636363"
PALETTE="$PALETTE_DEFAULT"

# ---------------------------
# Parse CLI
# ---------------------------
OPTS_SHORT="t:x:s:r:o:l:n:q:b:h"
OPTS_LONG="table-qza:,taxonomy-qza:,rep-seqs-qza:,rooted-tree-qza:,output-dir:,label:,top-n:,qiime-cmd:,biom-cmd:,help,level:,palette:,keep-intermediates"
PARSED=$(getopt --options="$OPTS_SHORT" --longoptions="$OPTS_LONG" --name "$0" -- "$@") || { usage; exit 2; }
eval set -- "$PARSED"
while true; do
  case "$1" in
    -t|--table-qza)        TABLE_QZA="$2"; shift 2 ;;
    -x|--taxonomy-qza)     TAX_QZA="$2"; shift 2 ;;
    -s|--rep-seqs-qza)     REP_SEQS_QZA="$2"; shift 2 ;;
    -r|--rooted-tree-qza)  ROOTED_TREE_QZA="$2"; shift 2 ;;
    -o|--output-dir)       OUTPUT_DIR="$2"; shift 2 ;;
    -l|--label)            LABEL="$2"; shift 2 ;;
    -n|--top-n)            TOP_N="$2"; shift 2 ;;
    --level)               LEVEL="$2"; shift 2 ;;
    -q|--qiime-cmd)        QIIME_CMD="$2"; shift 2 ;;
    -b|--biom-cmd)         BIOM_CMD="$2"; shift 2 ;;
    --palette)             PALETTE="$2"; shift 2 ;;
    --keep-intermediates)  KEEP_INTERM=1; shift 1 ;;
    -h|--help)             usage; exit 0 ;;
    --) shift; break ;;
    *) die "Unknown option: $1" ;;
  esac
done

# ---------------------------
# Validate
# ---------------------------
[ -n "$TABLE_QZA" ]       || { usage; die "Missing --table-qza"; }
[ -n "$TAX_QZA" ]         || { usage; die "Missing --taxonomy-qza"; }
[ -n "$REP_SEQS_QZA" ]    || { usage; die "Missing --rep-seqs-qza"; }
[ -n "$ROOTED_TREE_QZA" ] || { usage; die "Missing --rooted-tree-qza"; }
[ -n "$OUTPUT_DIR" ]      || { usage; die "Missing --output-dir"; }
[ -f "$TABLE_QZA" ]       || die "Feature table not found: $TABLE_QZA"
[ -f "$TAX_QZA" ]         || die "Taxonomy artifact not found: $TAX_QZA"
[ -f "$REP_SEQS_QZA" ]    || die "Rep-seqs artifact not found: $REP_SEQS_QZA"
[ -f "$ROOTED_TREE_QZA" ] || die "Rooted tree artifact not found: $ROOTED_TREE_QZA"
[[ "$TOP_N" =~ ^[0-9]+$ ]] || die "--top-n must be an integer"
[[ "$LEVEL" =~ ^[1-7]$ ]]  || die "--level must be 1..7 (1=Kingdom .. 7=Species)"
mkdir -p "$OUTPUT_DIR"/{itol,tables,seqs,tree,taxonomy}
command -v "$QIIME_CMD" >/dev/null 2>&1 || die "Cannot find '$QIIME_CMD' in PATH. Activate QIIME 2."
command -v "$BIOM_CMD"  >/dev/null 2>&1 || die "Cannot find '$BIOM_CMD' in PATH. Install biom-format."

# ---------------------------
# Capture environment
# ---------------------------
"$QIIME_CMD" info > "$OUTPUT_DIR/qiime-info.txt" || true

# ---------------------------
# 1) Collapse table to genus (or requested) level and compute top-N genera
# ---------------------------
COLLAPSED="$OUTPUT_DIR/tables/${LABEL}_topG${LEVEL}_collapsed.qza"
log "Collapsing table to taxonomy level L${LEVEL} to compute top-${TOP_N} groups"
"$QIIME_CMD" taxa collapse \
  --i-table "$TABLE_QZA" \
  --i-taxonomy "$TAX_QZA" \
  --p-level "$LEVEL" \
  --o-collapsed-table "$COLLAPSED"

# Export collapsed table to BIOM/TSV
TMP_COLLAPSE_DIR="$(mktemp -d "$OUTPUT_DIR/tables/collapse_export.XXXXXX")"
"$QIIME_CMD" tools export --input-path "$COLLAPSED" --output-path "$TMP_COLLAPSE_DIR" >/dev/null
[ -f "$TMP_COLLAPSE_DIR/feature-table.biom" ] || die "Collapsed BIOM not found; export failed."
"$BIOM_CMD" convert -i "$TMP_COLLAPSE_DIR/feature-table.biom" -o "$TMP_COLLAPSE_DIR/table.tsv" --to-tsv >/dev/null

# Sum rows across samples to rank genera
TOP_LIST="$OUTPUT_DIR/tables/${LABEL}_top${TOP_N}_genera.txt"
awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=2;i<=NF;i++) hdr[i]=$i; next}
     {
       sum=0; for(i=2;i<=NF;i++) sum+=$i;
       print $1, sum
     }' "$TMP_COLLAPSE_DIR/table.tsv" \
  | sed '1d' \
  | sort -k2,2nr \
  | head -n "$TOP_N" \
  | cut -f1 \
  > "$TOP_LIST"
[ -s "$TOP_LIST" ] || die "Top genera list is empty."

log "Top-$TOP_N groups (taxonomy strings) written to: $TOP_LIST"

# ---------------------------
# 2) Export taxonomy.qza → TSV and build ASV→Genus map; collect ASV IDs in top-N genera
# ---------------------------
TAX_EXPORT_DIR="$OUTPUT_DIR/taxonomy/export"
mkdir -p "$TAX_EXPORT_DIR"
"$QIIME_CMD" tools export --input-path "$TAX_QZA" --output-path "$TAX_EXPORT_DIR" >/dev/null
[ -f "$TAX_EXPORT_DIR/taxonomy.tsv" ] || die "taxonomy.tsv not found in export."

# Build ASV→Genus mapping (assumes SILVA/Greengenes-style ranks 'k__...; p__...; ...; g__Genus; s__...')
MAP_ASV_GENUS="$OUTPUT_DIR/taxonomy/asv_to_genus.tsv"
awk -F'\t' '
  NR==1 { next } # skip header
  {
    asv=$1; tax=$2;
    n=split(tax, a, /;[[:space:]]*/);
    genus="";
    for(i=1; i<=n; i++){
      if (a[i] ~ /^g__/) { genus=a[i]; gsub(/^g__/, "", genus); break; }
    }
    if (genus=="") genus="Unassigned";
    print asv "\t" genus;
  }' "$TAX_EXPORT_DIR/taxonomy.tsv" > "$MAP_ASV_GENUS"

# Normalize top genera names (strip rank prefixes to match mapping)
TOP_GENUS_NAMES="$OUTPUT_DIR/tables/${LABEL}_top${TOP_N}_genera_clean.txt"
sed -E 's/^.*g__([^;|]+).*$/\1/' "$TOP_LIST" > "$TOP_GENUS_NAMES"

# Collect ASV IDs belonging to top-N genera
KEEP_IDS="$OUTPUT_DIR/tables/${LABEL}_keep_feature_ids.txt"
awk 'NR==FNR {top[$1]=1; next} ($2 in top) {print $1}' "$TOP_GENUS_NAMES" "$MAP_ASV_GENUS" > "$KEEP_IDS"
[ -s "$KEEP_IDS" ] || die "No ASVs map to the selected top genera."

# Build a minimal QIIME 2-style metadata TSV for IDs
IDS_TSV="$OUTPUT_DIR/tables/${LABEL}_keep_feature_ids.tsv"
printf "FeatureID\n" > "$IDS_TSV"
cat "$KEEP_IDS" >> "$IDS_TSV"

# ---------------------------
# 3) Filter table to kept ASVs; filter sequences and prune tree accordingly
# ---------------------------
FILTERED_TABLE="$OUTPUT_DIR/tables/${LABEL}_filtered_table.qza"
log "Filtering feature table to ASVs from top-$TOP_N genera"
"$QIIME_CMD" feature-table filter-features \
  --i-table "$TABLE_QZA" \
  --m-metadata-file "$IDS_TSV" \
  --o-filtered-table "$FILTERED_TABLE"

FILTERED_SEQS="$OUTPUT_DIR/seqs/${LABEL}_filtered_rep-seqs.qza"
log "Filtering representative sequences using filtered table"
"$QIIME_CMD" feature-table filter-seqs \
  --i-data "$REP_SEQS_QZA" \
  --i-table "$FILTERED_TABLE" \
  --o-filtered-data "$FILTERED_SEQS"

FILTERED_TREE="$OUTPUT_DIR/tree/${LABEL}_filtered_rooted-tree.qza"
log "Pruning rooted phylogeny to kept ASVs"
"$QIIME_CMD" phylogeny filter-tree \
  --i-tree "$ROOTED_TREE_QZA" \
  --i-table "$FILTERED_TABLE" \
  --o-filtered-tree "$FILTERED_TREE"

# Export Newick
TMP_TREE_EXPORT="$(mktemp -d "$OUTPUT_DIR/tree/export.XXXXXX")"
"$QIIME_CMD" tools export --input-path "$FILTERED_TREE" --output-path "$TMP_TREE_EXPORT" >/dev/null
if [ -f "$TMP_TREE_EXPORT/tree.nwk" ]; then
  mv "$TMP_TREE_EXPORT/tree.nwk" "$OUTPUT_DIR/tree/${LABEL}_filtered_rooted-tree.nwk"
else
  log "WARNING: Expected tree.nwk not found during export."
fi

# ---------------------------
# 4) Build iTOL dataset files (color strip per genus)
# ---------------------------
ITOL_DIR="$OUTPUT_DIR/itol"
ITOL_STRIP="$ITOL_DIR/dataset_colorstrip.txt"
ITOL_LEGEND="$ITOL_DIR/legend.txt"
mkdir -p "$ITOL_DIR"

# Prepare color palette array
IFS=',' read -r -a COLORS <<< "$PALETTE"

# Map genus → color (cyclic)
GENUS_COLORS="$ITOL_DIR/genus_colors.tsv"
awk 'NR==FNR{pal[NR]=$1; n=NR; next}
     { i=NR; color=pal[((i-1)%n)+1]; print $0 "\t" color }' \
  <(printf "%s\n" "${COLORS[@]}") \
  "$TOP_GENUS_NAMES" > "$GENUS_COLORS"

# Header per iTOL "DATASET_COLORSTRIP" format
cat > "$ITOL_STRIP" <<HDR
DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL\t${LABEL}_top${TOP_N}_genera
COLOR\t#000000
FIELD_LABELS\tGenus
FIELD_COLORS\t#000000
LEGEND_SHAPES\t1\t1\t1\t1\t1
LEGEND_COLORS\t$(cut -f2 "$GENUS_COLORS" | tr '\n' ',' | sed 's/,$//')
LEGEND_LABELS\t$(cut -f1 "$GENUS_COLORS" | tr '\n' ',' | sed 's/,$//')
DATA
HDR

# Add one line per leaf (ASV): <leafId> <color> <label>
# Join kept IDs with ASV→Genus map and then with genus colors
awk 'NR==FNR {gen[$1]=$2; next} NR==FNR+FNR {col[$1]=$2; next} { if($1 in gen && gen[$1] in col){ print $1 "\t" col[gen[$1]] "\t" gen[$1] } }' \
  "$MAP_ASV_GENUS" "$GENUS_COLORS" "$KEEP_IDS" >> "$ITOL_STRIP"

# Simple legend file
paste "$GENUS_COLORS" | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2}' > "$ITOL_LEGEND"

# ---------------------------
# 5) Cleanup intermediates
# ---------------------------
if [ "$KEEP_INTERM" -eq 0 ]; then
  rm -rf "$TMP_COLLAPSE_DIR" "$TMP_TREE_EXPORT" || true
fi

log "Done.
- Filtered Newick: $OUTPUT_DIR/tree/${LABEL}_filtered_rooted-tree.nwk
- iTOL color strip: $ITOL_STRIP
- Legend: $ITOL_LEGEND
- Kept ASVs: $KEEP_IDS
"



