##### 09_Permdisp
##### By Carter Clinton



#!/usr/bin/env bash
# permdisp.sh — Test homogeneity of dispersion (PERMDISP) across groups (beta) for multiple distance matrices
# Maps original file: 07_Diversity_Analyses_Add-on.txt
# Analysis Stage: Phylogeny & Diversity
# Language: Bash

set -euo pipefail

log() { printf "[%s] %s\n" "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

usage() {
  cat <<'USAGE'
Run PERMDISP (beta-group-significance --p-method permdisp) for one or more
distance matrices against one or more metadata columns to assess homogeneity
of dispersion (a prerequisite check when interpreting PERMANOVA).

Required:
  -m, --metadata FILE            QIIME 2 metadata TSV
  -c, --columns LIST             Comma-separated metadata columns to test (e.g., Group,Sex)
  -o, --output-dir DIR           Output directory for results

Provide distance matrices via ONE of the following:
  (A) From a core-metrics directory/prefix (auto-discover standard matrices):
      -d, --core-metrics-dir DIR   Directory containing core metrics outputs
      -l, --label STR              Prefix used when generating core metrics (e.g., 'div' or 'runA')
         # Expected filenames in <DIR>:
         #   <label>_unweighted_unifrac_dm.qza
         #   <label>_weighted_unifrac_dm.qza
         #   <label>_jaccard_dm.qza
         #   <label>_bray_curtis_dm.qza

  (B) Explicit list of matrices (comma-separated .qza paths):
      --dm-list LIST              e.g., /path/uuf_dm.qza,/path/wuf_dm.qza

Optional:
  -q, --qiime-cmd PATH           'qiime' executable (default: qiime)
  --permutations INT             Number of permutations (default: 999)
  --pairwise                     Also perform pairwise post-hoc tests (default: off)
  --strata-col STR               Metadata column for stratified permutations (optional)
  --keep-names LIST              Comma-separated names to keep (filter which matrices to run; names must be one of:
                                 unweighted_unifrac, weighted_unifrac, jaccard, bray_curtis; default: all found)
  --name-aliases LIST            Optional aliases for matrices in --dm-list (comma-separated; same length as --dm-list).
  -h, --help                     Show help

Outputs (under <OUTPUT_DIR>):
  permdisp/<COLUMN>/<NAME>_permdisp.qzv
  qiime-info.txt

Examples:
  # Auto-discover matrices from core-metrics output
  beta_permdisp_tests.sh \
    -m <METADATA_TSV> \
    -c Group,Sex \
    -o <OUTPUT_DIR> \
    -d <CORE_METRICS_DIR> -l runA --permutations 999 --pairwise

  # Explicit list of matrices with custom names
  beta_permdisp_tests.sh \
    -m <METADATA_TSV> \
    -c Group \
    -o <OUTPUT_DIR> \
    --dm-list <UUF_DM_QZA>,<WUF_DM_QZA> \
    --name-aliases unweighted_unifrac,weighted_unifrac
USAGE
}

# Defaults
QIIME_CMD="${QIIME_CMD:-qiime}"
METADATA=""
COLUMNS=""
OUTPUT_DIR=""
CORE_DIR=""
LABEL=""
DM_LIST=""
NAME_ALIASES=""
KEEP_NAMES=""
PERMS=999
PAIRWISE=0
STRATA_COL=""

# Parse CLI
OPTS_SHORT="m:c:o:d:l:q:h"
OPTS_LONG="metadata:,columns:,output-dir:,core-metrics-dir:,label:,qiime-cmd:,help,dm-list:,name-aliases:,permutations:,pairwise,strata-col:,keep-names:"
PARSED=$(getopt --options="$OPTS_SHORT" --longoptions="$OPTS_LONG" --name "$0" -- "$@") || { usage; exit 2; }
eval set -- "$PARSED"
while true; do
  case "$1" in
    -m|--metadata)          METADATA="$2"; shift 2 ;;
    -c|--columns)           COLUMNS="$2"; shift 2 ;;
    -o|--output-dir)        OUTPUT_DIR="$2"; shift 2 ;;
    -d|--core-metrics-dir)  CORE_DIR="$2"; shift 2 ;;
    -l|--label)             LABEL="$2"; shift 2 ;;
    --dm-list)              DM_LIST="$2"; shift 2 ;;
    --name-aliases)         NAME_ALIASES="$2"; shift 2 ;;
    --permutations)         PERMS="$2"; shift 2 ;;
    --pairwise)             PAIRWISE=1; shift 1 ;;
    --strata-col)           STRATA_COL="$2"; shift 2 ;;
    --keep-names)           KEEP_NAMES="$2"; shift 2 ;;
    -q|--qiime-cmd)         QIIME_CMD="$2"; shift 2 ;;
    -h|--help)              usage; exit 0 ;;
    --) shift; break ;;
    *) die "Unknown option: $1" ;;
  esac
done

# Validate required
[ -n "$METADATA" ]   || { usage; die "Missing --metadata"; }
[ -n "$COLUMNS" ]    || { usage; die "Missing --columns"; }
[ -n "$OUTPUT_DIR" ] || { usage; die "Missing --output-dir"; }
[ -f "$METADATA" ]   || die "Metadata TSV not found: $METADATA"
mkdir -p "$OUTPUT_DIR/permdisp"
command -v "$QIIME_CMD" >/dev/null 2>&1 || die "Cannot find '$QIIME_CMD' in PATH. Activate QIIME 2."
[[ "$PERMS" =~ ^[0-9]+$ ]] || die "--permutations must be an integer"

# Ensure one source of DMs provided
if [ -n "$CORE_DIR" ] && [ -n "$DM_LIST" ]; then
  die "Provide either --core-metrics-dir OR --dm-list, not both."
fi
if [ -z "$CORE_DIR" ] && [ -z "$DM_LIST" ]; then
  die "Provide one of --core-metrics-dir or --dm-list."
fi

# Capture environment
"$QIIME_CMD" info > "$OUTPUT_DIR/qiime-info.txt" || true

declare -A DMS
# Load DMs from core-metrics dir
if [ -n "$CORE_DIR" ]; then
  [ -d "$CORE_DIR" ] || die "Core metrics directory not found: $CORE_DIR"
  [ -n "$LABEL" ] || die "When using --core-metrics-dir you must provide --label"
  DMS["unweighted_unifrac"]="$CORE_DIR/${LABEL}_unweighted_unifrac_dm.qza"
  DMS["weighted_unifrac"]="$CORE_DIR/${LABEL}_weighted_unifrac_dm.qza"
  DMS["jaccard"]="$CORE_DIR/${LABEL}_jaccard_dm.qza"
  DMS["bray_curtis"]="$CORE_DIR/${LABEL}_bray_curtis_dm.qza"
fi

# Load DMs from explicit list
if [ -n "$DM_LIST" ]; then
  IFS=',' read -r -a dm_arr <<< "$DM_LIST"
  if [ -n "$NAME_ALIASES" ]; then
    IFS=',' read -r -a nm_arr <<< "$NAME_ALIASES"
    [ "${#dm_arr[@]}" -eq "${#nm_arr[@]}" ] || die "--name-aliases must have the same number of items as --dm-list"
    for i in "${!dm_arr[@]}"; do
      name="${nm_arr[$i]}"
      path="${dm_arr[$i]}"
      DMS["$name"]="$path"
    done
  else
    # Assign generic names dm1, dm2, ... if none supplied
    idx=1
    for path in "${dm_arr[@]}"; do
      DMS["dm${idx}"]="$path"; idx=$((idx+1))
    done
  fi
fi

# Optionally keep only selected names
if [ -n "$KEEP_NAMES" ]; then
  IFS=',' read -r -a keep_arr <<< "$KEEP_NAMES"
  declare -A keep_set
  for n in "${keep_arr[@]}"; do keep_set["$n"]=1; done
  for n in "${!DMS[@]}"; do
    if [ -z "${keep_set[$n]+x}" ]; then unset "DMS[$n]"; fi
  done
fi

# Validate DM files
for name in "${!DMS[@]}"; do
  path="${DMS[$name]}"
  if [ ! -f "$path" ]; then
    log "WARNING: Distance matrix not found for '$name': $path (skipping)"
    unset "DMS[$name]"
  fi
done
[ "${#DMS[@]}" -gt 0 ] || die "No valid distance matrices available after validation."

# Iterate over metadata columns and matrices
IFS=',' read -r -a COLS <<< "$COLUMNS"
for col in "${COLS[@]}"; do
  COL_TRIM="$(echo "$col" | sed 's/^ *//;s/ *$//')"
  SAFE_COL="$(echo "$COL_TRIM" | tr ' ' '_' )"
  OUT_COL_DIR="$OUTPUT_DIR/permdisp/$SAFE_COL"
  mkdir -p "$OUT_COL_DIR"

  for name in "${!DMS[@]}"; do
    dm="${DMS[$name]}"
    base="$OUT_COL_DIR/${name}_permdisp"
    args=( --i-distance-matrix "$dm"
           --m-metadata-file "$METADATA"
           --m-metadata-column "$COL_TRIM"
           --p-permutations "$PERMS"
           --p-method permdisp
           --o-visualization "${base}.qzv" )
    [ "$PAIRWISE" -eq 1 ] && args+=( --p-pairwise )
    [ -n "$STRATA_COL" ] && args+=( --p-strata "$STRATA_COL" )

    log "PERMDISP on '$name' for column '$COL_TRIM' → ${base}.qzv"
    "$QIIME_CMD" diversity beta-group-significance "${args[@]}"
  done
done

log "Done. PERMDISP outputs in: $OUTPUT_DIR/permdisp"

