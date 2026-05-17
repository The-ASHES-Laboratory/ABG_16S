#!/usr/bin/env bash
# Ancient DNA Damage Profiling with DamageProfiler
# By Carter Clinton, Ph.D.
# adna_damage_profiler.sh — Map short reads and generate aDNA damage profiles with DamageProfiler
# Maps original file: aDNA_verification.txt
# Analysis Stage: aDNA Authentication
# Language: Bash

set -euo pipefail

# ---------------------------
# Logging / errors
# ---------------------------
log(){ printf "[%s] %s\n" "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "$*" >&2; }
die(){ log "ERROR: $*"; exit 1; }

usage(){
  cat <<'USAGE'
Align paired-end reads to a reference genome and compute **ancient DNA** damage profiles
using **DamageProfiler**. Designed for short, potentially damaged reads (e.g., 16S or
shotgun fragments). Defaults to BWA **aln** with seed disabled, suitable for aDNA.

Required:
  -r, --reference-fasta FILE         Reference FASTA
  -i, --reads-dir DIR                Directory containing FASTQ(.gz) files
  -o, --output-dir DIR               Destination directory

Optional:
  -l, --label STR                    Run label (used in output folder names; default: adna)
  --r1-suffix STR                    R1 suffix pattern (default: _R1_001.fastq.gz)
  --r2-suffix STR                    R2 suffix pattern (default: _R2_001.fastq.gz)
  --threads INT                      Threads for BWA/SAMtools (default: 4)
  --bwa-mode {aln,mem}               Align with 'aln' (aDNA-friendly; default) or 'mem'
  # aln tuning:
  --aln-seed-length INT              BWA aln -l (seed length; 1024 disables seeding; default: 1024)
  --aln-max-diff FLOAT               BWA aln -n (max fraction of differences; default: 0.01)
  --aln-opens INT                    BWA aln -o (gap opens; default: 2)
  --min-mapq INT                     Minimum MAPQ to keep (default: 30)
  --dedup                            Mark duplicates (samtools markdup)
  --damageprofiler-cmd PATH          'damageprofiler' executable (default: damageprofiler)
  --damageprofiler-jar FILE          Use JAR instead of CLI (requires 'java')
  --resume                           Skip samples with existing BAM + damageprofiler output
  -h, --help                         Show this help

Outputs (under <OUTPUT_DIR>/):
  bams/<sample>.bam(.bai)            Sorted, indexed BAM per sample
  damageprofiler/<sample>/           DamageProfiler outputs per sample
  manifests/run_manifest.tsv         Per-sample summary (paths, read counts)
  logs/                              Alignment/log files

Example:
  adna_damageprofiler_pipeline.sh \
    -r <REF_FASTA> -i <READS_DIR> -o <OUT_DIR> -l runB --threads 8 --dedup
USAGE
}

# ---------------------------
# Defaults
# ---------------------------
REF=""
READS_DIR=""
OUT_DIR=""
LABEL="adna"
R1_SUFFIX="_R1_001.fastq.gz"
R2_SUFFIX="_R2_001.fastq.gz"
THREADS=4
BWA_MODE="aln"
ALN_SEED=1024
ALN_N="0.01"
ALN_O=2
MIN_MAPQ=30
DEDUP=0
DP_CMD="${DAMAGEPROFILER_CMD:-damageprofiler}"
DP_JAR="${DAMAGEPROFILER_JAR:-}"
RESUME=0

# ---------------------------
# Parse CLI
# ---------------------------
OPTS_SHORT="r:i:o:l:h"
OPTS_LONG="reference-fasta:,reads-dir:,output-dir:,label:,help,r1-suffix:,r2-suffix:,threads:,bwa-mode:,aln-seed-length:,aln-max-diff:,aln-opens:,min-mapq:,dedup,damageprofiler-cmd:,damageprofiler-jar:,resume"
PARSED=$(getopt --options="$OPTS_SHORT" --longoptions="$OPTS_LONG" --name "$0" -- "$@") || { usage; exit 2; }
eval set -- "$PARSED"
while true; do
  case "$1" in
    -r|--reference-fasta) REF="$2"; shift 2 ;;
    -i|--reads-dir)       READS_DIR="$2"; shift 2 ;;
    -o|--output-dir)      OUT_DIR="$2"; shift 2 ;;
    -l|--label)           LABEL="$2"; shift 2 ;;
    --r1-suffix)          R1_SUFFIX="$2"; shift 2 ;;
    --r2-suffix)          R2_SUFFIX="$2"; shift 2 ;;
    --threads)            THREADS="$2"; shift 2 ;;
    --bwa-mode)           BWA_MODE="$2"; shift 2 ;;
    --aln-seed-length)    ALN_SEED="$2"; shift 2 ;;
    --aln-max-diff)       ALN_N="$2"; shift 2 ;;
    --aln-opens)          ALN_O="$2"; shift 2 ;;
    --min-mapq)           MIN_MAPQ="$2"; shift 2 ;;
    --dedup)              DEDUP=1; shift 1 ;;
    --damageprofiler-cmd) DP_CMD="$2"; shift 2 ;;
    --damageprofiler-jar) DP_JAR="$2"; shift 2 ;;
    --resume)             RESUME=1; shift 1 ;;
    -h|--help)            usage; exit 0 ;;
    --) shift; break ;;
    *) die "Unknown option: $1" ;;
  esac
done

# ---------------------------
# Validate
# ---------------------------
[ -n "$REF" ]       || { usage; die "Missing --reference-fasta"; }
[ -n "$READS_DIR" ] || { usage; die "Missing --reads-dir"; }
[ -n "$OUT_DIR" ]   || { usage; die "Missing --output-dir"; }
[ -f "$REF" ]       || die "Reference FASTA not found: $REF"
[ -d "$READS_DIR" ] || die "Reads directory not found: $READS_DIR"
[[ "$THREADS" =~ ^[0-9]+$ ]]  || die "--threads must be integer"
[[ "$ALN_SEED" =~ ^[0-9]+$ ]] || die "--aln-seed-length must be integer"
[[ "$ALN_O" =~ ^[0-9]+$ ]]    || die "--aln-opens must be integer"
# tools
command -v bwa >/dev/null 2>&1      || die "bwa not found in PATH"
command -v samtools >/dev/null 2>&1 || die "samtools not found in PATH"
if [ -n "$DP_JAR" ]; then
  command -v java >/dev/null 2>&1 || die "java not found in PATH (required for --damageprofiler-jar)"
else
  command -v "$DP_CMD" >/dev/null 2>&1 || die "damageprofiler CLI not found (set --damageprofiler-cmd or use --damageprofiler-jar)"
fi

# ---------------------------
# Prepare output structure
# ---------------------------
BAM_DIR="$OUT_DIR/bams"
DP_DIR="$OUT_DIR/damageprofiler"
LOG_DIR="$OUT_DIR/logs"
MAN_DIR="$OUT_DIR/manifests"
mkdir -p "$BAM_DIR" "$DP_DIR" "$LOG_DIR" "$MAN_DIR"

# Index reference (bwa + samtools)
if [ ! -f "${REF}.bwt" ]; then
  log "Indexing reference with bwa index"
  bwa index "$REF"
fi
[ -f "${REF}.fai" ] || samtools faidx "$REF"

# ---------------------------
# Iterate over R1 files and process pairs
# ---------------------------
manifest="$MAN_DIR/run_manifest.tsv"
printf "sample\tR1\tR2\tbam\taligned_reads\tmapq_filter\tdedup\tdamageprofiler_dir\n" > "$manifest"

shopt -s nullglob
for R1 in "$READS_DIR"/*"$R1_SUFFIX" "$READS_DIR"/*"${R1_SUFFIX%.gz}"; do
  [ -e "$R1" ] || continue
  sample=$(basename "$R1")
  sample=${sample%"$R1_SUFFIX"}
  sample=${sample%"${R1_SUFFIX%.gz}"}
  R2="${R1/$R1_SUFFIX/$R2_SUFFIX}"
  [ -f "$R2" ] || R2="${R1/${R1_SUFFIX%.gz}/${R2_SUFFIX%.gz}}"
  [ -f "$R2" ] || { log "WARNING: R2 not found for $R1 (skipping)"; continue; }

  BAM="$BAM_DIR/${sample}.bam"
  BAI="${BAM}.bai"
  DP_OUT="$DP_DIR/${sample}"

  if [ "$RESUME" -eq 1 ] && [ -f "$BAI" ] && [ -d "$DP_OUT" ]; then
    log "Resume: Skipping existing sample $sample"
    printf "%s\t%s\t%s\t%s\tNA\t%u\t%u\t%s\n" "$sample" "$R1" "$R2" "$BAM" "$MIN_MAPQ" "$DEDUP" "$DP_OUT" >> "$manifest"
    continue
  fi

  log "Aligning sample: $sample"

  if [ "$BWA_MODE" = "aln" ]; then
    # aDNA-friendly: disable seeding, relaxed mismatch (params configurable)
    SAI1="$BAM_DIR/${sample}.R1.sai"
    SAI2="$BAM_DIR/${sample}.R2.sai"
    bwa aln -l "$ALN_SEED" -n "$ALN_N" -o "$ALN_O" -t "$THREADS" "$REF" "$R1" > "$SAI1" 2> "$LOG_DIR/${sample}.aln.R1.log"
    bwa aln -l "$ALN_SEED" -n "$ALN_N" -o "$ALN_O" -t "$THREADS" "$REF" "$R2" > "$SAI2" 2> "$LOG_DIR/${sample}.aln.R2.log"
    bwa sampe "$REF" "$SAI1" "$SAI2" "$R1" "$R2" \
      2> "$LOG_DIR/${sample}.sampe.log" \
      | samtools view -b -q "$MIN_MAPQ" - \
      | samtools sort -@ "$THREADS" -o "$BAM" -
    rm -f "$SAI1" "$SAI2"
  else
    # mem mode (not aDNA-optimized by default)
    bwa mem -t "$THREADS" "$REF" "$R1" "$R2" 2> "$LOG_DIR/${sample}.mem.log" \
      | samtools view -b -q "$MIN_MAPQ" - \
      | samtools sort -@ "$THREADS" -o "$BAM" -
  fi

  if [ "$DEDUP" -eq 1 ]; then
    # Mark or remove duplicates (keep marked BAM for transparency)
    TMP_BAM="$BAM_DIR/${sample}.tmp.bam"
    samtools sort -n -@ "$THREADS" -o "$TMP_BAM" "$BAM"
    samtools fixmate -m "$TMP_BAM" "$BAM_DIR/${sample}.fixmate.bam"
    samtools sort -@ "$THREADS" -o "$BAM_DIR/${sample}.pos.bam" "$BAM_DIR/${sample}.fixmate.bam"
    samtools markdup -@ "$THREADS" "$BAM_DIR/${sample}.pos.bam" "$BAM"
    rm -f "$TMP_BAM" "$BAM_DIR/${sample}.fixmate.bam" "$BAM_DIR/${sample}.pos.bam"
  fi

  samtools index "$BAM"

  # Count aligned reads (primary alignments)
  aligned=$(samtools view -c -F 0x900 "$BAM" || echo "NA")

  # DamageProfiler
  mkdir -p "$DP_OUT"
  log "Running DamageProfiler for $sample"
  if [ -n "$DP_JAR" ]; then
    # Common jar usage pattern:
    #   java -jar DamageProfiler.jar -i <bam> -r <ref> -o <outdir> -t <threads>
    java -jar "$DP_JAR" -i "$BAM" -r "$REF" -o "$DP_OUT" -t "$THREADS" > "$LOG_DIR/${sample}.damageprofiler.log" 2>&1
  else
    # CLI wrapper usage (may vary slightly by version)
    "$DP_CMD" -i "$BAM" -r "$REF" -o "$DP_OUT" -t "$THREADS" > "$LOG_DIR/${sample}.damageprofiler.log" 2>&1
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%u\t%u\t%s\n" "$sample" "$R1" "$R2" "$BAM" "$aligned" "$MIN_MAPQ" "$DEDUP" "$DP_OUT" >> "$manifest"
done
shopt -u nullglob

log "Done. BAMs in $BAM_DIR, damage profiles in $DP_DIR"
log "Manifest: $manifest"
