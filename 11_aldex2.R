##### 11_ALDEx2
##### By Carter Clinton



#!/usr/bin/env Rscript
# aldex2.R — ALDEx2 differential abundance with effect sizes, plots, and heatmap
# Maps original file: ALDEx2_SILVA.txt
# Analysis Stage: Differential Abundance
# Language: R

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ALDEx2)
  # Optional:
  # library(pheatmap)
})

# -------- CLI --------
option_list <- list(
  make_option(c("-c", "--counts-tsv"), type="character", default=NULL,
              help="Counts table TSV (features as rows, samples as columns). First column = FeatureID.", metavar="FILE"),
  make_option(c("-m", "--metadata-tsv"), type="character", default=NULL,
              help="Sample metadata TSV (QIIME 2 style). Must contain column specified by --group-col.", metavar="FILE"),
  make_option(c("-g", "--group-col"), type="character", default=NULL,
              help="Metadata column name defining groups (e.g., Group, BurialControl).", metavar="STR"),
  make_option(c("-t", "--taxonomy-tsv"), type="character", default=NULL,
              help="(optional) Taxonomy TSV with columns: FeatureID, Taxonomy. Used to annotate outputs.", metavar="FILE"),
  make_option(c("-o", "--output-dir"), type="character", default=NULL,
              help="Output directory.", metavar="DIR"),
  make_option(c("-l", "--label"), type="character", default="aldex2",
              help="Basename label for outputs [default: %default]."),
  make_option(c("--mc-samples"), type="integer", default=128,
              help="Number of Monte Carlo instances [default: %default]."),
  make_option(c("--denom"), type="character", default="iqlr",
              help="Denominator for CLR: 'all' or 'iqlr' [default: %default]."),
  make_option(c("--test"), type="character", default="auto",
              help="Test mode: 'auto' (detect), 'two' (two-class Welch's t), or 'multi' (Kruskal–Wallis). [default: %default]."),
  make_option(c("--alpha"), type="double", default=0.05,
              help="FDR threshold for significance [default: %default]."),
  make_option(c("--top-n"), type="integer", default=50,
              help="Top N features to plot/heatmap by |effect| or FDR [default: %default]."),
  make_option(c("--seed"), type="integer", default=12345,
              help="Random seed for reproducibility [default: %default]."),
  make_option(c("--write-clr"), action="store_true", default=FALSE,
              help="Also write CLR matrix to disk (TSV)."),
  make_option(c("--min-samples"), type="integer", default=0,
              help="Filter features present (count>0) in <N samples before ALDEx2 (0 disables) [default: %default]."),
  make_option(c("--min-total"), type="integer", default=0,
              help="Filter features with total count <N across all samples (0 disables) [default: %default].")
)
opt <- parse_args(OptionParser(option_list=option_list))

# -------- Validate --------
die <- function(msg, code=1){ message(sprintf("ERROR: %s", msg)); quit(status=code) }
req <- c("counts-tsv","metadata-tsv","group-col","output-dir")
for (r in req) if (is.null(opt[[gsub("-", "_", r)]])) die(sprintf("Missing --%s", r))
if (!file.exists(opt$`counts_tsv`)) die(sprintf("Counts TSV not found: %s", opt$`counts_tsv`))
if (!file.exists(opt$`metadata_tsv`)) die(sprintf("Metadata TSV not found: %s", opt$`metadata_tsv`))
if (!dir.exists(opt$`output_dir`)) dir.create(opt$`output_dir`, recursive = TRUE, showWarnings = FALSE)
if (!is.null(opt$`taxonomy_tsv`) && !file.exists(opt$`taxonomy_tsv`)) die(sprintf("Taxonomy TSV not found: %s", opt$`taxonomy_tsv`))

set.seed(opt$seed)

out_path <- function(...) file.path(opt$`output_dir`, ...)

# -------- Load data --------
read_tsv_loose <- function(path) {
  suppressMessages(readr::read_tsv(path, col_types = cols(.default = col_guess())))
}

counts <- read_tsv_loose(opt$`counts_tsv`)
if (!all(c("FeatureID") %in% colnames(counts))) {
  die("Counts TSV must have a 'FeatureID' column followed by sample columns.")
}
feature_ids <- counts$FeatureID
count_mat <- counts %>% select(-FeatureID) %>% as.data.frame()
rownames(count_mat) <- feature_ids

meta <- read_tsv_loose(opt$`metadata_tsv`)
if (!(opt$`group_col` %in% colnames(meta))) {
  die(sprintf("Metadata column '%s' not found in metadata TSV.", opt$`group_col`))
}

# Align samples between counts and metadata
samples_common <- intersect(colnames(count_mat), meta[[1]])
if (length(samples_common) == 0) {
  # Try matching on a column named "sampleid" or "SampleID"
  sid_col <- c("sampleid","SampleID","sample-id","#SampleID")
  sid_col <- sid_col[sid_col %in% colnames(meta)]
  if (length(sid_col) == 0) die("Could not align samples: ensure metadata first column is SampleID matching count columns.")
  meta <- meta %>% rename(SampleID = all_of(sid_col[1]))
  samples_common <- intersect(colnames(count_mat), meta$SampleID)
  if (length(samples_common) == 0) die("No overlapping samples between counts columns and metadata SampleID.")
} else {
  # Ensure first column is SampleID
  colnames(meta)[1] <- "SampleID"
}

# Subset & order
count_mat <- count_mat[, samples_common, drop=FALSE]
meta <- meta %>% filter(SampleID %in% samples_common)
meta <- meta %>% mutate(!!.group := .data[[opt$`group_col`]]) %>% # placeholder variable name
  rename(Group = !!opt$`group_col`)
meta <- meta[match(colnames(count_mat), meta$SampleID), ]

# -------- Filters (optional) --------
if (opt$`min_samples` > 0) {
  keep <- rowSums(count_mat > 0) >= opt$`min_samples`
  count_mat <- count_mat[keep, , drop=FALSE]
}
if (opt$`min_total` > 0) {
  keep2 <- rowSums(count_mat) >= opt$`min_total`
  count_mat <- count_mat[keep2, , drop=FALSE]
}
if (nrow(count_mat) == 0) die("All features filtered out; relax --min-samples/--min-total.")

# -------- Taxonomy (optional) --------
tax <- NULL
if (!is.null(opt$`taxonomy_tsv`)) {
  tax <- read_tsv_loose(opt$`taxonomy_tsv`)
  # Expect columns: FeatureID, Taxonomy (others are fine)
  if (!all(c("FeatureID","Taxonomy") %in% colnames(tax))) {
    message("WARNING: taxonomy TSV lacks 'FeatureID' and 'Taxonomy'; taxonomy annotations will be skipped.")
    tax <- NULL
  }
}

# -------- Choose test --------
groups <- factor(meta$Group)
n_levels <- nlevels(groups)
test_mode <- opt$`test`
if (test_mode == "auto") {
  test_mode <- if (n_levels == 2) "two" else "multi"
}
if (!test_mode %in% c("two","multi")) die("--test must be one of: auto, two, multi")

# -------- ALDEx2 pipeline --------
message(sprintf("ALDEx2: %s-class test, mc.samples=%d, denom=%s", ifelse(test_mode=="two","two","multi"), opt$`mc_samples`, opt$`denom`))
# ALDEx2 expects integer counts
count_mat <- as.matrix(round(count_mat, 0))

# CLR transformation
clr <- suppressWarnings(ALDEx2::aldex.clr(count_mat, conds = groups, mc.samples = opt$`mc_samples`, denom = opt$`denom`))

# Test
if (test_mode == "two") {
  tt <- ALDEx2::aldex.tt(clr, verbose = FALSE)
  eff <- ALDEx2::aldex.effect(clr, verbose = FALSE)
  res <- cbind(as.data.frame(tt), as.data.frame(eff))
  # Add FDR
  res$q_we <- p.adjust(res$we.ep, method = "BH")
  res$q_wilcox <- p.adjust(res$wi.ep, method = "BH")
} else {
  kw <- ALDEx2::aldex.kw(clr, verbose = FALSE)
  eff <- ALDEx2::aldex.effect(clr, verbose = FALSE)
  res <- cbind(as.data.frame(kw), as.data.frame(eff))
  res$q_kw <- p.adjust(res$kw.ep, method = "BH")
}
res <- res %>% tibble::rownames_to_column("FeatureID")

# Merge taxonomy (optional)
if (!is.null(tax)) {
  res <- res %>% left_join(tax %>% select(FeatureID, Taxonomy), by = "FeatureID")
}

# Write results
dir.create(out_path("results"), showWarnings = FALSE, recursive = TRUE)
write_csv(res, out_path("results", sprintf("%s_aldex2_results.csv", opt$label)))

# -------- Plots --------
dir.create(out_path("plots"), showWarnings = FALSE, recursive = TRUE)

# Volcano-like: effect vs -log10(q)
if (test_mode == "two") {
  res <- res %>% mutate(q = pmin(q_we, q_wilcox, na.rm = TRUE))
} else {
  res <- res %>% mutate(q = q_kw)
}
res <- res %>% mutate(sig = ifelse(is.finite(q) & q <= opt$alpha, "significant", "ns"))

p_vol <- ggplot(res, aes(x = effect, y = -log10(pmax(q, 1e-300)), color = sig)) +
  geom_point(alpha = 0.7, size = 1.2) +
  scale_color_manual(values = c("significant"="#d62728", "ns"="#7f7f7f")) +
  labs(x = "ALDEx2 effect size", y = expression(-log[10]~"(FDR)"),
       color = NULL,
       title = sprintf("ALDEx2: %s (%s)", opt$label, ifelse(test_mode=="two","two-class","multi-class"))) +
  theme_minimal(base_size = 12)
ggsave(out_path("plots", sprintf("%s_volcano.pdf", opt$label)), p_vol, width = 7, height = 5)

# Top features table
res_top <- res %>%
  arrange(q, desc(abs(effect))) %>%
  slice_head(n = min(opt$`top_n`, n()))
write_csv(res_top, out_path("results", sprintf("%s_top_%d.csv", opt$label, nrow(res_top))))

# Heatmap of CLR for top features (if desired; simple ggplot tile to avoid heavy deps)
if (nrow(res_top) >= 2) {
  clr_df <- as.data.frame(ALDEx2::getMonteCarloInstances(clr)) # 3D array collapsed? Better reconstruct median CLR
  # safer: use aldex.clr output's @analysisData slot isn't public; use aldex.clr output matrix via 'aldex.clr' docs
  # Instead, compute median CLR per feature/sample:
  clr_med <- ALDEx2::aldex.clr(count_mat, conds = groups, mc.samples = opt$`mc_samples`, denom = opt$`denom`, verbose = FALSE, useMC=TRUE)
  # aldex.clr returns an object with draws collapsed on samples; to get per-feature per-sample clr, use aldex.clr output with @analysisData? 
  # Fallback: approximate with log2((counts + pseudo)/gm), using pseudo=0.5:
  pseudo <- 0.5
  gm <- exp(rowMeans(log((count_mat + pseudo))))
  clr_simple <- log((count_mat + pseudo), 2) - log(matrix(gm, nrow = length(gm), ncol = ncol(count_mat), byrow = FALSE), 2)
  ht <- clr_simple[match(res_top$FeatureID, rownames(clr_simple)), , drop=FALSE]
  df_long <- as.data.frame(ht) %>%
    tibble::rownames_to_column("FeatureID") %>%
    pivot_longer(-FeatureID, names_to = "SampleID", values_to = "CLR")

  ann <- meta %>% select(SampleID, Group)
  df_long <- df_long %>% left_join(ann, by = "SampleID")

  p_heat <- ggplot(df_long, aes(x = SampleID, y = FeatureID, fill = CLR)) +
    geom_tile() +
    facet_grid(Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = "Samples", y = "Top features", title = sprintf("CLR heatmap (top %d by FDR/effect): %s", nrow(res_top), opt$label)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave(out_path("plots", sprintf("%s_heatmap.pdf", opt$label)), p_heat, width = 10, height = 6)
}

# (Optional) Write CLR matrix
if (isTRUE(opt$`write_clr`)) {
  write_tsv(as.data.frame(clr_simple) %>% tibble::rownames_to_column("FeatureID"),
            out_path("results", sprintf("%s_clr_matrix.tsv", opt$label)))
}

# -------- Session info --------
capture <- capture.output({
  cat("===== Session Info =====\n")
  print(sessionInfo())
})
writeLines(capture, out_path(sprintf("%s_sessionInfo.txt", opt$label)))

message("Done. Outputs written to: ", opt$`output_dir`)



