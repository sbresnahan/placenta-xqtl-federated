#!/usr/bin/env Rscript
# xqtl_normalize.R -- normalize one phenotype modality within one ancestry stratum.
#
#   xqtl_normalize.R --bed expression.bed.gz [--counts expression.counts.bed.gz] \
#       --modality expression --samples EUR.samples.tsv [--exclude EUR_expression.outliers.tsv] \
#       --out EUR_expression
#
# Ported from the atlas (05_qtl_mapping/normalize_expression_hcp.R, vst mode, and
# normalize_modalities.R). At a single site the stratum is one cohort, so "within-cohort INT"
# is INT over the stratum.
#   expression          TPM > 0.1 in > 25% of samples, drop no-variance genes; TMM -> VST on
#                       counts (terminal, no INT); connectivity outliers (signed bicor adjacency,
#                       z < -3) are removed and written to OUT.outliers.tsv
#   isoform_usage,      TPM-style filter (> 0.1 in > 25%); isoform_expression uses TMM -> VST
#   isoform_expression  when --counts is given
#   all others          detected (non-NA) in >= 40% of samples
#   then (non-VST)      drop no-variance / < 2 observed features; median impute; quantile
#                       normalization (ties averaged); rank-based inverse normal transform
# Technical replicates: columns are renamed rnaseq_id -> genotype_id and replicate columns are
# averaged after normalization (atlas order). Without --counts, expression falls back to
# QN + INT (reported in OUT.norm_summary.tsv).
#
# Outputs: OUT.bed.gz (columns = genotype IDs), OUT.phenotype_groups.txt (if phenotype IDs are
# GENE__feature), OUT.features.tsv (every input feature: coordinates, group, detection fraction,
# filter status), OUT.norm_summary.tsv, OUT.outliers.tsv (connectivity outliers; expression only, header-only otherwise).
suppressPackageStartupMessages({ library(optparse); library(data.table) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--bed"), make_option("--counts", default = NA_character_),
  make_option("--modality"), make_option("--samples"),
  make_option("--exclude", default = NA_character_),
  make_option("--tpm-min", dest = "tpm_min", type = "double", default = 0.1),
  make_option("--tpm-min-prop", dest = "tpm_min_prop", type = "double", default = 0.25),
  make_option("--min-detect-prop", dest = "min_detect_prop", type = "double", default = 0.4),
  make_option("--outlier-z", dest = "outlier_z", type = "double", default = -3),
  make_option("--min-samples", dest = "min_samples", type = "integer", default = 10L),
  make_option("--out")
)))
for (k in c("bed", "modality", "samples", "out")) if (is.null(opt[[k]])) stop("xqtl_normalize.R: --", k, " is required")
mod <- opt$modality
TPM_MODALITIES <- c("expression", "isoform_usage", "isoform_expression")

# ---- helpers (atlas) --------------------------------------------------------------------
quantile_normalize_rows <- function(mat) {           # samples x features
  t_mat <- t(mat)
  mean_dist <- rowMeans(apply(t_mat, 2, sort))
  result <- matrix(0, nrow(t_mat), ncol(t_mat))
  for (j in seq_len(ncol(t_mat))) {
    o <- order(t_mat[, j]); result[o, j] <- mean_dist
    runs <- rle(t_mat[o, j])
    if (any(runs$lengths > 1)) {
      ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1
      for (g in which(runs$lengths > 1)) result[o[starts[g]:ends[g]], j] <- mean(mean_dist[starts[g]:ends[g]])
    }
  }
  dimnames(result) <- dimnames(t_mat); t(result)
}
inverse_normal_transform <- function(mat) {
  r <- apply(mat, 2, function(x) qnorm(rank(x, ties.method = "average", na.last = "keep") / (sum(!is.na(x)) + 1)))
  r <- matrix(r, nrow(mat), ncol(mat), dimnames = dimnames(mat)); r
}
tmm_vst <- function(counts_mat) {                    # samples x features
  g <- t(counts_mat); lib <- colSums(g)
  if (any(!is.finite(lib) | lib <= 0)) stop("tmm_vst: sample(s) with zero library size: ",
                                            paste(head(names(lib)[!is.finite(lib) | lib <= 0], 10), collapse = ", "))
  eff <- lib * edgeR::calcNormFactors(g, method = "TMM")
  dds <- DESeq2::DESeqDataSetFromMatrix(round(g), colData = data.frame(row.names = colnames(g)), design = ~ 1)
  DESeq2::sizeFactors(dds) <- eff / exp(mean(log(eff)))
  t(SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE, nsub = min(1000L, nrow(g)))))
}
bicor_samples <- function(X) {
  med <- apply(X, 1, median, na.rm = TRUE)
  mad1 <- apply(X, 1, function(x) mad(x, constant = 1, na.rm = TRUE)); mad1[!is.finite(mad1) | mad1 == 0] <- NA
  U <- (X - med) / (9 * mad1); W <- (1 - U^2)^2; W[!is.finite(W) | abs(U) >= 1] <- 0
  Y <- (X - med) * W; Y[!is.finite(Y)] <- 0
  d2 <- rowSums(Y^2); bc <- tcrossprod(Y) / sqrt(outer(d2, d2)); bc[!is.finite(bc)] <- 0
  dimnames(bc) <- list(rownames(X), rownames(X)); bc
}
connectivity_outliers <- function(mat, z) {
  adj <- (1 + bicor_samples(mat)) / 2; diag(adj) <- 0
  conn <- rowSums(adj, na.rm = TRUE); s <- sd(conn)
  if (!is.finite(s) || s == 0) return(character(0))
  rownames(mat)[(conn - mean(conn)) / s < z]
}
read_bed <- function(path) {
  b <- if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", header = TRUE, check.names = FALSE)
       else fread(path, sep = "\t", header = TRUE, check.names = FALSE)
  setnames(b, 1:4, c("chr", "start", "end", "phenotype_id")); b
}

# ---- inputs -------------------------------------------------------------------------------
smp <- fread(opt$samples, colClasses = "character")
excl <- character(0)
if (!is.na(opt$exclude) && file.exists(opt$exclude)) excl <- fread(opt$exclude, colClasses = "character")[[1]]
bed <- read_bed(opt$bed)
cols <- setdiff(intersect(smp$rnaseq_id, names(bed)), excl)
missing <- setdiff(smp$rnaseq_id, c(names(bed), excl))
if (length(missing)) message(sprintf("%s: %d stratum sample(s) absent from the %s BED", opt$out, length(missing), mod))
if (length(cols) < opt$min_samples) stop(sprintf("xqtl_normalize.R: only %d samples for %s (need %d)", length(cols), mod, opt$min_samples))
X <- t(as.matrix(bed[, ..cols])); colnames(X) <- bed$phenotype_id
n_in <- ncol(X)

no_var <- apply(X, 2, function(x) length(unique(x[!is.na(x)])) <= 1)
few_obs <- colSums(!is.na(X)) < 2
if (mod %in% TPM_MODALITIES) {
  low <- colMeans(X > opt$tpm_min, na.rm = TRUE) <= opt$tpm_min_prop
} else {
  low <- colMeans(!is.na(X)) < opt$min_detect_prop
}
low[is.na(low)] <- TRUE
detect <- if (mod %in% TPM_MODALITIES) colMeans(X > opt$tpm_min, na.rm = TRUE) else colMeans(!is.na(X))
feat <- data.table(phenotype_id = colnames(X), detection_fraction = round(detect, 4),
                   status = fifelse(low, "low_detection", fifelse(no_var, "no_variance", fifelse(few_obs, "too_few_observed", "kept"))))
X <- X[, !(low | no_var | few_obs), drop = FALSE]
if (ncol(X) == 0) stop("xqtl_normalize.R: no ", mod, " features pass the filters")

use_vst <- !is.na(opt$counts) && mod %in% c("expression", "isoform_expression")
if (!is.na(opt$counts) && !use_vst) stop("--counts is only used for expression and isoform_expression")
if (use_vst) {
  cb <- read_bed(opt$counts)
  gi <- match(colnames(X), cb$phenotype_id)
  if (anyNA(gi)) stop("counts BED lacks ", sum(is.na(gi)), " of the filtered features")
  if (length(setdiff(cols, names(cb)))) stop("counts BED lacks some samples")
  C <- t(as.matrix(cb[gi, ..cols])); colnames(C) <- colnames(X)
  cv <- apply(C, 2, var); okv <- is.finite(cv) & cv > 0
  C <- C[, okv, drop = FALSE]
  N <- tmm_vst(C); method <- "TMM-VST"
} else {
  for (j in seq_len(ncol(X))) { na <- is.na(X[, j]); if (any(na)) X[na, j] <- median(X[!na, j]) }
  N <- inverse_normal_transform(quantile_normalize_rows(X)); method <- "QN-INT"
}
rownames(N) <- cols

outliers <- character(0)
if (mod == "expression") {
  outliers <- connectivity_outliers(N, opt$outlier_z)
  if (length(outliers)) N <- N[setdiff(rownames(N), outliers), , drop = FALSE]
}
# always written (header only when no outliers or modality != expression)
fwrite(data.table(rnaseq_id = outliers), paste0(opt$out, ".outliers.tsv"), sep = "\t")

# ---- replicate averaging (rnaseq_id -> genotype_id) ---------------------------------------
gid <- setNames(smp$genotype_id, smp$rnaseq_id)[rownames(N)]
if (anyNA(gid)) stop("samples without genotype_id: ", paste(rownames(N)[is.na(gid)], collapse = ", "))
n_rep <- sum(duplicated(gid))
if (n_rep) {
  N <- rowsum(N, gid, reorder = TRUE) / as.vector(table(gid)[sort(unique(gid))])
} else {
  rownames(N) <- gid; N <- N[order(rownames(N)), , drop = FALSE]
}

meta <- bed[match(colnames(N), bed$phenotype_id), .(chr, start, end, phenotype_id)]
out <- cbind(meta, as.data.table(t(N)))
setorder(out, chr, start, end, phenotype_id)
setnames(out, "chr", "#chr")
fwrite(out, paste0(opt$out, ".bed.gz"), sep = "\t", quote = FALSE, compress = "gzip")
if (any(grepl("__", out$phenotype_id, fixed = TRUE))) {
  fwrite(data.table(out$phenotype_id, sub("__.*$", "", out$phenotype_id)),
         paste0(opt$out, ".phenotype_groups.txt"), sep = "\t", col.names = FALSE)
}
feat[status == "kept" & !(phenotype_id %in% out$phenotype_id), status := "zero_count_variance"]
fm <- bed[, .(phenotype_id, chr, start, end)]
feat <- merge(fm, feat, by = "phenotype_id", sort = FALSE)
feat[, group_id := sub("__.*$", "", phenotype_id)]
if ("novel" %in% names(bed)) feat[, novel := bed$novel[match(phenotype_id, bed$phenotype_id)]]
fwrite(feat, paste0(opt$out, ".features.tsv"), sep = "\t", quote = FALSE)
fwrite(data.table(modality = mod, method = method, features_in = n_in, features_kept = nrow(out),
                  samples_in = length(cols), connectivity_outliers = length(outliers),
                  replicates_averaged = n_rep, individuals = nrow(N)),
       paste0(opt$out, ".norm_summary.tsv"), sep = "\t")
message(sprintf("%s: %s %d -> %d features, %d individuals (%s; %d outlier(s), %d replicate column(s) averaged)",
                opt$out, mod, n_in, nrow(out), nrow(N), method, length(outliers), n_rep))
