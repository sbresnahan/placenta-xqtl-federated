#!/usr/bin/env Rscript
# Ported verbatim from Placenta_xQTL_atlas/05_qtl_mapping/compute_qvalues.R (commit 22be872).
# Only change: the hard-coded seadragon lab library is replaced by the optional
# XQTL_R_LIBS env var (colon-separated), prepended to .libPaths() when set.
local({ x <- Sys.getenv("XQTL_R_LIBS"); if (nzchar(x)) .libPaths(c(strsplit(x, ":")[[1]], .libPaths())) })
# =============================================================================
# compute_qvalues.R — Storey q-values via the R qvalue package (GTEx convention)
# =============================================================================
# File-based bridge called by 27_run_tensorqtl.py. This replaces tensorQTL's
# calculate_qvalues() rpy2 path, which is fragile on seadragon (conda-R /
# libstdc++ conflicts). Statistically identical: tensorQTL's calculate_qvalues
# is itself a wrapper around qvalue::qvalue with default lambda.
#
# Usage: Rscript compute_qvalues.R <in.tsv> <out.tsv>
#   in.tsv : single column 'pval_beta' (permutation-calibrated p-values)
#   out.tsv: same rows + 'qval' column. Prints pi0 to stdout.
# Exit status 1 (and no output file) on any failure — the caller hard-errors.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript compute_qvalues.R <in.tsv> <out.tsv>")
}
in_file <- args[1]
out_file <- args[2]

ok <- tryCatch({
  suppressPackageStartupMessages(library(qvalue))

  df <- read.delim(in_file)
  if (!("pval_beta" %in% colnames(df))) {
    stop("input TSV must have a 'pval_beta' column")
  }
  p <- as.numeric(df$pval_beta)
  p[!is.finite(p)] <- 1          # NA/NaN -> 1 (conservative; matches BH fallback)
  p <- pmin(pmax(p, 0), 1)       # clamp to [0,1]
  if (length(unique(p)) < 2) {
    stop("fewer than 2 unique p-values — qvalue cannot estimate pi0")
  }

  # qvalue <= 2.38 pi0est() bug workaround: the default lambda grid
  # (0.05..0.95 by 0.05) assumes at least one p-value >= max(lambda).
  # pi0est computes pi0(lambda) via tabulate(findInterval(p, lambda)) WITHOUT
  # an nbins argument, so when max(p) <= 0.95 the tabulated vector is short,
  # the reverse-indexing yields NA, and smooth.spline aborts with
  # "missing or infinite values in inputs are not allowed". This bites
  # strong-signal subsets (e.g. chr1 HCP-optimization arms) where every
  # group's pval_beta can fall below 0.95. Remedy — the one qvalue's own
  # (devel) error message recommends: restrict lambda to the observed range.
  # Falls back to a fixed-lambda Storey pi0 when fewer than 4 grid points
  # remain (smoother minimum), and finally to pi0 = 1 (BH-equivalent,
  # conservative) when max(p) < 0.05.
  lam <- seq(0.05, 0.95, 0.05)
  if (max(p) <= max(lam)) {
    lam <- lam[lam < max(p)]
    cat(sprintf(paste0("  NOTE: max(pval_beta) = %.4f <= 0.95 — restricting ",
                       "qvalue lambda grid to %d point(s)\n"),
                max(p), length(lam)))
  }
  # Nested fallback: if the pi0 smoother still fails (e.g. "estimated pi0
  # <= 0" on heavily discrete inputs), retry with lambda = 0, which forces
  # pi0 = 1 — BH-equivalent q-values, conservative, and cannot fail (the
  # single-lambda pi0est path uses no spline). A pathological k then simply
  # deflates its own eGene count instead of aborting the whole grid.
  q <- tryCatch({
    if (length(lam) >= 4) {
      qvalue(p, lambda = lam)      # smoother pi0 (GTEx convention)
    } else if (length(lam) >= 1) {
      qvalue(p, lambda = max(lam)) # fixed-lambda Storey pi0 (no spline)
    } else {
      qvalue(p, lambda = 0)        # pi0 = 1: conservative, BH-equivalent
    }
  }, error = function(e) {
    cat(sprintf(paste0("  NOTE: qvalue smoother failed (%s); falling back ",
                       "to pi0 = 1 (BH-equivalent, conservative)\n"),
                conditionMessage(e)))
    qvalue(p, lambda = 0)
  })
  df$qval <- q$qvalues
  cat(sprintf("  R qvalue: pi0 = %.4f (estimated signal fraction 1-pi0 = %.4f)\n",
              q$pi0, 1 - q$pi0))
  write.table(df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  TRUE
}, error = function(e) {
  cat(sprintf("ERROR in compute_qvalues.R: %s\n", conditionMessage(e)), file = stderr())
  FALSE
})

if (!ok) quit(save = "no", status = 1)
