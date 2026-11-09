#!/usr/bin/env Rscript
# Optional XQTL_R_LIBS (colon-separated) is prepended to the R library path.
# =============================================================================
# qu_correct_salmon.R
# -----------------------------------------------------------------------------
# Apply quantification-uncertainty (QU) correction to Salmon transcript
# quantifications, per:
#   Chen, Y., Lun, A.T.L., Zhou, X., et al. (2023).
#   "Dividing out quantification uncertainty allows efficient assessment of
#    differential transcript expression with edgeR."
#   Nucleic Acids Research, gkad1167.  https://doi.org/10.1093/nar/gkad1167
#
# Method
# -------
# Salmon (and kallisto) emit bootstrap resamples of transcript counts. The
# read-to-transcript ambiguity (RTA) in those bootstraps induces extra-Poisson
# variation (an "overdispersion" sigma^2_t per transcript). The edgeR function
# catchSalmon fits a quasi-Poisson model to the bootstraps and estimates
# sigma^2_t. Dividing each transcript's count by its overdispersion
#   z_ti = y_ti / sigma^2_t
# yields "scaled counts" that follow the standard negative-binomial
# mean-variance relationship, so downstream gene-level tools (edgeR, etc.) work
# at full statistical efficiency.
#
# IMPORTANT: catchSalmon is in the **edgeR** Bioconductor package (NOT fishpond;
# fishpond provides Swish, a different method). The paper states this
# explicitly.
#
# What this script does
# ---------------------
# 1. Reads the per-sample Salmon directories (<salmon-dir>/<sample>), in the
#    order given by <samples.txt>.
# 2. edgeR::catchSalmon(dirs) -> per-transcript RTA overdispersion.
# 3. Builds tx2gene (transcript_id -> gene_id) from the reference GTF via
#    rtracklayer::import (standard, fast, C-backed GTF parser).
# 4. tximport(..., type="salmon", tx2gene=tx2g, txOut=TRUE, dropInfReps=TRUE)
#    to import the quant.sf files (inferential replicates dropped AFTER
#    catchSalmon has consumed the bootstraps).
# 5. QU correction: counts <- counts / overdispersion (clamped to >= 1, as the
#    paper's estimator does via max(1, ...)).
# 6. Recomputes abundance as CPM from the scaled counts (edgeR::cpm), matching
#    the reference workflow.
# 7. Writes an ADJUSTED quant.sf per sample to <out-dir>/<sample>/quant.sf,
#    with the same 5 columns Salmon uses (Name, Length, EffectiveLength, TPM,
#    NumReads):
#       Name            = transcript_id
#       NumReads        = QU-corrected (scaled) counts
#       TPM             = CPM recomputed from scaled counts (per-sample
#                         scaled_count / sum(scaled_count) * 1e6). NOTE: this
#                         is count-proportional CPM, NOT length-normalized
#                         TPM — the column keeps the Salmon name only for
#                         drop-in quant.sf compatibility. Downstream,
#                         assemble_bed.py is called with --units
#                         tpm_from_counts, which recomputes true
#                         length-normalized TPM from NumReads x
#                         EffectiveLength, so the CPM-scale column is not
#                         used for the isoform BEDs.
#       Length          = carried from the original quant.sf
#       EffectiveLength = carried from the original quant.sf
#    This makes <out-dir> a drop-in replacement for the original Salmon dir
#    from assemble_bed.py's perspective (it reads Name/TPM/NumReads, plus
#    EffectiveLength for the tpm_from_counts mode).
#
# Scope
# -----
# QU correction is applied to ISOFORM (transcript) counts only. Gene-level
# expression is produced by assemble_bed.py from the ORIGINAL Salmon dir
# (gene-level counts are largely unambiguous, so RTA is negligible there).
#
# Usage
# -----
#   Rscript qu_correct_salmon.R \
#       --samples    <samples.txt> \
#       --salmon-dir <cohort_dir>/intermediate/expression \
#       --ref-anno   <HPLRv2.0.annotated.PANTRY.gtf> \
#       --out-dir    <cohort_dir>/intermediate/expression_qu
#
# Requires in the seadragon R library: edgeR, tximport, rtracklayer (and their
# deps). Errors clearly if a package is missing.
# =============================================================================

local({ x <- Sys.getenv("XQTL_R_LIBS"); if (nzchar(x)) .libPaths(c(strsplit(x, ":")[[1]], .libPaths())) })

suppressPackageStartupMessages(library(optparse))

# =============================================================================
# Helper: build tx2gene (transcript_id -> gene_id) from a GTF via rtracklayer
# =============================================================================
# rtracklayer::import is the standard, fast (C-backed) GTF parser. We import
# the GTF, subset to transcript features, extract gene_id + transcript_id, and
# deduplicate. This replaces an earlier base-R line-by-line parser that was
# slower and more fragile; rtracklayer is near-universal in any Bioconductor
# environment that already has edgeR + tximport.
.build_tx2gene <- function(gtf_file) {
  if (!file.exists(gtf_file)) stop("GTF not found: ", gtf_file)
  gr <- rtracklayer::import(gtf_file)            # GRanges; mcols hold attributes
  gr <- gr[gr$type == "transcript"]              # subset to transcript features
  if (length(gr) == 0) {
    stop("No transcript features with type=='transcript' found in GTF: ", gtf_file)
  }
  tx2g <- data.frame(
    transcript_id = gr$transcript_id,
    gene_id       = gr$gene_id,
    stringsAsFactors = FALSE
  )
  tx2g <- unique(tx2g)                           # a transcript may span lines
  # Drop any rows missing either ID (shouldn't happen for a well-formed GTF):
  tx2g <- tx2g[!is.na(tx2g$transcript_id) & !is.na(tx2g$gene_id) &
               nzchar(tx2g$transcript_id) & nzchar(tx2g$gene_id), ]
  if (nrow(tx2g) == 0) {
    stop("No transcript_id/gene_id attributes on transcript features in GTF: ", gtf_file)
  }
  tx2g
}

option_list <- list(
  make_option("--samples", type = "character", default = NULL,
              help = "File with one sample ID per line (no header)"),
  make_option("--salmon-dir", type = "character", default = NULL,
              help = "Directory containing per-sample Salmon output subdirs (each with quant.sf + aux_info/)"),
  make_option("--ref-anno", type = "character", default = NULL,
              help = "Reference GTF with transcript_id / gene_id attributes (used to build tx2gene)"),
  make_option("--out-dir", type = "character", default = NULL,
              help = "Output directory: adjusted quant.sf written to <out-dir>/<sample>/quant.sf"),
  make_option("--chunk-size", type = "integer", default = 10,
              help = "Process samples in chunks of this size (0 = all at once). Reduces peak memory ~proportionally. [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("samples", "salmon-dir", "ref-anno", "out-dir")
missing <- required[vapply(required, function(k) is.null(opt[[k]]), logical(1))]
if (length(missing) > 0) {
  stop(sprintf("Missing required argument(s): %s. Required: --samples, --salmon-dir, --ref-anno, --out-dir",
               paste0("--", missing, collapse = ", ")))
}

samples_file   <- opt[["samples"]]
salmon_dir     <- opt[["salmon-dir"]]
ref_anno_file  <- opt[["ref-anno"]]
out_dir        <- opt[["out-dir"]]
chunk_size     <- opt[["chunk-size"]]

# ---- Load required packages ------------------------------------------------
for (pkg in c("edgeR", "tximport", "rtracklayer")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required R package '%s' not found in .libPaths(): %s\n",
                 pkg, paste(.libPaths(), collapse = ", ")))
  }
}
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(rtracklayer))

# ---- Read samples ----------------------------------------------------------
samples <- readLines(samples_file)
samples <- trimws(samples)
samples <- samples[nzchar(samples)]
if (length(samples) == 0) stop("No samples found in samples file: ", samples_file)

dirs <- file.path(salmon_dir, samples)
missing_dirs <- dirs[!dir.exists(dirs)]
if (length(missing_dirs) > 0) {
  stop("Missing Salmon output directories:\n  ", paste(missing_dirs, collapse = "\n  "))
}
files <- file.path(dirs, "quant.sf")
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing quant.sf files:\n  ", paste(missing_files, collapse = "\n  "))
}
names(files) <- samples

# Pre-flight: catchSalmon needs each sample's aux_info/ (meta_info.json +
# bootstrap resamples). A missing aux_info surfaces deep inside catchSalmon
# as a cryptic jsonlite "lexical error: invalid char in json text" (fromJSON
# parses the nonexistent path as literal JSON text) — fail fast with a clear
# sample list instead.
meta_files <- file.path(dirs, "aux_info", "meta_info.json")
missing_meta <- meta_files[!file.exists(meta_files)]
has_boot <- vapply(dirs, function(d) {
  bd <- file.path(d, "aux_info", "bootstrap")
  (dir.exists(bd) && length(list.files(bd)) > 0) ||
    file.exists(file.path(d, "aux_info", "eq_classes.txt"))
}, logical(1))
if (length(missing_meta) > 0 || any(!has_boot)) {
  bad <- unique(c(dirname(dirname(missing_meta)), dirs[!has_boot]))
  stop("catchSalmon pre-flight failed — ", length(bad), " sample(s) lack ",
       "aux_info/meta_info.json and/or bootstrap resamples (aux_info/ is ",
       "required for RTA overdispersion estimation):\n  ",
       paste(utils::head(bad, 20), collapse = "\n  "),
       if (length(bad) > 20) sprintf("\n  ... and %d more", length(bad) - 20) else "",
       "\nIf aux_info/ was cleaned to save space, the affected samples must ",
       "be re-quantified with Salmon --numBootstraps before QU correction.")
}

cat("================================================================\n")
cat("qu_correct_salmon.R\n")
cat("  samples    :", length(samples), "\n")
cat("  salmon-dir :", salmon_dir, "\n")
cat("  ref-anno   :", ref_anno_file, "\n")
cat("  out-dir    :", out_dir, "\n")
cat("  chunk-size :", chunk_size, "\n")
cat("================================================================\n")

# ---- Build tx2gene from the reference GTF (once, before chunk loop) --------
cat("[1/5] Building tx2gene from GTF (rtracklayer::import)\n")
tx2gene <- .build_tx2gene(ref_anno_file)
cat("       tx2gene rows:", nrow(tx2gene), "\n")

# Carry Length / EffectiveLength from the original quant.sf files (per-sample).
# Build a transcript -> (Length, EffectiveLength) table from the first sample's
# quant.sf (these are annotation-derived and identical across samples for the
# same index). Read once before the chunk loop.
first_qsf <- read.delim(files[1], stringsAsFactors = FALSE)
length_tbl <- first_qsf[, c("Name", "Length", "EffectiveLength")]
rownames(length_tbl) <- length_tbl$Name
rm(first_qsf)

# ---- Split samples into chunks ---------------------------------------------
# catchSalmon + tximport hold all bootstrap replicates in memory; chunking
# bounds peak RSS ~proportionally to chunk-size. The RTA overdispersion is a
# transcript-level property (read-to-transcript ambiguity), so estimating it
# from chunk-size samples x 20 bootstraps is statistically stable.
if (chunk_size > 0 && length(samples) > chunk_size) {
  chunk_idx_list <- split(seq_along(samples), ceiling(seq_along(samples) / chunk_size))
} else {
  chunk_idx_list <- list(seq_along(samples))
}
n_chunks <- length(chunk_idx_list)
cat("Processing ", length(samples), " samples in ", n_chunks, " chunk(s)\n", sep = "")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (ci in seq_along(chunk_idx_list)) {
  idx <- chunk_idx_list[[ci]]
  chunk_samples <- samples[idx]
  chunk_dirs    <- dirs[idx]
  chunk_files   <- files[idx]

  cat(sprintf("[chunk %d/%d] %d samples (%s ... %s)\n",
              ci, n_chunks, length(chunk_samples),
              chunk_samples[1], chunk_samples[length(chunk_samples)]))

  # ---- catchSalmon: estimate RTA overdispersion from bootstraps -----------
  cat("  [2/5] edgeR::catchSalmon\n")
  s <- edgeR::catchSalmon(chunk_dirs)
  overdisp <- s$annotation$Overdispersion
  overdisp <- pmax(overdisp, 1)
  cat("       overdispersion: min=", round(min(overdisp), 3),
      " median=", round(median(overdisp), 3),
      " max=", round(max(overdisp), 3), "\n", sep = "")

  # ---- tximport: import quant.sf ------------------------------------------
  cat("  [3/5] tximport\n")
  txi <- tximport::tximport(chunk_files, type = "salmon", tx2gene = tx2gene,
                            txOut = TRUE, dropInfReps = TRUE)

  # Align row order of overdispersion to txi$counts rows.
  if (!all(rownames(txi$counts) == rownames(s$counts))) {
    overdisp <- overdisp[match(rownames(txi$counts), rownames(s$counts))]
    if (any(is.na(overdisp))) {
      n_na <- sum(is.na(overdisp))
      warning(n_na, " transcripts in tximport output not found in catchSalmon output; ",
              "setting their overdispersion to 1 (no scaling).")
      overdisp[is.na(overdisp)] <- 1
    }
  }

  # ---- QU correction: divide counts by overdispersion ---------------------
  cat("  [4/5] QU correction\n")
  scaled_counts <- txi$counts / overdisp

  # ---- Write adjusted quant.sf per sample ---------------------------------
  cat("  [5/5] Writing adjusted quant.sf\n")
  for (i in seq_along(chunk_samples)) {
    sm <- chunk_samples[i]
    sm_out_dir <- file.path(out_dir, sm)
    dir.create(sm_out_dir, showWarnings = FALSE, recursive = TRUE)
    out_file <- file.path(sm_out_dir, "quant.sf")

    numreads <- scaled_counts[, i]
    tpm <- numreads / sum(numreads) * 1e6

    tx_ids <- rownames(scaled_counts)
    len_df <- length_tbl[tx_ids, c("Length", "EffectiveLength")]

    qsf <- data.frame(
      Name            = tx_ids,
      Length          = len_df$Length,
      EffectiveLength = len_df$EffectiveLength,
      TPM             = tpm,
      NumReads        = numreads,
      stringsAsFactors = FALSE
    )
    qsf$TPM[is.nan(qsf$TPM) | is.infinite(qsf$TPM)] <- 0
    qsf$NumReads[is.nan(qsf$NumReads) | is.infinite(qsf$NumReads)] <- 0

    write.table(qsf, file = out_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)
  }

  cat(sprintf("  chunk %d/%d done (%d samples written)\n", ci, n_chunks, length(chunk_samples)))

  # Free all large objects before the next chunk.
  rm(s, txi, scaled_counts, overdisp)
  gc()
}

cat("       wrote ", length(samples), " adjusted quant.sf files total\n", sep = "")
cat("================================================================\n")
cat("qu_correct_salmon.R DONE\n")
