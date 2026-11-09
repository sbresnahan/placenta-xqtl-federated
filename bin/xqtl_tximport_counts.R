#!/usr/bin/env Rscript
# Optional XQTL_R_LIBS (colon-separated) is prepended to the R library path.
# =============================================================================
# tximport_counts.R
# -----------------------------------------------------------------------------
# Build count-scale expression BEDs for the TMM -> VST normalization chain
# (PsychENCODE/isoTWAS convention; Bhattacharya et al., Nat Genet 2023,
# PMC10703692). These BEDs are the count input that stage 5 normalizes with
# edgeR TMM -> DESeq2 VST (terminal VST; no per-gene INT for the two
# count-based expression modalities).
#
# Outputs (written to <out-dir>, i.e. <cohort_dir>/output/unnorm/):
#   expression_counts.bed
#       Gene-level counts via tximport(countsFromAbundance="lengthScaledTPM")
#       from the ORIGINAL Salmon quants (<salmon-dir>). This is the
#       PsychENCODE repo convention: lengthScaledTPM gene counts are on a
#       count-like scale that is comparable across genes (gene length divided
#       out), suitable for TMM -> VST.
#   isoform_expression_counts.bed
#       Transcript-level counts via tximport(txOut=TRUE) from the
#       QU-CORRECTED Salmon quants (<qu-dir>, i.e. intermediate/expression_qu).
#       With txOut=TRUE and countsFromAbundance="no", txi$counts are the
#       (QU-corrected) NumReads — the count-like values for the isoform
#       abundance modality.
#
# BED schema (matches assemble_bed.py):
#   #chr, start, end, phenotype_id, <one column per sample, samples.txt order>
#   Coordinates are the GENE TSS (0-based half-open): for + strand the TSS is
#   the gene start, for - strand the gene end; start = TSS - 1, end = TSS.
#   Isoform phenotype_id = gene_id__transcript_id (gene TSS coordinates).
#
# ALL features are emitted — no detection filter here. The stage-5
# normalization scripts apply the devBrain detection filter on the TPM
# matrices (TPM > 0.1 in > 25% of stratum samples) and then subset these
# count matrices to the kept features, so the VST feature set is identical
# to the TPM-path feature set.
#
# These BEDs are pooling intermediates only: they are NOT copied to output/,
# NOT bgzipped/tabix-indexed, and are consumed by
# pool_expression_within_ancestry.py / pool_modalities_within_ancestry.py
# (modalities "expression_counts" / "isoform_expression_counts").
#
# Usage
# -----
#   Rscript tximport_counts.R \
#       --samples    <samples.txt> \
#       --salmon-dir <cohort_dir>/intermediate/expression \
#       --qu-dir     <cohort_dir>/intermediate/expression_qu \
#       --ref-anno   <HPLRv2.0.annotated.PANTRY.gtf> \
#       --out-dir    <cohort_dir>/output/unnorm
#
# Requires in the seadragon R library: tximport, rtracklayer (and deps).
# edgeR/DESeq2 are NOT needed at this step (they enter at stage 5).
# =============================================================================

local({ x <- Sys.getenv("XQTL_R_LIBS"); if (nzchar(x)) .libPaths(c(strsplit(x, ":")[[1]], .libPaths())) })

suppressPackageStartupMessages(library(optparse))

# =============================================================================
# Helper: build tx2gene (transcript_id -> gene_id) from a GTF via rtracklayer
# (same implementation as qu_correct_salmon.R — duplicated per repo convention)
# =============================================================================
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

# =============================================================================
# Helper: gene TSS table from a GTF, matching assemble_bed.py::load_tss
#   + strand: TSS = gene start;  - strand: TSS = gene end
#   BED start = TSS - 1, end = TSS (0-based half-open); sorted by chr, start
# =============================================================================
.build_gene_tss <- function(gtf_file) {
  gr <- rtracklayer::import(gtf_file)
  gr <- gr[gr$type == "gene"]
  if (length(gr) == 0) {
    stop("No gene features with type=='gene' found in GTF: ", gtf_file)
  }
  tss <- ifelse(as.character(GenomicRanges::strand(gr)) == "+",
                GenomicRanges::start(gr), GenomicRanges::end(gr))
  df <- data.frame(
    chr     = as.character(GenomicRanges::seqnames(gr)),
    start   = as.integer(tss - 1),
    end     = as.integer(tss),
    gene_id = gr$gene_id,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$gene_id) & nzchar(df$gene_id), ]
  df <- df[order(df$chr, df$start), ]
  rownames(df) <- NULL
  df
}

# =============================================================================
# Helper: write a BED (coords + count matrix), %g-style float formatting to
# match the pandas float_format='%g' used by assemble_bed.py
# =============================================================================
.write_bed <- function(coords_df, counts_mat, out_file) {
  stopifnot(nrow(coords_df) == nrow(counts_mat))
  # MODIFIED: trimws() -- formatC pads every cell to a common width ("     0")
  mat_fmt <- matrix(trimws(formatC(counts_mat, format = "g", digits = 6)),
                    nrow = nrow(counts_mat))
  out <- data.frame(coords_df, mat_fmt, check.names = FALSE,
                    stringsAsFactors = FALSE)
  colnames(out) <- c("#chr", "start", "end", "phenotype_id", colnames(counts_mat))
  write.table(out, file = out_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}

option_list <- list(
  make_option("--samples", type = "character", default = NULL,
              help = "File with one sample ID per line (no header)"),
  make_option("--salmon-dir", type = "character", default = NULL,
              help = "Directory with per-sample ORIGINAL Salmon output subdirs (gene-level counts source)"),
  make_option("--qu-dir", type = "character", default = NULL,
              help = "Directory with per-sample QU-CORRECTED quant.sf subdirs (isoform-level counts source)"),
  make_option("--ref-anno", type = "character", default = NULL,
              help = "Reference GTF with gene + transcript features (tx2gene and gene TSS)"),
  make_option("--out-dir", type = "character", default = NULL,
              help = "Output directory for expression_counts.bed and isoform_expression_counts.bed")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("samples", "salmon-dir", "qu-dir", "ref-anno", "out-dir")
missing <- required[vapply(required, function(k) is.null(opt[[k]]), logical(1))]
if (length(missing) > 0) {
  stop(sprintf("Missing required argument(s): %s. Required: --samples, --salmon-dir, --qu-dir, --ref-anno, --out-dir",
               paste0("--", missing, collapse = ", ")))
}

samples_file  <- opt[["samples"]]
salmon_dir    <- opt[["salmon-dir"]]
qu_dir        <- opt[["qu-dir"]]
ref_anno_file <- opt[["ref-anno"]]
out_dir       <- opt[["out-dir"]]

# ---- Load required packages ------------------------------------------------
for (pkg in c("tximport", "rtracklayer")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required R package '%s' not found in .libPaths(): %s\n",
                 pkg, paste(.libPaths(), collapse = ", ")))
  }
}
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(rtracklayer))

# ---- Read samples ----------------------------------------------------------
samples <- readLines(samples_file)
samples <- trimws(samples)
samples <- samples[nzchar(samples)]
if (length(samples) == 0) stop("No samples found in samples file: ", samples_file)

files_orig <- file.path(salmon_dir, samples, "quant.sf")
files_qu   <- file.path(qu_dir, samples, "quant.sf")
names(files_orig) <- samples
names(files_qu)   <- samples
for (tag in c("original", "QU-corrected")) {
  ff <- if (tag == "original") files_orig else files_qu
  miss <- ff[!file.exists(ff)]
  if (length(miss) > 0) {
    stop("Missing ", tag, " quant.sf files:\n  ",
         paste(utils::head(miss, 20), collapse = "\n  "),
         if (length(miss) > 20) sprintf("\n  ... and %d more", length(miss) - 20) else "")
  }
}

cat("================================================================\n")
cat("tximport_counts.R\n")
cat("  samples    :", length(samples), "\n")
cat("  salmon-dir :", salmon_dir, "\n")
cat("  qu-dir     :", qu_dir, "\n")
cat("  ref-anno   :", ref_anno_file, "\n")
cat("  out-dir    :", out_dir, "\n")
cat("================================================================\n")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Annotation tables (once) ----------------------------------------------
cat("[1/4] Building tx2gene + gene TSS table from GTF\n")
tx2gene <- .build_tx2gene(ref_anno_file)
tss_df  <- .build_gene_tss(ref_anno_file)
cat("       tx2gene rows:", nrow(tx2gene), " gene TSS rows:", nrow(tss_df), "\n")

# ---- Gene-level counts from ORIGINAL quants (lengthScaledTPM) --------------
cat("[2/4] tximport gene-level (countsFromAbundance='lengthScaledTPM') from original quants\n")
txi_gene <- tximport::tximport(files_orig, type = "salmon", tx2gene = tx2gene,
                               txOut = FALSE, dropInfReps = TRUE,
                               countsFromAbundance = "lengthScaledTPM")
counts_gene <- txi_gene$counts
rm(txi_gene)

# Keep genes present in both the count matrix and the TSS annotation;
# order by the (sorted) TSS table.
keep_g <- tss_df$gene_id %in% rownames(counts_gene)
if (sum(keep_g) == 0) stop("No overlap between tximport gene IDs and GTF gene IDs")
tss_g <- tss_df[keep_g, ]
counts_gene <- counts_gene[tss_g$gene_id, , drop = FALSE]
cat("       genes:", nrow(counts_gene), " x samples:", ncol(counts_gene), "\n")

cat("[3/4] Writing expression_counts.bed\n")
.write_bed(data.frame(chr = tss_g$chr, start = tss_g$start,
                      end = tss_g$end, phenotype_id = tss_g$gene_id,
                      stringsAsFactors = FALSE),
           counts_gene, file.path(out_dir, "expression_counts.bed"))
rm(counts_gene, tss_g)
gc()

# ---- Isoform-level counts from QU-corrected quants (txOut) -----------------
cat("[4/4] tximport transcript-level (txOut=TRUE) from QU-corrected quants\n")
txi_tx <- tximport::tximport(files_qu, type = "salmon", txOut = TRUE,
                             dropInfReps = TRUE,
                             countsFromAbundance = "no")
counts_tx <- txi_tx$counts
rm(txi_tx)

# Map transcripts to genes; keep transcripts whose gene has a TSS annotation.
tx2g_f <- tx2gene[tx2gene$transcript_id %in% rownames(counts_tx) &
                  tx2gene$gene_id %in% tss_df$gene_id, ]
if (nrow(tx2g_f) == 0) stop("No overlap between tximport transcript IDs and GTF tx2gene")
# Order by gene TSS order, then transcript_id within gene (deterministic).
gene_rank <- match(tx2g_f$gene_id, tss_df$gene_id)
tx2g_f <- tx2g_f[order(gene_rank, tx2g_f$transcript_id), ]
counts_tx <- counts_tx[tx2g_f$transcript_id, , drop = FALSE]
tss_tx <- tss_df[match(tx2g_f$gene_id, tss_df$gene_id), ]
cat("       transcripts:", nrow(counts_tx), " x samples:", ncol(counts_tx), "\n")

.write_bed(data.frame(chr = tss_tx$chr, start = tss_tx$start,
                      end = tss_tx$end,
                      phenotype_id = paste0(tx2g_f$gene_id, "__", tx2g_f$transcript_id),
                      stringsAsFactors = FALSE),
           counts_tx, file.path(out_dir, "isoform_expression_counts.bed"))

cat("================================================================\n")
cat("tximport_counts.R DONE\n")
cat("  wrote:", file.path(out_dir, "expression_counts.bed"), "\n")
cat("  wrote:", file.path(out_dir, "isoform_expression_counts.bed"), "\n")
