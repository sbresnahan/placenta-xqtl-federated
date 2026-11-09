#!/usr/bin/env Rscript
# Optional XQTL_R_LIBS (colon-separated) is prepended to the R library path.
# Script 3: Run MuSiC deconvolution on a cohort's bulk TPM matrix.
#
# Input:
#   - placenta_music_reference.rds  (SingleCellExperiment, from Script 1)
#   - <cohort>_bulk_tpm_symbols.tsv (genes x samples TPM, from Script 2)
#
# Output:
#   - <cohort>_cell_proportions_full.tsv       (samples x 27 cell types)
#   - <cohort>_cell_proportions_collapsed.tsv  (samples x 8 collapsed types)
#   - <cohort>_maternal_flag.tsv               (samples with >10% maternal fraction)
#
# Usage:
#   Rscript run_music_deconvolution.R \
#     --reference  placenta_music_reference.rds \
#     --bulk       cohort1_bulk_tpm_symbols.tsv \
#     --cohort     cohort1 \
#     --output-dir /mnt/results/deconvolution \
#     --maternal-threshold 0.10

local({ x <- Sys.getenv("XQTL_R_LIBS"); if (nzchar(x)) .libPaths(c(strsplit(x, ":")[[1]], .libPaths())) })

suppressPackageStartupMessages({
  library(optparse)
  library(MuSiC)
  library(SingleCellExperiment)
  library(Matrix)
  library(data.table)
})

option_list <- list(
  make_option("--reference", type="character", default=NULL,
              help="Path to placenta_music_reference.rds [required]"),
  make_option("--bulk", type="character", default=NULL,
              help="Path to <cohort>_bulk_tpm_symbols.tsv [required]"),
  make_option("--cohort", type="character", default="cohort",
              help="Cohort name (used in output filenames)"),
  make_option("--output-dir", type="character", default=".",
              help="Directory for output TSVs"),
  make_option("--maternal-threshold", type="double", default=0.10,
              help="Maternal fraction threshold for flagging (default 0.10)"),
  make_option("--normalize", action="store_true", type="logical", default=FALSE,
              help="Enable MuSiC normalize=TRUE (cross-platform normalization)")
)
opt <- parse_args(OptionParser(option_list=option_list))

# optparse stores hyphenated option names with hyphens preserved;
# access via [[ ]] bracket notation.
ref_path       <- opt[["reference"]]
bulk_path      <- opt[["bulk"]]
cohort         <- opt[["cohort"]]
outdir         <- opt[["output-dir"]]
mat_threshold  <- opt[["maternal-threshold"]]
do_normalize   <- opt[["normalize"]]

if (is.null(ref_path) || is.null(bulk_path)) {
  stop("--reference and --bulk are required")
}

dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

cat("=== MuSiC deconvolution:", cohort, "===\n")
cat("Reference:", ref_path, "\n")
cat("Bulk:     ", bulk_path, "\n")
cat("Output dir:", outdir, "\n\n")

# --- Load reference SCE ---
cat("Loading reference SCE...\n")
sce <- readRDS(ref_path)
cat("  Reference:", paste(dim(sce), collapse=" x "),
    "| types:", nlevels(sce$cell.type),
    "| subjects:", nlevels(sce$subject), "\n\n")

# --- Load bulk TPM matrix ---
cat("Loading bulk TPM matrix...\n")
bulk <- fread(bulk_path, header=TRUE, sep="\t", data.table=FALSE)
genes <- bulk[[1]]
bulk <- as.matrix(bulk[, -1])
rownames(bulk) <- genes
cat("  Bulk:", paste(dim(bulk), collapse=" x "), "(genes x samples)\n")

# --- Intersect genes with reference ---
ref_genes <- rownames(sce)
shared <- intersect(rownames(bulk), ref_genes)
cat("  Gene overlap:", length(shared), "/", length(ref_genes),
    "(", round(100 * length(shared) / length(ref_genes), 1), "% of reference )\n\n")

if (length(shared) < 100) {
  stop("Insufficient gene overlap (<100 genes). Check gene symbol mapping.")
}

# Subset to shared genes
bulk_sub <- bulk[shared, , drop=FALSE]
sce_sub <- sce[shared, ]

# --- Run MuSiC ---
cat("Running music_prop...\n")
est <- music_prop(
  bulk.mtx  = bulk_sub,
  sc.sce    = sce_sub,
  clusters  = "cell.type",
  samples   = "subject",
  select.ct = levels(sce_sub$cell.type),
  verbose   = TRUE,
  normalize = do_normalize
)
cat("  Done. Est prop dim:", paste(dim(est$Est.prop.weighted), collapse=" x "), "\n\n")

# --- Extract proportions (samples x cell types) ---
prop_full <- est$Est.prop.weighted
# Verify proportions sum to ~1
rowsums <- rowSums(prop_full)
cat("Proportion row sums: min=", round(min(rowsums), 4),
    " max=", round(max(rowsums), 4),
    " mean=", round(mean(rowsums), 4), "\n\n")

# --- Write full 27-type proportions ---
full_path <- file.path(outdir, paste0(cohort, "_cell_proportions_full.tsv"))
prop_full_df <- data.frame(sample = rownames(prop_full), prop_full, check.names=FALSE)
write.table(prop_full_df, full_path, sep="\t", row.names=FALSE, quote=FALSE)
cat("Wrote", full_path, "\n")

# --- Collapse to 8 cell-type groups ---
collapse_map <- list(
  Syncytiotrophoblast = c("Fetal Syncytiotrophoblast"),
  Cytotrophoblast = c("Fetal Cytotrophoblasts", "Fetal Proliferative Cytotrophoblasts"),
  EVT = c("Fetal Extravillous Trophoblasts"),
  Endothelial = c("Fetal Endothelial Cells"),
  Fibroblast_Stromal = c("Fetal Fibroblasts", "Fetal Mesenchymal Stem Cells"),
  Hofbauer = c("Fetal Hofbauer Cells"),
  Fetal_Immune = c(
    "Fetal CD14+ Monocytes", "Fetal B Cells", "Fetal CD8+ Activated T Cells",
    "Fetal Naive CD4+ T Cells", "Fetal Naive CD8+ T Cells",
    "Fetal Memory CD4+ T Cells", "Fetal Natural Killer T Cells",
    "Fetal GZMB+ Natural Killer", "Fetal GZMK+ Natural Killer",
    "Fetal Plasmacytoid Dendritic Cells", "Fetal Nucleated Red Blood Cells"
  ),
  Maternal = c(
    "Maternal B Cells", "Maternal CD14+ Monocytes",
    "Maternal CD8+ Activated T Cells", "Maternal FCGR3A+ Monocytes",
    "Maternal Naive CD4+ T Cells", "Maternal Naive CD8+ T Cells",
    "Maternal Natural Killer Cells", "Maternal Plasma Cells"
  )
)

prop_collapsed <- matrix(0, nrow=nrow(prop_full), ncol=length(collapse_map),
                         dimnames=list(rownames(prop_full), names(collapse_map)))
for (grp in names(collapse_map)) {
  cols <- collapse_map[[grp]]
  cols_present <- intersect(cols, colnames(prop_full))
  if (length(cols_present) > 0) {
    prop_collapsed[, grp] <- rowSums(prop_full[, cols_present, drop=FALSE])
  }
}

# --- Compute maternal fraction and flag ---
maternal_types <- collapse_map[["Maternal"]]
maternal_types_present <- intersect(maternal_types, colnames(prop_full))
maternal_frac <- rowSums(prop_full[, maternal_types_present, drop=FALSE])
flagged <- names(maternal_frac)[maternal_frac > mat_threshold]

cat("\nMaternal fraction: min=", round(min(maternal_frac), 4),
    " median=", round(median(maternal_frac), 4),
    " max=", round(max(maternal_frac), 4), "\n")
cat("Samples flagged (>", mat_threshold*100, "% maternal):",
    length(flagged), "/", length(maternal_frac), "\n\n")

# --- Write collapsed proportions ---
collapsed_path <- file.path(outdir, paste0(cohort, "_cell_proportions_collapsed.tsv"))
prop_collapsed_df <- data.frame(sample = rownames(prop_collapsed), prop_collapsed, check.names=FALSE)
write.table(prop_collapsed_df, collapsed_path, sep="\t", row.names=FALSE, quote=FALSE)
cat("Wrote", collapsed_path, "\n")

# --- Write maternal flag ---
flag_path <- file.path(outdir, paste0(cohort, "_maternal_flag.tsv"))
flag_df <- data.frame(
  sample = names(maternal_frac),
  maternal_fraction = round(maternal_frac, 4),
  flagged = maternal_frac > mat_threshold
)
write.table(flag_df, flag_path, sep="\t", row.names=FALSE, quote=FALSE)
cat("Wrote", flag_path, "\n")

# --- Summary ---
cat("\n=== Summary:", cohort, "===\n")
cat("Samples:", nrow(prop_full), "\n")
cat("Full cell types:", ncol(prop_full), "\n")
cat("Collapsed cell types:", ncol(prop_collapsed), "\n")
cat("Gene overlap:", length(shared), "/", length(ref_genes), "\n")
cat("Flagged samples:", length(flagged), "\n")
cat("\nCollapsed proportion medians:\n")
print(round(apply(prop_collapsed, 2, median), 4))
cat("\n=== Deconvolution complete:", cohort, "===\n")
