#!/usr/bin/env Rscript
# compare_atlas_phenotypes.R -- concordance of pipeline phenotype BEDs with atlas phenotype BEDs for
# the same samples (development check; the atlas BEDs are private and never committed).
#
#   Rscript tests/scripts/compare_atlas_phenotypes.R PIPELINE_PHENOTYPES_DIR ATLAS_DIR OUT_TSV
#
# PIPELINE_PHENOTYPES_DIR  <outdir>/phenotypes of a phenotyping run (unnormalized BEDs)
# ATLAS_DIR                atlas <modality>.chr22.bed.gz files (tests/fixtures/extract_atlas_refs_chr22.sh)
#
# Per modality: feature counts and ID overlap; median per-feature Spearman correlation across samples
# (valid against atlas values that were rank-normal transformed, since that is monotone within a
# feature); for modalities on a natural scale, median per-sample Spearman across shared features;
# and a sample-identity check (each pipeline sample's best-matching atlas sample).
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 3) stop("usage: compare_atlas_phenotypes.R PIPELINE_PHENOTYPES_DIR ATLAS_DIR OUT_TSV")
ours_dir <- a[1]; atlas_dir <- a[2]; out <- a[3]

# pipeline modality, atlas file stem, natural scale (per-sample correlation meaningful), ID key
mods <- data.table(
  ours  = c("expression", "isoform_expression", "isoform_usage", "alt_tss", "alt_polya",
            "splicing", "intron_retention", "rna_editing", "stability"),
  atlas = c("expression", "isoform_expression", "isoforms", "alt_TSS", "alt_polyA",
            "splicing", "intron_retention", "RNA_editing", "stability"),
  natural = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE))

# Matching keys: intron retention on gene + intron coordinates and strand (atlas IDs are
# IR_<gene>_<chr>_<s>_<e>_<strand>_<n>; see match_ir); isoforms on the transcript ID (some atlas isoform IDs carry
# HPLRG gene IDs where the annotation's gene_id is an Ensembl ID); everything else on the full ID
key <- function(id, mod) {
  if (mod == "intron_retention") return(paste(sub("__.*", "", id), regmatches(id, regexpr("chr[^_]+_[0-9]+_[0-9]+_[+-]", id))))
  if (mod %in% c("isoform_expression", "isoform_usage")) return(sub(".*__", "", id))
  id
}

# The atlas measured intron retention with MAJIQ (intron segments between annotated exons), the
# pipeline with IRFinder-S (whole introns). Each atlas intron is matched to the shortest pipeline
# intron of the same gene and strand that contains it; atlas rows without a match are dropped.
match_ir <- function(o, at) {
  parse <- function(k) {
    p <- regmatches(k, regexec("^(\\S+) (chr[^_]+)_([0-9]+)_([0-9]+)_([+-])$", k))
    data.table(id = k, gene = vapply(p, `[`, "", 2), chr = vapply(p, `[`, "", 3),
               s = as.numeric(vapply(p, `[`, "", 4)), e = as.numeric(vapply(p, `[`, "", 5)),
               strand = vapply(p, `[`, "", 6))
  }
  po <- parse(rownames(o)); pa <- parse(rownames(at))
  hit <- vapply(seq_len(nrow(pa)), function(i) {
    c <- po[gene == pa$gene[i] & chr == pa$chr[i] & strand == pa$strand[i] &
            s <= pa$s[i] & e >= pa$e[i] - 1]
    if (nrow(c) == 0) NA_character_ else c$id[which.min(c$e - c$s)]
  }, "")
  keep <- !is.na(hit) & !duplicated(hit)
  at <- at[keep, , drop = FALSE]
  rownames(at) <- hit[keep]
  at
}

read_bed <- function(f, mod) {
  d <- fread(cmd = paste("gzip -cd", shQuote(f)), na.strings = c("", "NA", "nan"))
  ids <- d[[4]]
  m <- as.matrix(d[, -(1:4)])
  storage.mode(m) <- "double"
  rownames(m) <- vapply(ids, key, "", mod = mod)
  m[!duplicated(rownames(m)), , drop = FALSE]
}

sp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok], method = "spearman")
}

res <- rbindlist(lapply(seq_len(nrow(mods)), function(i) {
  m <- mods[i]
  fo <- file.path(ours_dir, paste0(m$ours, ".bed.gz"))
  fa <- file.path(atlas_dir, paste0(m$atlas, ".chr22.bed.gz"))
  if (!file.exists(fo) || !file.exists(fa)) {
    return(data.table(modality = m$ours, n_pipeline = NA, n_atlas = NA, n_shared = NA))
  }
  o <- read_bed(fo, m$ours); at <- read_bed(fa, m$ours)
  if (m$ours == "intron_retention") at <- match_ir(o, at)
  s <- intersect(colnames(o), colnames(at))
  f <- intersect(rownames(o), rownames(at))
  if (length(f) < 2) {
    return(data.table(modality = m$ours, n_pipeline = nrow(o), n_atlas = nrow(at), n_shared = length(f)))
  }
  o <- o[f, s, drop = FALSE]; at <- at[f, s, drop = FALSE]
  # TPMs: rescale both to each sample's total over the shared features, so the comparison does not
  # depend on how much of the transcriptome each quantification saw (the fixture reads are chr22 only)
  if (m$ours %in% c("expression", "isoform_expression")) {
    o <- sweep(o, 2, colSums(o, na.rm = TRUE), "/"); at <- sweep(at, 2, colSums(at, na.rm = TRUE), "/")
  }
  per_feat <- vapply(f, function(k) sp(o[k, ], at[k, ]), 0)
  per_samp <- if (m$natural) vapply(s, function(j) sp(o[, j], at[, j]), 0) else NA_real_
  # sample identity: best atlas match for each pipeline sample
  # (each feature's values replaced by centred within-feature ranks, so scale does not matter)
  cen <- function(x) t(apply(x, 1, function(r) { q <- rank(r, na.last = "keep"); q - mean(q, na.rm = TRUE) }))
  # on features that vary in both data sets and are not mostly zero (ties at 0 dominate otherwise)
  use <- is.finite(per_feat) & rowMeans(at > 0, na.rm = TRUE) >= 0.5
  cc <- suppressWarnings(cor(cen(o[use, , drop = FALSE]), cen(at[use, , drop = FALSE]),
                             use = "pairwise.complete.obs", method = "spearman"))
  ident <- sum(colnames(at)[apply(cc, 1, which.max)] == rownames(cc))
  data.table(modality = m$ours, n_pipeline = nrow(read_bed(fo, m$ours)), n_atlas = nrow(read_bed(fa, m$ours)) ,
             n_shared = length(f), n_feat_tested = sum(is.finite(per_feat)),
             median_feature_rho = round(median(per_feat, na.rm = TRUE), 3),
             frac_feature_rho_gt_0.8 = round(mean(per_feat > 0.8, na.rm = TRUE), 3),
             median_sample_rho = round(median(per_samp, na.rm = TRUE), 3),
             samples = length(s), samples_self_matched = ident)
}), fill = TRUE)
fwrite(res, out, sep = "\t", na = "NA")
print(res)
