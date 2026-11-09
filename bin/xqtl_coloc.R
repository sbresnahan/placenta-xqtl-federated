#!/usr/bin/env Rscript
# Colocalization of central xQTL loci with GWAS summary statistics (central mode).
#
# Each task (one xQTL locus x one ancestry-matched GWAS trait) is prepared by
# `xqtl_central.py coloc-prep`: <task>.tsv (variant_id, beta_x, se_x, eaf_x, beta_g, se_g, n_g, pos; GWAS
# alleles already aligned to the xQTL effect allele) and <task>.ld.bin (float64 little-endian, nsnps x nsnps,
# the ancestry-matched N-weighted xQTL LD from file C, also used for the GWAS side).
#
# Methods
#   SuSiE-coloc : coloc::runsusie on both traits with the shared LD, then coloc::coloc.susie; when either trait
#                 has no credible set the task falls back to coloc::coloc.abf (single causal variant), and the
#                 `method` column says which was used. A pair colocalizes when PP.H4 >= --pp-h4.
#   ColocBoost  : optional (--colocboost), colocboost(sumstat = list(xQTL, GWAS), LD = R); colocalization
#                 confidence sets from get_colocboost_summary()$cos_summary.
#
# Usage:
#   xqtl_coloc.R --tasks coloc_tasks.tsv --task-dir coloc_in --pp-h4 0.7 --L 10 \
#                [--colocboost --M 500] [--task-ids a,b,c] --out-prefix shard1
# Writes <prefix>.coloc.tsv (one row per task, best credible-set pair), <prefix>.coloc_pairs.tsv (all pairs)
# and, with --colocboost, <prefix>.colocboost.tsv.

local({ x <- Sys.getenv("XQTL_R_LIBS"); if (nzchar(x)) .libPaths(c(strsplit(x, ":")[[1]], .libPaths())) })

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(coloc)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--tasks", type = "character"),
  make_option("--task-dir", type = "character", default = "."),
  make_option("--task-ids", type = "character", default = "", help = "comma-separated subset (default: all prepared)"),
  make_option("--pp-h4", type = "double", default = 0.7),
  make_option("--L", type = "integer", default = 10),
  make_option("--colocboost", action = "store_true", default = FALSE),
  make_option("--M", type = "integer", default = 500),
  make_option("--out-prefix", type = "character")
)))
stopifnot(!is.null(opt$tasks), !is.null(opt[["out-prefix"]]))

tasks <- fread(opt$tasks, sep = "\t", colClasses = list(character = c("prop_cases")))
tasks <- tasks[status == "prepared"]
if (nzchar(opt[["task-ids"]])) tasks <- tasks[task_id %in% strsplit(opt[["task-ids"]], ",")[[1]]]

pp_cols <- c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")

n_cs <- function(s) if (is.null(s$sets$cs)) 0L else length(s$sets$cs)

one_task <- function(t) {
  x <- fread(file.path(opt[["task-dir"]], paste0(t$task_id, ".tsv")), sep = "\t")
  p <- nrow(x)
  R <- matrix(readBin(file.path(opt[["task-dir"]], paste0(t$task_id, ".ld.bin")), "double", n = p * p,
                      size = 8, endian = "little"), p, p)
  R <- (R + t(R)) / 2
  dimnames(R) <- list(x$variant_id, x$variant_id)
  maf <- pmin(x$eaf_x, 1 - x$eaf_x)
  maf[!is.finite(maf) | maf <= 0] <- 1e-4
  # xQTL phenotypes are rank-normalized within site, so sdY = 1 on the meta scale.
  dx <- list(beta = setNames(x$beta_x, x$variant_id), varbeta = setNames(x$se_x^2, x$variant_id),
             N = round(t$n_xqtl), sdY = 1, MAF = setNames(maf, x$variant_id), snp = x$variant_id,
             position = x$pos, LD = R, type = "quant")
  gtype <- if (identical(t$trait_type, "cc") || identical(t$trait_type, "binary")) "cc" else "quant"
  ng <- round(stats::median(x$n_g, na.rm = TRUE))
  dg <- list(beta = setNames(x$beta_g, x$variant_id), varbeta = setNames(x$se_g^2, x$variant_id),
             N = ng, MAF = setNames(maf, x$variant_id), snp = x$variant_id, position = x$pos, LD = R,
             type = gtype)
  if (gtype == "cc") dg$s <- as.numeric(t$prop_cases)
  base <- data.table(task_id = t$task_id, locus_id = t$locus_id, modality = t$modality,
                     phenotype_id = t$phenotype_id, ancestry = t$ancestry, trait_id = t$trait_id, nsnps = p)
  res <- tryCatch({
    sx <- suppressWarnings(runsusie(dx, L = min(opt$L, p), repeat_until_convergence = FALSE, maxit = 1000))
    sg <- suppressWarnings(runsusie(dg, L = min(opt$L, p), repeat_until_convergence = FALSE, maxit = 1000))
    kx <- n_cs(sx); kg <- n_cs(sg)
    if (kx > 0 && kg > 0) {
      cs <- coloc.susie(sx, sg)
      pr <- as.data.table(cs$summary)
      setnames(pr, c("hit1", "hit2", paste0(pp_cols, ".abf")), c("hit_xqtl", "hit_gwas", pp_cols))
      pr[, method := "susie"]
    } else {
      ab <- suppressWarnings(coloc.abf(dx[setdiff(names(dx), "LD")], dg[setdiff(names(dg), "LD")]))
      pr <- as.data.table(as.list(setNames(ab$summary[paste0(pp_cols, ".abf")], pp_cols)))
      pr[, `:=`(hit_xqtl = names(which.max(abs(dx$beta / sqrt(dx$varbeta)))),
                hit_gwas = names(which.max(abs(dg$beta / sqrt(dg$varbeta)))), method = "abf")]
    }
    pr[, `:=`(n_cs_xqtl = kx, n_cs_gwas = kg, status = "ok")]
    pr
  }, error = function(e) data.table(method = NA_character_, status = paste0("error: ", conditionMessage(e))))
  pairs <- cbind(base, res)
  cb <- NULL
  if (opt$colocboost) {
    cb <- tryCatch({
      suppressPackageStartupMessages(library(colocboost))
      ss <- list(data.frame(z = x$beta_x / x$se_x, n = round(t$n_xqtl), variant = x$variant_id),
                 data.frame(z = x$beta_g / x$se_g, n = ng, variant = x$variant_id))
      fit <- suppressMessages(suppressWarnings(colocboost(sumstat = ss, LD = R, M = opt$M,
                                                         outcome_names = c("xqtl", "gwas"))))
      s <- suppressMessages(get_colocboost_summary(fit)$cos_summary)
      if (is.null(s) || nrow(s) == 0) data.table(cos_id = NA_character_, status = "no_cos")
      else cbind(as.data.table(lapply(as.data.frame(s), function(v) if (is.list(v)) vapply(v, paste, "", collapse = ";") else v)),
                 status = "ok")
    }, error = function(e) data.table(cos_id = NA_character_, status = paste0("error: ", conditionMessage(e))))
    cb <- cbind(base, cb)
  }
  list(pairs = pairs, cb = cb)
}

out <- lapply(seq_len(nrow(tasks)), function(i) {
  t0 <- Sys.time()
  r <- one_task(as.list(tasks[i]))
  message(sprintf("[coloc] %s: %s (%.1fs)", tasks$task_id[i], r$pairs$status[1],
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  r
})

pairs <- rbindlist(lapply(out, `[[`, "pairs"), fill = TRUE)
pair_cols <- c("task_id", "locus_id", "modality", "phenotype_id", "ancestry", "trait_id", "nsnps", "method",
               "n_cs_xqtl", "n_cs_gwas", "hit_xqtl", "hit_gwas", pp_cols, "idx1", "idx2", "status")
if (nrow(pairs) == 0) pairs <- as.data.table(setNames(replicate(length(unique(pair_cols)), character(), simplify = FALSE),
                                                    unique(pair_cols)))
for (c in setdiff(unique(pair_cols), names(pairs))) pairs[, (c) := NA]
pairs <- pairs[, unique(pair_cols), with = FALSE]
if (nrow(pairs)) pairs[, coloc := !is.na(PP.H4) & PP.H4 >= opt[["pp-h4"]]]
fwrite(pairs, paste0(opt[["out-prefix"]], ".coloc_pairs.tsv"), sep = "\t", na = "NA")
best <- if (nrow(pairs)) pairs[order(task_id, -PP.H4)][, .SD[1], by = task_id] else pairs
fwrite(best, paste0(opt[["out-prefix"]], ".coloc.tsv"), sep = "\t", na = "NA")
if (opt$colocboost) {
  cb <- rbindlist(lapply(out, `[[`, "cb"), fill = TRUE)
  fwrite(cb, paste0(opt[["out-prefix"]], ".colocboost.tsv"), sep = "\t", na = "NA")
}
