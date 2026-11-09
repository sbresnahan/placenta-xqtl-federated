#!/usr/bin/env Rscript
# check_central_ivw.R -- independent check of the central nominal meta-analysis and the FUSION weight files.
#
#   check_central_ivw.R EXTRACT_DIR CENTRAL_DIR
#
# EXTRACT_DIR is written by `check_central.py extract`: ivw_inputs.tsv holds the raw file A rows of every
# contribution (read straight from the two study exports) for a set of meta rows, ivw_meta.tsv the pipeline's
# meta rows. Fixed-effect and DerSimonian-Laird random-effects IVW, Cochran's Q, I^2 and tau^2 are recomputed
# here per scope (ALL = every contribution; <ANC> = that ancestry's contributions) and compared with the
# pipeline (relative tolerance 1e-6). The FUSION .wgt.RDat files under CENTRAL_DIR/twas/fusion are loaded
# and checked for the objects and shapes FUSION.assoc_test.R reads. Exit status 1 on any failure.
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: check_central_ivw.R EXTRACT_DIR CENTRAL_DIR")
fails <- 0L
ok <- function(cond, msg) {
  cat(if (isTRUE(cond)) "PASS  " else "FAIL  ", msg, "\n", sep = "")
  if (!isTRUE(cond)) fails <<- fails + 1L
}
close_enough <- function(a, b, rel = 1e-6, abs_tol = 1e-12) {
  both_na <- is.na(a) & is.na(b)
  both_na | (!is.na(a) & !is.na(b) & abs(a - b) <= abs_tol + rel * pmax(abs(a), abs(b)))
}

ivw <- function(b, s, contrib) {
  w <- 1 / s^2
  k <- length(b)
  bf <- sum(w * b) / sum(w)
  sf <- sqrt(1 / sum(w))
  q <- sum(w * (b - bf)^2)
  cc <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (k >= 2 && cc > 0) max(0, (q - (k - 1)) / cc) else 0
  wr <- 1 / (s^2 + tau2)
  br <- sum(wr * b) / sum(wr)
  sr <- sqrt(1 / sum(wr))
  list(k = k, beta_fe = bf, se_fe = sf, p_fe = 2 * pnorm(-abs(bf / sf)),
       beta_re = br, se_re = sr, p_re = 2 * pnorm(-abs(br / sr)), tau2 = tau2,
       q_het = if (k >= 2) q else NA_real_,
       i2 = if (k < 2) NA_real_ else if (q > 0) max(0, (q - (k - 1)) / q) else 0,
       contrib_mask = sum(2^contrib))
}

# ---- nominal IVW ------------------------------------------------------------------------------------
inp <- fread(file.path(args[1], "ivw_inputs.tsv"))
met <- fread(file.path(args[1], "ivw_meta.tsv"))
stats <- c("beta_fe", "se_fe", "p_fe", "beta_re", "se_re", "p_re", "tau2", "q_het", "i2")
setkey(inp, unit, phenotype_id, variant_id)
bad <- list()
for (i in seq_len(nrow(met))) {
  m <- met[i]
  x <- inp[.(m$unit, m$phenotype_id, m$variant_id), nomatch = NULL]
  if (m$scope != "ALL") x <- x[ancestry == m$scope]
  if (nrow(x) == 0) {
    bad[[length(bad) + 1]] <- data.table(row = i, stat = "inputs", pipeline = m$k, recomputed = 0)
    next
  }
  r <- ivw(x$beta, x$se, x$contrib)
  if (r$k != m$k || r$contrib_mask != m$contrib_mask) {
    bad[[length(bad) + 1]] <- data.table(row = i, stat = "k/contrib_mask", pipeline = m$k, recomputed = r$k)
  }
  for (s in stats) {
    if (!close_enough(m[[s]], r[[s]])) {
      bad[[length(bad) + 1]] <- data.table(row = i, stat = s, pipeline = m[[s]], recomputed = r[[s]])
    }
  }
}
bad <- rbindlist(bad)
het <- met[!is.na(tau2) & tau2 > 0, .N]
ok(nrow(bad) == 0,
   sprintf("IVW recomputed for %d meta rows (%s; %d with tau^2 > 0): FE/RE beta, se, p, Q, I^2, tau^2 agree",
           nrow(met), paste(met[, .N, by = scope][, sprintf("%s %d", scope, N)], collapse = ", "), het))
if (nrow(bad)) print(head(merge(bad, met[, .(row = .I, unit, scope, phenotype_id, variant_id)], by = "row"), 20))
ok(het > 0, sprintf("the recomputed rows include between-study heterogeneity (%d rows with tau^2 > 0)", het))
ok(any(met[, k] == 2) && any(met[scope == "ALL", k] == 4),
   sprintf("rows span k = %d..%d contributions", min(met$k), max(met$k)))

# ---- FUSION weight files ----------------------------------------------------------------------------
wfiles <- Sys.glob(file.path(args[2], "twas", "fusion", "*", "WEIGHTS", "*.wgt.RDat"))
problems <- character(0)
for (f in wfiles) {
  e <- new.env()
  load(f, envir = e)
  why <- character(0)
  if (!all(c("wgt.matrix", "snps") %in% ls(e))) {
    why <- "missing wgt.matrix or snps"
  } else {
    w <- e$wgt.matrix
    if (!is.matrix(w) || nrow(w) != nrow(e$snps)) why <- c(why, "wgt.matrix rows != snps rows")
    if (!identical(colnames(w), "sushie")) why <- c(why, "column is not 'sushie'")
    if (anyNA(w)) why <- c(why, "NA weights")
    if (ncol(e$snps) != 6) why <- c(why, "snps does not have 6 columns")
    if (!identical(rownames(w), as.character(e$snps[[2]]))) why <- c(why, "wgt.matrix rownames != snps IDs")
  }
  if (length(why)) problems <- c(problems, sprintf("%s: %s", basename(f), paste(why, collapse = "; ")))
}
ok(length(wfiles) > 0 && length(problems) == 0,
   sprintf("%d FUSION .wgt.RDat files load with wgt.matrix ('sushie') aligned to a 6-column snps table",
           length(wfiles)))
if (length(problems)) writeLines(head(problems, 20))

cat(sprintf("check_central_ivw.R: %d failure(s)\n", fails))
quit(status = if (fails) 1L else 0L)
