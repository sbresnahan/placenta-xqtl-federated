#!/usr/bin/env Rscript
# Convert central SuShiE prediction weights to FUSION .wgt.RDat files (one per feature and ancestry).
#
# Input (written by `xqtl_central.py twas-inputs`):
#   <twas-dir>/fusion/<ANC>/sushie_weights.txt   CHROM POS A1 A2 TargetID ES   (A1 = effect allele = ALT;
#                                                ES = SuShiE posterior effect, standardized genotype scale)
#   <twas-dir>/otters/<ANC>/exp_anno.txt         CHROM GeneStart GeneEnd TargetID
#   <twas-dir>/otters/<ANC>/eQTLSumStats.txt     ... TargetID N   (N.tot = median N of the feature)
# Output per ancestry:
#   <twas-dir>/fusion/<ANC>/WEIGHTS/<TargetID>.wgt.RDat with wgt.matrix (column "sushie"), snps (6 PLINK .bim
#   columns: CHR, SNP, cM, BP, A1, A2), cv.performance (rsq/pval rows; NA - summary-statistics fine-mapping has
#   no cross-validation), hsq (NA, NA), hsq.pv (NA), N.tot;
#   <twas-dir>/fusion/<ANC>/<ANC>.pos (PANEL WGT ID CHR P0 P1), the FUSION weights list.
# Because cv.performance is NA, run FUSION.assoc_test.R with --force_model sushie.
#
# Usage: xqtl_fusion_wgt.R --twas-dir twas [--panel-prefix placenta_xqtl]

local({ x <- Sys.getenv("XQTL_R_LIBS"); if (nzchar(x)) .libPaths(c(strsplit(x, ":")[[1]], .libPaths())) })

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--twas-dir", type = "character", default = "twas"),
  make_option("--panel-prefix", type = "character", default = "placenta_xqtl")
)))
td <- opt[["twas-dir"]]

safe <- function(x) gsub("[^A-Za-z0-9._+-]", "_", x)

for (wf in Sys.glob(file.path(td, "fusion", "*", "sushie_weights.txt"))) {
  d <- dirname(wf)
  anc <- basename(d)
  w <- fread(wf, sep = "\t", colClasses = list(character = c("CHROM", "A1", "A2", "TargetID")))
  if (nrow(w) == 0) next
  anno_f <- file.path(td, "otters", anc, "exp_anno.txt")
  anno <- if (file.exists(anno_f)) fread(anno_f, sep = "\t", colClasses = list(character = c("CHROM", "TargetID"))) else
    data.table(CHROM = character(), GeneStart = integer(), GeneEnd = integer(), TargetID = character())
  sst_f <- file.path(td, "otters", anc, "eQTLSumStats.txt")
  ntot <- if (file.exists(sst_f)) fread(sst_f, sep = "\t", select = c("TargetID", "N"))[, .(N = median(N)), by = TargetID] else
    data.table(TargetID = character(), N = numeric())
  dir.create(file.path(d, "WEIGHTS"), showWarnings = FALSE)
  pos <- list()
  for (tid in unique(w$TargetID)) {
    x <- w[TargetID == tid][order(POS)]
    snp_id <- paste0("chr", x$CHROM, ":", x$POS, ":", x$A2, ":", x$A1)
    wgt.matrix <- matrix(x$ES, ncol = 1, dimnames = list(snp_id, "sushie"))
    snps <- data.frame(V1 = x$CHROM, V2 = snp_id, V3 = 0, V4 = x$POS, V5 = x$A1, V6 = x$A2, stringsAsFactors = FALSE)
    cv.performance <- matrix(NA_real_, 2, 1, dimnames = list(c("rsq", "pval"), "sushie"))
    hsq <- c(NA_real_, NA_real_)
    hsq.pv <- NA_real_
    N.tot <- if (tid %in% ntot$TargetID) ntot[TargetID == tid, N] else NA_real_
    fn <- paste0(safe(tid), ".wgt.RDat")
    save(wgt.matrix, snps, cv.performance, hsq, hsq.pv, N.tot, file = file.path(d, "WEIGHTS", fn))
    a <- anno[TargetID == tid]
    pos[[tid]] <- data.table(PANEL = paste0(opt[["panel-prefix"]], "_", anc), WGT = file.path("WEIGHTS", fn),
                             ID = tid, CHR = x$CHROM[1],
                             P0 = if (nrow(a)) a$GeneStart[1] else min(x$POS),
                             P1 = if (nrow(a)) a$GeneEnd[1] else max(x$POS))
  }
  fwrite(rbindlist(pos), file.path(d, paste0(anc, ".pos")), sep = "\t")
  message(sprintf("[fusion_wgt] %s: %d weight files", anc, length(pos)))
}
