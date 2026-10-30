#!/usr/bin/env Rscript
# xqtl_hcp.R -- hidden covariates with prior (HCP) for one normalized phenotype BED.
#
#   xqtl_hcp.R --bed EUR_expression.bed.gz --qc-metrics qc_metrics.tsv --samples EUR.samples.tsv \
#              --k 10 --out EUR_expression.hcp_k10.tsv [--max-phenotypes 20000 --seed 20261029]
#
# Atlas hcp_from_matrix.R: QC metrics (keyed by rnaseq_id, column 'sample') are averaged over
# technical replicates per genotype_id, median-imputed, zero-variance metrics dropped and
# metrics with |r| > 0.9 pruned iteratively (most partners first); phenotypes with variance
# < 1e-10 are dropped and at most --max-phenotypes are sampled (seeded). Both matrices are
# centred and scaled to unit norm; Rhcpp::hcp(k, lambda1 = 0.5, lambda2 = 1, lambda3 = 1,
# iter = 100). Output: rows HCP_1..HCP_k, columns = genotype IDs (k = 0 writes the header only).
suppressPackageStartupMessages({ library(optparse); library(data.table) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--bed"), make_option("--qc-metrics", dest = "qc"), make_option("--samples"),
  make_option("--k", type = "integer"), make_option("--out"),
  make_option("--max-phenotypes", dest = "max_ph", type = "integer", default = 20000L),
  make_option("--seed", type = "integer", default = 20261029L),
  make_option("--lambda1", type = "double", default = 0.5),
  make_option("--lambda2", type = "double", default = 1),
  make_option("--lambda3", type = "double", default = 1),
  make_option("--qc-cor-threshold", dest = "qc_cor", type = "double", default = 0.9)
)))
for (k in c("bed", "qc", "samples", "k", "out")) if (is.null(opt[[k]])) stop("xqtl_hcp.R: --", k, " is required")

rd <- function(p) if (grepl("\\.gz$", p)) fread(cmd = paste("gzip -dc", shQuote(p)), check.names = FALSE) else fread(p, check.names = FALSE)
bed <- rd(opt$bed)
ids <- names(bed)[-(1:4)]
if (opt$k == 0) {
  writeLines(paste(c("covariate", ids), collapse = "\t"), opt$out)
  quit(save = "no", status = 0)
}
if (opt$max_ph > 0 && nrow(bed) > opt$max_ph) {
  set.seed(opt$seed); bed <- bed[sort(sample(nrow(bed), opt$max_ph))]
}
Y <- t(as.matrix(bed[, ..ids])); rownames(Y) <- ids
v <- apply(Y, 2, var, na.rm = TRUE)
Y <- Y[, is.finite(v) & v >= 1e-10, drop = FALSE]
if (opt$k >= min(dim(Y))) stop(sprintf("xqtl_hcp.R: k = %d is not smaller than min(samples, phenotypes) = %d", opt$k, min(dim(Y))))

qc <- fread(opt$qc, check.names = FALSE)
scol <- if ("sample" %in% names(qc)) "sample" else names(qc)[1]
smp <- fread(opt$samples, colClasses = "character")
qc[[scol]] <- as.character(qc[[scol]])
qc <- merge(qc, smp[, .(rnaseq_id, genotype_id)], by.x = scol, by.y = "rnaseq_id")
metrics <- setdiff(names(qc), c(scol, "genotype_id"))
metrics <- metrics[vapply(qc[, ..metrics], is.numeric, logical(1))]
Z <- qc[, lapply(.SD, mean, na.rm = TRUE), by = genotype_id, .SDcols = metrics]
Z <- as.data.frame(Z[match(ids, genotype_id), ..metrics]); rownames(Z) <- ids
for (m in metrics) { x <- Z[[m]]; x[!is.finite(x)] <- NA; med <- median(x, na.rm = TRUE); x[is.na(x)] <- if (is.na(med)) 0 else med; Z[[m]] <- x }
Z <- Z[, vapply(Z, function(x) sd(x) > 0, logical(1)), drop = FALSE]
if (ncol(Z) > 1) {
  cm <- cor(as.matrix(Z)); diag(cm) <- 0
  while (any(abs(cm) > opt$qc_cor) && ncol(cm) > 1) {
    worst <- names(which.max(colSums(abs(cm) > opt$qc_cor)))
    cm <- cm[rownames(cm) != worst, colnames(cm) != worst, drop = FALSE]
  }
  Z <- Z[, colnames(cm), drop = FALSE]
}
if (ncol(Z) == 0) stop("xqtl_hcp.R: no usable QC metrics")
stdz <- function(m) { m <- scale(as.matrix(m), scale = FALSE); s <- sqrt(colSums(m^2)); s[s == 0] <- 1; sweep(m, 2, s, "/") }
suppressPackageStartupMessages(library(Rhcpp))
res <- hcp(Z = stdz(Z), Y = stdz(Y), k = opt$k, lambda1 = opt$lambda1, lambda2 = opt$lambda2,
           lambda3 = opt$lambda3, iter = 100, stand = FALSE, log = FALSE, fast = FALSE, verbose = FALSE)
W <- res$W
colnames(W) <- sprintf("HCP_%d", seq_len(ncol(W)))
out <- data.table(covariate = colnames(W), t(W))
setnames(out, c("covariate", ids))
fwrite(out, opt$out, sep = "\t", quote = FALSE)
message(sprintf("%s: k = %d, %d phenotypes x %d samples, %d QC metrics (%s)", opt$out, opt$k, ncol(Y), nrow(Y),
                ncol(Z), paste(colnames(Z), collapse = ", ")))
