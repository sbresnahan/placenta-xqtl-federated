#!/usr/bin/env Rscript
# xqtl_ancestry.R -- PC-AiR on the merged cohort + 1000 Genomes panel and ancestry assignment.
# Called by xqtl_ancestry.sh; see that script for the overall logic.
#
# Assignment: nearest superpopulation centroid (panel samples) on PCs 1-10 (atlas PC_AiR.R).
# Confidence: posterior probability of the assigned label from a QDA classifier fitted to the
# panel PCs 1-10 (equal priors); samples below --prob-min are "unassigned". Ancestry groups with >= --min-n unique genotyped
# individuals become mapping strata.
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(SNPRelate); library(GENESIS)
  library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--bfile"), make_option("--superpop"), make_option("--metadata"),
  make_option("--cohort-psam", dest = "cohort_psam"),
  make_option("--min-n", dest = "min_n", type = "integer", default = 50L),
  make_option("--prob-min", dest = "prob_min", type = "double", default = 0.8),
  make_option("--groups", default = "AFR,AMR,EAS,EUR,SAS"),
  make_option("--threads", type = "integer", default = 1L),
  make_option("--npcs", type = "integer", default = 10L),
  make_option("--out")
)))
for (k in c("bfile", "superpop", "metadata", "cohort_psam", "out")) {
  if (is.null(opt[[k]])) stop("xqtl_ancestry.R: --", gsub("_", "-", k), " is required")
}
set.seed(20261029)
groups <- strsplit(opt$groups, ",", fixed = TRUE)[[1]]

sp <- fread(opt$superpop, colClasses = "character")
stopifnot(all(c("sample", "superpop") %in% names(sp)))
meta <- fread(opt$metadata, colClasses = list(character = c("rnaseq_id", "genotype_id")))
psam <- fread(opt$cohort_psam, colClasses = "character")
setnames(psam, 1, "IID")
cohort_ids <- psam$IID

gds_file <- "anc_pruned.gds"
snpgdsBED2GDS(paste0(opt$bfile, ".bed"), paste0(opt$bfile, ".fam"), paste0(opt$bfile, ".bim"),
              gds_file, family = FALSE, verbose = FALSE)
gds <- snpgdsOpen(gds_file)
ids <- read.gdsn(index.gdsn(gds, "sample.id"))
panel_ids <- intersect(ids, sp$sample)
if (length(panel_ids) < 50) stop("xqtl_ancestry.R: fewer than 50 panel samples in the merged data")

king <- snpgdsIBDKING(gds, type = "KING-robust", num.thread = opt$threads, verbose = FALSE)
kin <- king$kinship
dimnames(kin) <- list(king$sample.id, king$sample.id)
pc <- pcair(gds, kinobj = kin, divobj = kin, unrel.set = panel_ids,
            eigen.cnt = max(20L, opt$npcs), num.cores = opt$threads, verbose = FALSE)
snpgdsClose(gds)
npcs <- min(opt$npcs, ncol(pc$vectors))
V <- pc$vectors[, seq_len(npcs), drop = FALSE]
colnames(V) <- paste0("PC", seq_len(npcs))

pcs_out <- data.table(sample = rownames(pc$vectors), pc$vectors[, seq_len(min(20, ncol(pc$vectors)))])
setnames(pcs_out, c("sample", paste0("PC", seq_len(ncol(pcs_out) - 1))))
pcs_out[, set := fifelse(sample %in% panel_ids, "1kGP", "cohort")]

lab <- setNames(sp$superpop, sp$sample)[panel_ids]
keep_ref <- panel_ids[lab %in% groups]
Vref <- V[keep_ref, , drop = FALSE]; lref <- lab[keep_ref]
cent <- t(sapply(split(seq_along(lref), lref), function(i) colMeans(Vref[i, , drop = FALSE])))
coh <- intersect(cohort_ids, rownames(V))
Vc <- V[coh, , drop = FALSE]
d_cent <- as.matrix(dist(rbind(cent, Vc)))[rownames(Vc), rownames(cent), drop = FALSE]
nearest <- colnames(d_cent)[apply(d_cent, 1, which.min)]
# Posterior from a Gaussian classifier (QDA) fitted to the panel PCs, equal priors
# (in the spirit of gnomAD's ancestry probabilities). Admixed groups (AMR) have broad
# covariance, so genuine members of tight clusters keep high posteriors.
fit <- MASS::qda(Vref, grouping = factor(lref), prior = rep(1 / nrow(cent), nrow(cent)))
post <- predict(fit, Vc)$posterior
conf <- post[cbind(seq_along(coh), match(nearest, colnames(post)))]
anc <- data.table(genotype_id = coh, nearest_superpop = nearest, posterior = round(conf, 3))
anc[, ancestry := fifelse(posterior >= opt$prob_min, nearest_superpop, "unassigned")]
anc <- cbind(anc, round(Vc[, seq_len(min(5, npcs)), drop = FALSE], 5))

# Per RNA-seq sample (replicates inherit their genotype's ancestry)
m <- merge(meta[, .(rnaseq_id, genotype_id)], anc, by = "genotype_id", all.x = TRUE, sort = FALSE)
m[is.na(ancestry), ancestry := "no_genotypes"]
fwrite(m[order(rnaseq_id)], paste0(opt$out, "_ancestry.tsv"), sep = "\t", quote = FALSE, na = "NA")

st <- m[, .(n_rnaseq = .N, n_individuals = uniqueN(genotype_id)), by = ancestry]
st[, status := fifelse(ancestry %in% c("unassigned", "no_genotypes"), ancestry,
                       fifelse(n_individuals >= opt$min_n, "mapped", "below_min_n"))]
setorder(st, -n_individuals)
fwrite(st[, .(ancestry, n_individuals, n_rnaseq, status)], paste0(opt$out, "_strata.tsv"), sep = "\t", quote = FALSE)
for (a in st[status == "mapped", ancestry]) {
  fwrite(m[ancestry == a, .(rnaseq_id, genotype_id)][order(rnaseq_id)],
         paste0("stratum_", a, ".samples.tsv"), sep = "\t", quote = FALSE)
}
fwrite(pcs_out[set == "1kGP"], paste0(opt$out, "_pcair_pcs.tsv"), sep = "\t", quote = FALSE)
# Cohort PCs are written separately (sample-level; never exported)
fwrite(pcs_out[set == "cohort"], paste0(opt$out, "_pcair_pcs_cohort.tsv"), sep = "\t", quote = FALSE)

pd <- rbind(data.table(PC1 = V[keep_ref, 1], PC2 = V[keep_ref, 2], group = paste0("1kGP ", lref), set = "1kGP"),
            data.table(PC1 = Vc[, 1], PC2 = Vc[, 2], group = anc$ancestry, set = "cohort"))
p <- ggplot(pd, aes(PC1, PC2)) +
  geom_point(data = pd[set == "1kGP"], aes(colour = group), alpha = 0.25, size = 1) +
  geom_point(data = pd[set == "cohort"], aes(fill = group), shape = 21, size = 1.8, colour = "black", stroke = 0.2) +
  scale_colour_brewer(palette = "Set2", name = "Reference") +
  scale_fill_brewer(palette = "Dark2", name = "Cohort assignment") +
  labs(title = "Genetic ancestry (PC-AiR)") + theme_bw(base_size = 10) +
  theme(text = element_text(family = "Liberation Sans"))
ggsave(paste0(opt$out, "_pca.png"), p, width = 7, height = 5, dpi = 150)
print(as.data.frame(st))
