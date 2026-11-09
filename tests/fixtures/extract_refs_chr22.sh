#!/usr/bin/env bash
# extract_refs_chr22.sh -- build the editing/ and music/ parts of the chr22 phenotyping reference
# fixture (tests/fixtures/refs_chr22) from the public reference record.
#
#   editing/rediportal_hg38.bed.gz         chr22 rows of rediportal_hg38.bed.gz (BED6)
#   music/placenta_music_reference.rds     mini MuSiC reference: genes on chr22 of the HPLRv2
#                                          normalized GTF, <= 30 cells per cell type x donor (seed 42)
#
# Usage:
#   tests/fixtures/extract_refs_chr22.sh RECORD_DIR [OUT_DIR]
# RECORD_DIR holds rediportal_hg38.bed.gz, placenta_music_reference.rds and
# HPLRv2.0.annotated.PANTRY.gtf.gz (a download of the reference record, or the image's refs
# directory). OUT_DIR defaults to tests/fixtures/refs_chr22. Requires gzip and Rscript with
# SingleCellExperiment and Matrix (scripts/dev_env.sh provides both); RSCRIPT overrides Rscript.
# Afterwards re-run the MD5SUMS step of tests/fixtures/make_refs_fixture.sh.
set -euo pipefail

REC=${1:?usage: extract_refs_chr22.sh RECORD_DIR [OUT_DIR]}
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${2:-${HERE}/refs_chr22}
RSCRIPT=${RSCRIPT:-Rscript}
locate() { for f in "$REC/$1" "$REC"/*/"$1"; do [ -s "$f" ] && { echo "$f"; return 0; }; done
           echo "missing $1 in $REC" >&2; exit 1; }
REDI=$(locate rediportal_hg38.bed.gz)
MUSIC=$(locate placenta_music_reference.rds)
GTF=$(locate HPLRv2.0.annotated.PANTRY.gtf.gz)
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
mkdir -p "$OUT/editing" "$OUT/music"

zcat "$REDI" | awk -F'\t' '/^#/ || $1 == "chr22"' | gzip -9 -n > "$OUT/editing/rediportal_hg38.bed.gz"
echo "editing: $(zcat "$OUT/editing/rediportal_hg38.bed.gz" | grep -vc '^#') chr22 sites"

zcat "$GTF" | awk -F'\t' '$1 == "chr22" && $3 == "gene"' | grep -o 'gene_name "[^"]*"' | cut -d'"' -f2 \
    | sort -u > "$TMP/chr22_gene_names.txt"
$RSCRIPT - "$MUSIC" "$TMP/chr22_gene_names.txt" "$OUT/music/placenta_music_reference.rds" <<'EOF'
a <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages({ library(SingleCellExperiment); library(Matrix) })
sce <- readRDS(a[1])
genes <- readLines(a[2])
set.seed(42)
grp <- paste(sce$cell.type, sce$subject, sep = "|")
idx <- unlist(lapply(split(seq_len(ncol(sce)), grp), function(i) if (length(i) > 30) sample(i, 30) else i))
mini <- sce[rownames(sce) %in% genes, sort(idx)]
for (v in c("cell.type", "subject")) if (is.factor(mini[[v]])) mini[[v]] <- droplevels(mini[[v]])
saveRDS(mini, a[3], compress = "xz")
cat("music: mini reference with", nrow(mini), "genes x", ncol(mini), "cells,",
    length(unique(mini$cell.type)), "cell types\n")
EOF
