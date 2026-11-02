#!/usr/bin/env bash
# extract_atlas_refs_chr22.sh -- MDACC-side helper (run once on seadragon) that collects the chr22
# pieces of the atlas phenotyping references, plus the atlas SNUH (cohort3) chr22 phenotypes of the
# test samples, into one archive (atlas_chr22_refs.tar) for M4 fixture building and validation.
#
#   refs/   chr22 subsets of the atlas reference_dir:
#             HPLRv2.0.annotated.PANTRY.chr22.gtf.gz     normalized HPLRv2 GTF
#             txrevise.{grp_1,grp_2}.{upstream,contained,downstream}.chr22.gff3.gz
#             rediportal_hg38.chr22.bed.gz               preprocessed REDIportal sites
#             TABLE1_hg38_v3.chr22.head2000.txt.gz       raw REDIportal rows (unit test of preprocessing)
#             edit_sites_to_genes.chr22.tsv.gz, exonic.chr22.gtf.gz, intronic.chr22.gtf.gz,
#             HPLRv2.chr22.refFlat.gz                    atlas-derived files, to check the ports
#             placenta_music_reference.mini.rds          MuSiC SCE: chr22 gene symbols, <= 30 cells
#                                                        per cell type x pseudo-subject (seed 42)
#   atlas/  cohort3 output/*.bed.gz rows on chr22, only the test-sample columns, phenotype_groups
#           (chr22 phenotypes), and MuSiC proportions of the test samples -- used to check that the
#           federated phenotyping reproduces the atlas; never committed to the public repository.
#
# Usage:
#   eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
#   conda activate samtools-1.16.1          # for bgzip/tabix
#   bash extract_atlas_refs_chr22.sh [-o OUT_DIR] [ACC ...]
# Default accessions: the 8 runs in snuh_chr22_fastq.tar. R runs in the atlas RStudio 4.3.1
# container (override with RSCRIPT="..."). Missing inputs are reported and skipped.
set -euo pipefail

BASE=${BASE:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY}
REFDIR=${BASE}/reference
COHORT_OUT=${BASE}/cohort3/output
DECONV=${BASE}/deconvolution
REDI_RAW=${REDI_RAW:-/rsrch5/home/epi/bhattacharya_lab/data/REDIportal/TABLE1_hg38_v3.txt.gz}
RSCRIPT=${RSCRIPT:-singularity exec --bind /rsrch5 --bind /rsrch9 /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript}
export R_LIBS_USER=${R_LIBS_USER:-/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1}
OUT_DIR=${PWD}
DEFAULT_ACCS=(SRR18553497 SRR18553498 SRR18553501 SRR18553511 SRR18553512 SRR18553513 SRR18553520 SRR18553521)

while getopts "o:h" opt; do
    case $opt in
        o) OUT_DIR=$OPTARG ;;
        h) sed -n '2,27p' "$0"; exit 0 ;;
        *) exit 2 ;;
    esac
done
shift $((OPTIND - 1))
if [ "$#" -gt 0 ]; then ACCS=("$@"); else ACCS=("${DEFAULT_ACCS[@]}"); fi
command -v bgzip > /dev/null || { echo "bgzip not on PATH (conda activate samtools-1.16.1)" >&2; exit 1; }

STAGE=${OUT_DIR}/atlas_chr22
mkdir -p "$STAGE/refs" "$STAGE/atlas"
MISSING=()
have() { if [ -s "$1" ]; then return 0; fi; echo "WARNING: missing $1 -- skipped" >&2; MISSING+=("$1"); return 1; }
cat_any() { case $1 in *.gz) zcat "$1" ;; *) cat "$1" ;; esac; }
first() { for f in "$@"; do if [ -s "$f" ]; then echo "$f"; return 0; fi; done; echo "$1"; }

echo "[$(date +%T)] reference subsets"
# GTF-like files: keep header comments and chr22 rows
for pair in "HPLRv2.0.annotated.PANTRY.gtf:HPLRv2.0.annotated.PANTRY.chr22.gtf.gz" \
            "exonic.gtf:exonic.chr22.gtf.gz" "intronic.gtf:intronic.chr22.gtf.gz"; do
    src=$(first "${REFDIR}/${pair%%:*}" "${REFDIR}/${pair%%:*}.gz")
    if have "$src"; then
        cat_any "$src" | awk -F'\t' '/^#/ || $1 == "chr22"' | gzip -9 > "$STAGE/refs/${pair##*:}"
    fi
done
for g in grp_1 grp_2; do
    for p in upstream contained downstream; do
        src=$(first "${REFDIR}/txrevise/txrevise.${g}.${p}.gff3" "${REFDIR}/txrevise/txrevise.${g}.${p}.gff3.gz")
        if have "$src"; then
            cat_any "$src" | awk -F'\t' '/^#/ || $1 == "chr22"' | gzip -9 > "$STAGE/refs/txrevise.${g}.${p}.chr22.gff3.gz"
        fi
    done
done
src=$(first "${REFDIR}/rediportal_hg38.bed" "${REFDIR}/rediportal_hg38.bed.gz")
if have "$src"; then
    cat_any "$src" | awk -F'\t' '/^#/ || $1 == "chr22"' | gzip -9 > "$STAGE/refs/rediportal_hg38.chr22.bed.gz"
fi
src=$(first "${REFDIR}/edit_sites_to_genes.tsv" "${REFDIR}/edit_sites_to_genes.tsv.gz")
if have "$src"; then
    # keep the header and rows mentioning chr22 in any column (site IDs carry the contig)
    cat_any "$src" | awk 'NR == 1 || /chr22[^0-9]/' | gzip -9 > "$STAGE/refs/edit_sites_to_genes.chr22.tsv.gz"
fi
src=$(first "${REFDIR}/HPLRv2.refFlat" "${REFDIR}/HPLRv2.refFlat.gz")
if have "$src"; then
    cat_any "$src" | awk -F'\t' '$3 == "chr22"' | gzip -9 > "$STAGE/refs/HPLRv2.chr22.refFlat.gz"
fi
if have "$REDI_RAW"; then
    # header + first 2000 chr22 rows (pipefail-safe: awk reads to the end)
    zcat "$REDI_RAW" | awk -F'\t' 'NR == 1 { print; next } $2 == "chr22" && n < 2000 { print; n++ }' \
        | gzip -9 > "$STAGE/refs/TABLE1_hg38_v3.chr22.head2000.txt.gz"
fi

echo "[$(date +%T)] mini MuSiC reference"
if have "${DECONV}/placenta_music_reference.rds" && [ -s "$STAGE/refs/HPLRv2.0.annotated.PANTRY.chr22.gtf.gz" ]; then
    zcat "$STAGE/refs/HPLRv2.0.annotated.PANTRY.chr22.gtf.gz" \
        | awk -F'\t' '$3 == "gene"' | grep -o 'gene_name "[^"]*"' | cut -d'"' -f2 | sort -u > "$STAGE/refs/chr22_gene_names.txt"
    $RSCRIPT - "${DECONV}/placenta_music_reference.rds" "$STAGE/refs/chr22_gene_names.txt" \
        "$STAGE/refs/placenta_music_reference.mini.rds" <<'EOF'
suppressPackageStartupMessages({ library(SingleCellExperiment); library(Matrix) })
a <- commandArgs(trailingOnly = TRUE)
sce <- readRDS(a[1])
genes <- readLines(a[2])
set.seed(42)
grp <- paste(sce$cell.type, sce$subject, sep = "|")
idx <- unlist(lapply(split(seq_len(ncol(sce)), grp), function(i) if (length(i) > 30) sample(i, 30) else i))
mini <- sce[rownames(sce) %in% genes, sort(idx)]
mini$cell.type <- droplevels(mini$cell.type)
mini$subject <- droplevels(mini$subject)
saveRDS(mini, a[3], compress = "xz")
cat("mini MuSiC reference:", nrow(mini), "genes x", ncol(mini), "cells,", nlevels(mini$cell.type), "cell types\n")
EOF
fi

echo "[$(date +%T)] atlas cohort3 chr22 phenotypes for: ${ACCS[*]}"
ACC_LIST=$(IFS=,; echo "${ACCS[*]}")
for bed in expression isoforms isoform_expression alt_TSS alt_polyA splicing intron_retention RNA_editing stability; do
    src=${COHORT_OUT}/${bed}.bed.gz
    have "$src" || continue
    # header row starts with '#chr'; keep the 4 coordinate/ID columns plus the test samples
    zcat "$src" | awk -F'\t' -v OFS='\t' -v accs="$ACC_LIST" '
        NR == 1 {
            n = split(accs, want, ","); for (k = 1; k <= n; k++) w[want[k]] = 1
            for (i = 1; i <= NF; i++) if (i <= 4 || ($i in w)) keep[++m] = i
            if (m - 4 != n) { print "ERROR: " FILENAME " has " m - 4 " of " n " test samples" > "/dev/stderr"; exit 1 }
        }
        NR == 1 || $1 == "chr22" { line = $keep[1]; for (j = 2; j <= m; j++) line = line OFS $keep[j]; print line }
    ' | bgzip > "$STAGE/atlas/${bed}.chr22.bed.gz"
    # phenotype groups (phenotype_id <tab> group) of the chr22 phenotypes
    grp=${COHORT_OUT}/${bed}.phenotype_groups.txt
    if [ -s "$grp" ]; then
        awk -F'\t' 'NR == FNR { if (FNR > 1) id[$4] = 1; next } ($1 in id)' \
            <(zcat "$STAGE/atlas/${bed}.chr22.bed.gz") "$grp" | gzip -9 > "$STAGE/atlas/${bed}.phenotype_groups.chr22.txt.gz"
    fi
done
for f in "${DECONV}/cohort3_cell_proportions_full.tsv" "${DECONV}/cohort3_cell_proportions_collapsed.tsv"; do
    have "$f" || continue
    awk -F'\t' -v accs="$ACC_LIST" 'BEGIN { n = split(accs, a, ","); for (k = 1; k <= n; k++) w[a[k]] = 1 } NR == 1 || ($1 in w)' "$f" \
        > "$STAGE/atlas/$(basename "$f")"
done

{
    echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "host: $(hostname)"
    echo "samples: ${ACCS[*]}"
    echo "missing inputs: ${MISSING[*]:-none}"
} > "$STAGE/versions.txt"
( cd "$STAGE" && find . -type f ! -name MD5SUMS -print0 | sort -z | xargs -0 md5sum > MD5SUMS )
tar -C "$OUT_DIR" -cf "$OUT_DIR/atlas_chr22_refs.tar" atlas_chr22
echo "Wrote $OUT_DIR/atlas_chr22_refs.tar"
ls -lh "$STAGE/refs" "$STAGE/atlas"
[ "${#MISSING[@]}" -eq 0 ] || echo "NOTE: ${#MISSING[@]} input(s) missing; see versions.txt"
