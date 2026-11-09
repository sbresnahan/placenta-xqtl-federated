#!/usr/bin/env bash
# make_refs_fixture.sh -- rebuild the chr22 phenotyping reference fixture (tests/fixtures/refs_chr22).
#
# STAND-IN ANNOTATION: the committed fixture uses GENCODE v45 chr22, not HPLRv2 (HPLRv2 is public
# only through the MDACC reference record). Genes of type protein_coding or lncRNA and their
# transcripts tagged "basic" are kept and written in the normalized form the atlas uses
# (prep_PANTRY_gtf.R output: gene lines, gene_biotype / transcript_biotype / gene_name attributes).
# txrevise events are then built from it exactly as in the atlas reference prep (PANTRY txrevise
# scripts: extractTranscriptTags.py -> prepareAnnotations.R -> constructEvents.R, one batch -> merge).
# PANTRY's preprocess_gtf.py is not needed: the normalized GTF already carries the biotype fields.
#
# Requirements (development only, not in the image): python3, R with optparse, dplyr, purrr,
# rtracklayer, GenomicFeatures and txrevise (github.com/kauralasoo/txrevise @ badb87f) on
# R_LIBS / XQTL_TXREVISE_RLIB, and a PANTRY checkout (github.com/PejLab/Pantry @ 8da7787).
#
# Usage: tests/fixtures/make_refs_fixture.sh [GENCODE_GTF_GZ] [PANTRY_DIR]
#   GENCODE_GTF_GZ defaults to a download of gencode.v45.annotation.gtf.gz (sha256-checked).
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT=${HERE}/refs_chr22
GENCODE_URL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
GENCODE_SHA256=56a166e4a43c934a1910c0d09dec786ac048862fc65130731d41f3c7f451ae68
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

GTF_GZ=${1:-}
PANTRY=${2:-${TMP}/pantry}
if [ -z "$GTF_GZ" ]; then
    GTF_GZ=${TMP}/gencode.v45.annotation.gtf.gz
    curl -fsSL -o "$GTF_GZ" "$GENCODE_URL"
fi
got=$(sha256sum "$GTF_GZ" | cut -d' ' -f1)
if [ "$got" != "$GENCODE_SHA256" ]; then
    echo "ERROR: $GTF_GZ sha256 $got != pinned $GENCODE_SHA256" >&2
    exit 1
fi
if [ ! -d "$PANTRY/phenotyping/scripts/txrevise" ]; then
    git clone -q https://github.com/PejLab/Pantry.git "$PANTRY"
    git -C "$PANTRY" checkout -q 8da77879ba7a6c12c5f7a2e95954ae144f68130c
fi
if [ -n "${XQTL_TXREVISE_RLIB:-}" ]; then export R_LIBS="${XQTL_TXREVISE_RLIB}${R_LIBS:+:$R_LIBS}"; fi

mkdir -p "$OUT/annotation" "$OUT/txrevise"
GTF=${TMP}/gencode.v45.chr22.PANTRY.gtf

# ---- 1. chr22 normalized GTF ---------------------------------------------------------------
python3 - "$GTF_GZ" "$GTF" <<'EOF'
import gzip, re, sys

src, out = sys.argv[1], sys.argv[2]
KEEP_TYPES = {"protein_coding", "lncRNA"}
KEEP_ATTR = ["gene_id", "transcript_id", "gene_name", "gene_biotype", "transcript_biotype",
             "transcript_name", "exon_number", "tag"]
attr_re = re.compile(r'(\S+) ("[^"]*"|\S+);')

def parse(a):
    return [(k, v.strip('"')) for k, v in attr_re.findall(a)]

def fmt(pairs):
    out = []
    for key in KEEP_ATTR:
        for k, v in pairs:
            if k == key:
                out.append(f'{k} "{v}";' if k != "exon_number" else f"{k} {v};")
    return " ".join(out)

rows = []
with gzip.open(src, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if c[0] != "chr22" or c[2] not in ("gene", "transcript", "exon", "CDS"):
            continue
        a = [("gene_biotype" if k == "gene_type" else "transcript_biotype" if k == "transcript_type" else k, v)
             for k, v in parse(c[8])]
        d = dict(a)
        if d.get("gene_biotype") not in KEEP_TYPES:
            continue
        tags = {v for k, v in a if k == "tag"}
        if c[2] != "gene" and "basic" not in tags:
            continue
        if c[2] == "gene":   # gene lines carry only gene-level attributes, as in prep_PANTRY_gtf.R
            a = [(k, v) for k, v in a if k in ("gene_id", "gene_biotype", "gene_name")]
        rows.append((c[:8], a, d["gene_id"]))

genes_with_tx = {g for c, a, g in rows if c[2] == "transcript"}
with open(out, "w") as o:
    o.write("##gff-version 2\n##source: GENCODE v45 chr22 (protein_coding + lncRNA, basic transcripts), "
            "normalized like prep_PANTRY_gtf.R -- test fixture stand-in for HPLRv2\n")
    for c, a, g in rows:
        if g in genes_with_tx:
            o.write("\t".join(c) + "\t" + fmt(a) + "\n")
EOF
gzip -9 -n -c "$GTF" > "$OUT/annotation/gencode.v45.chr22.PANTRY.gtf.gz"

# ---- 2. txrevise events (PANTRY scripts, as in atlas 001/002/003_refprep) -----------------
TX=${TMP}/txrevise
mkdir -p "$TX/batch"
python3 "$PANTRY/phenotyping/scripts/txrevise/extractTranscriptTags.py" \
    --gtf <(gzip -c "$GTF") > "$TX/transcript_tags.txt"
Rscript "$PANTRY/phenotyping/scripts/txrevise/prepareAnnotations.R" \
    --gtf "$GTF" --tags "$TX/transcript_tags.txt" --out "$TX/txrevise_annotations.rds"
# constructEvents.R is single-threaded (~2 s per gene): run NBATCH gene batches in parallel, then
# concatenate them in batch order (as the atlas merge of its batch array does)
NBATCH=${NBATCH:-$(nproc)}
seq 1 "$NBATCH" | xargs -P "$NBATCH" -I{} Rscript "$PANTRY/phenotyping/scripts/txrevise/constructEvents.R" \
    --annot "$TX/txrevise_annotations.rds" --batch "{} $NBATCH" --out "$TX/batch" > "$TX/construct.log" 2>&1 \
    || { tail -30 "$TX/construct.log" >&2; exit 1; }
for g in grp_1 grp_2; do
    for p in upstream contained downstream; do
        for i in $(seq 1 "$NBATCH"); do
            f="$TX/batch/txrevise.${g}.${p}.${i}_${NBATCH}.gff3"
            [ -f "$f" ] || { echo "missing txrevise batch output $f" >&2; tail -30 "$TX/construct.log" >&2; exit 1; }
            grep -v '^#' "$f" || true
        done | gzip -9 -n > "$OUT/txrevise/txrevise.${g}.${p}.gff3.gz"
    done
done

# ---- 3. README + checksums ------------------------------------------------------------------
cat > "$OUT/README.txt" <<'EOF'
chr22 phenotyping reference fixture (tests only; rebuilt by tests/fixtures/make_refs_fixture.sh).
annotation/  GENCODE v45 chr22 stand-in for HPLRv2, normalized like the atlas (prep_PANTRY_gtf.R form)
txrevise/    txrevise events from that annotation (PANTRY txrevise scripts, atlas settings)
editing/     REDIportal hg38 sites on chr22 (BED6), cut from the MDACC reference record
music/       mini MuSiC reference: 254 chr22 genes, <= 30 cells per cell type x donor (3,857 cells)
             editing/ and music/ are made by tests/fixtures/extract_atlas_refs_chr22.sh (seadragon)
The genome is not stored here: the tests cut chr22 out of the image's GRCh38 FASTA.
EOF
( cd "$OUT" && find . -type f ! -name MD5SUMS | sort | xargs md5sum > MD5SUMS )
du -sh "$OUT"
