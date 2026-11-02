#!/usr/bin/env bash
# extract_snuh_chr22.sh -- MDACC-side helper (run once on seadragon) that cuts the chr22 read
# pairs out of the atlas SNUH (cohort3, PRJNA820329) STAR BAMs and archives them as FASTQ for the
# phenotyping test fixture (tests/fixtures/fastq; see tests/fixtures/make_fastq_fixture.sh).
#
# What is kept: primary alignments (-F 0x900) of pairs whose two mates both align to chr22.
# Pairs with a mate elsewhere or unaligned are dropped (STAR writes no unmapped reads by default).
# The BAMs hold the atlas's Trim Galore-trimmed reads, so read lengths vary (<= 101 bp).
# samtools fastq restores reverse-strand reads to their sequenced orientation.
# The BAMs are only read: their index is written to the scratch directory, not next to the BAM.
#
# Usage (login node is fine; ~2-5 min per sample):
#   eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
#   conda activate samtools-1.16.1
#   bash extract_snuh_chr22.sh [-b BAM_DIR] [-o OUT_DIR] [-t THREADS] [ACC ...]
#
# Default accessions: the 8 runs below. Missing BAMs are reported and skipped.
# Output: OUT_DIR/snuh_chr22_fastq.tar containing snuh_chr22/<ACC>_R1.fastq.gz,
#         <ACC>_R2.fastq.gz, manifest.tsv and versions.txt. Upload the .tar.
set -euo pipefail

BAM_DIR=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/cohort3/intermediate/star_out
OUT_DIR=${PWD}
THREADS=4
REGION=chr22
DEFAULT_ACCS=(SRR18553497 SRR18553498 SRR18553501 SRR18553511 SRR18553512 SRR18553513 SRR18553520 SRR18553521)

while getopts "b:o:t:h" opt; do
    case $opt in
        b) BAM_DIR=$OPTARG ;;
        o) OUT_DIR=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        h) sed -n '2,20p' "$0"; exit 0 ;;
        *) exit 2 ;;
    esac
done
shift $((OPTIND - 1))
if [ "$#" -gt 0 ]; then ACCS=("$@"); else ACCS=("${DEFAULT_ACCS[@]}"); fi

command -v samtools > /dev/null || { echo "samtools not on PATH (conda activate samtools-1.16.1)" >&2; exit 1; }

STAGE=${OUT_DIR}/snuh_chr22
TMP=$(mktemp -d "${OUT_DIR}/.extract_tmp.XXXXXX")
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$STAGE"

printf 'sample\tbam\tbam_bytes\tbam_mtime\tpairs_written\tR1_md5\tR2_md5\n' > "$STAGE/manifest.tsv"
{
    echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "host: $(hostname)"
    echo "region: ${REGION}; filter: -F 0x900, both mates on ${REGION}"
    samtools --version | head -2
} > "$STAGE/versions.txt"

n_ok=0
for acc in "${ACCS[@]}"; do
    bam=${BAM_DIR}/${acc}.Aligned.sortedByCoord.out.bam
    if [ ! -s "$bam" ]; then
        echo "WARNING: $bam not found -- skipped" >&2
        continue
    fi
    echo "[$(date +%T)] $acc"
    # index into scratch (read-only use of the atlas BAM)
    samtools index -@ "$THREADS" -o "$TMP/${acc}.bai" "$bam"
    if [ "$(samtools view -H "$bam" | awk -v c="SN:${REGION}" '$1 == "@SQ" && $2 == c { n++ } END { print n + 0 }')" -ne 1 ]; then
        echo "ERROR: contig $REGION not in $bam header" >&2
        exit 1
    fi
    # primary alignments in chr22 whose mate is also on chr22 (RNEXT '=' ), then pair up by name
    samtools view -h -F 0x900 -X "$bam" "$TMP/${acc}.bai" "$REGION" \
        | awk -F'\t' '/^@/ || $7 == "="' \
        | samtools sort -n -@ "$THREADS" -m 1G -T "$TMP/${acc}.nsort" -o "$TMP/${acc}.byname.bam" -
    samtools fastq -@ "$THREADS" -n -c 9 \
        -1 "$STAGE/${acc}_R1.fastq.gz" -2 "$STAGE/${acc}_R2.fastq.gz" \
        -0 /dev/null -s /dev/null "$TMP/${acc}.byname.bam" 2> "$TMP/${acc}.fastq.log"
    pairs=$(( $(zcat "$STAGE/${acc}_R1.fastq.gz" | wc -l) / 4 ))
    pairs2=$(( $(zcat "$STAGE/${acc}_R2.fastq.gz" | wc -l) / 4 ))
    if [ "$pairs" -ne "$pairs2" ] || [ "$pairs" -eq 0 ]; then
        echo "ERROR: $acc R1/R2 read counts differ or are zero ($pairs vs $pairs2)" >&2
        exit 1
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$acc" "$(basename "$bam")" "$(stat -c %s "$bam")" \
        "$(stat -c %y "$bam" | cut -d. -f1)" "$pairs" \
        "$(md5sum < "$STAGE/${acc}_R1.fastq.gz" | cut -d' ' -f1)" \
        "$(md5sum < "$STAGE/${acc}_R2.fastq.gz" | cut -d' ' -f1)" >> "$STAGE/manifest.tsv"
    rm -f "$TMP/${acc}".*
    echo "    ${pairs} pairs"
    n_ok=$((n_ok + 1))
done

[ "$n_ok" -gt 0 ] || { echo "ERROR: no BAMs processed" >&2; exit 1; }
tar -C "$OUT_DIR" -cf "$OUT_DIR/snuh_chr22_fastq.tar" snuh_chr22
echo "Wrote $OUT_DIR/snuh_chr22_fastq.tar ($n_ok samples):"
cat "$STAGE/manifest.tsv"
ls -lh "$OUT_DIR/snuh_chr22_fastq.tar"
