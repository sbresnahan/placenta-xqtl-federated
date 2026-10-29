#!/usr/bin/env bash
# MDACC only, run once per reference release on seadragon: package the public reference
# files derived from the atlas for upload to a public Zenodo record.
#
#   scripts/make_public_reference_record.sh --env scripts/packaging/seadragon.env --out DIR [--dry-run] [--threads N]
#
# Writes DIR/ (flat, ready to upload as one Zenodo record):
#   HPLRv2.0.annotated.gtf.gz, HPLRv2.0.annotated.PANTRY.gtf.gz, HPLRv2_classification.annotated.txt.gz
#   txrevise.{grp_1,grp_2}.{upstream,contained,downstream}.gff3.gz
#   placenta_music_reference.rds
#   1kGP_hg38_pca_panel.{pgen,pvar.zst,psam}, 1kg_superpop.tsv
#   rediportal_hg38.bed.gz                (only with INCLUDE_REDIPORTAL=yes)
#   atlas_genome_contigs.tsv              (name, length, MD5 of the atlas genome FASTA)
#   files.tsv, SHA256SUMS, README.md, build_log.txt
# Then: upload everything in DIR to Zenodo (public), publish, and pin the record with
#   scripts/pin_reference_record.py --record <id> --files DIR/files.tsv
# Tools: gzip, sha256sum, plink2 (>= 2.0.0a6), samtools. On seadragon: module load plink samtools
# Compression uses gzip -n so the outputs (and checksums) do not depend on file dates.
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
envf="" out="" dry=0 threads=8
while [ $# -gt 0 ]; do
  case "$1" in
    --env) envf=$2; shift 2 ;;
    --out) out=$2; shift 2 ;;
    --dry-run) dry=1; shift ;;
    --threads) threads=$2; shift 2 ;;
    -h|--help) sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown argument $1" >&2; exit 2 ;;
  esac
done
[ -n "$envf" ] && [ -n "$out" ] || { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//' >&2; exit 2; }
# shellcheck disable=SC1090
source "$envf"
: "${PANEL_MAF:=0.05}" "${PANEL_PRUNE:=1000 100 0.2}" "${INCLUDE_REDIPORTAL:=no}"

log() { printf '[public_record] %s\n' "$*"; }
size() { if [ -e "$1" ]; then du -sh "$1" 2>/dev/null | cut -f1; else echo MISSING; fi; }

# ---- inputs ---------------------------------------------------------------------------
declare -a inputs=(
  "ATLAS_HPLR_GTF|$ATLAS_HPLR_GTF"
  "ATLAS_NORMALIZED_GTF|$ATLAS_NORMALIZED_GTF"
  "ATLAS_CLASSIFICATION|$ATLAS_CLASSIFICATION"
  "ATLAS_MUSIC_RDS|$ATLAS_MUSIC_RDS"
  "ATLAS_GENOME_FASTA|$ATLAS_GENOME_FASTA"
  "KGP_PFILE.pgen|$KGP_PFILE.pgen"
  "KGP_PFILE.pvar|$KGP_PFILE.pvar"
  "KGP_PFILE.psam|$KGP_PFILE.psam"
  "KGP_SUPERPOP|$KGP_SUPERPOP"
)
for g in grp_1 grp_2; do for p in upstream contained downstream; do
  inputs+=("txrevise.$g.$p|$ATLAS_TXREVISE_DIR/txrevise.$g.$p.gff3")
done; done
[ "$INCLUDE_REDIPORTAL" = yes ] && inputs+=("ATLAS_REDIPORTAL_BED|$ATLAS_REDIPORTAL_BED")

missing=0
log "inputs:"
for kv in "${inputs[@]}"; do
  k=${kv%%|*} v=${kv#*|}
  s=$(size "$v"); [ "$s" = MISSING ] && missing=$((missing + 1))
  printf '  %-22s %-8s %s\n' "$k" "$s" "$v"
done
log "PCA panel: biallelic autosomal SNPs, MAF >= $PANEL_MAF, --indep-pairwise $PANEL_PRUNE"
log "REDIportal sites included: $INCLUDE_REDIPORTAL"
tools_ok=1
for t in gzip sha256sum plink2 samtools; do
  if command -v "$t" >/dev/null; then log "tool  $t: $(command -v "$t")"; else log "tool  $t: NOT FOUND"; tools_ok=0; fi
done
if [ "$dry" -eq 1 ]; then
  log "dry run: would write to $out"
  log "outputs: HPLRv2 GTFs + classification (gz), 6 txrevise GFF3 (gz), MuSiC .rds, 1kGP panel"
  log "         (pgen/pvar.zst/psam + labels), atlas_genome_contigs.tsv$([ "$INCLUDE_REDIPORTAL" = yes ] && echo ', rediportal_hg38.bed.gz'),"
  log "         files.tsv, SHA256SUMS, README.md, build_log.txt"
  [ "$missing" -eq 0 ] && [ "$tools_ok" -eq 1 ] && log "dry run: READY" || log "dry run: NOT READY ($missing missing input(s), tools ok=$tools_ok)"
  exit 0
fi
[ "$missing" -eq 0 ] || { log "ERROR: $missing input(s) missing"; exit 1; }
[ "$tools_ok" -eq 1 ] || { log "ERROR: required tools missing"; exit 1; }

mkdir -p "$out"; out=$(cd "$out" && pwd)
tmp=$(mktemp -d "$out/.tmp.XXXX"); trap 'rm -rf "$tmp"' EXIT
exec > >(tee "$out/build_log.txt") 2>&1
log "started; pipeline commit $(git -C "$here" rev-parse --short HEAD 2>/dev/null || echo unknown)"

gz() { gzip -n -9 -c "$1" > "$2"; log "wrote $(basename "$2")"; }
gz "$ATLAS_HPLR_GTF" "$out/HPLRv2.0.annotated.gtf.gz"
gz "$ATLAS_NORMALIZED_GTF" "$out/HPLRv2.0.annotated.PANTRY.gtf.gz"
gz "$ATLAS_CLASSIFICATION" "$out/HPLRv2_classification.annotated.txt.gz"
for g in grp_1 grp_2; do for p in upstream contained downstream; do
  gz "$ATLAS_TXREVISE_DIR/txrevise.$g.$p.gff3" "$out/txrevise.$g.$p.gff3.gz"
done; done
cp "$ATLAS_MUSIC_RDS" "$out/placenta_music_reference.rds"; log "wrote placenta_music_reference.rds"
if [ "$INCLUDE_REDIPORTAL" = yes ]; then gz "$ATLAS_REDIPORTAL_BED" "$out/rediportal_hg38.bed.gz"; fi

# ---- 1kGP PCA panel ---------------------------------------------------------------------
# keep the labelled samples (one-column #IID file, matched against the psam IIDs)
iidcol=$(head -1 "$KGP_PFILE.psam" | awk '{print ($1 == "#FID") ? 2 : 1}')
awk -v c="$iidcol" 'NR == FNR { if (FNR > 1) ids[$c] = 1; next } ($1 in ids) { print $1 }' \
    "$KGP_PFILE.psam" "$KGP_SUPERPOP" | sort -u | { echo "#IID"; cat; } > "$tmp/keep.txt"
log "1kGP samples with superpopulation labels: $(($(wc -l < "$tmp/keep.txt") - 1))"
# shellcheck disable=SC2086
plink2 --pfile "$KGP_PFILE" --keep "$tmp/keep.txt" --autosome --snps-only just-acgt --max-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1000 --rm-dup exclude-all \
    --maf "$PANEL_MAF" --indep-pairwise $PANEL_PRUNE --threads "$threads" --out "$tmp/prune"
plink2 --pfile "$KGP_PFILE" --keep "$tmp/keep.txt" --autosome --snps-only just-acgt --max-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1000 --rm-dup exclude-all \
    --extract "$tmp/prune.prune.in" --sort-vars --make-pgen vzs pvar-cols=-xheader,-info \
    --threads "$threads" --out "$out/1kGP_hg38_pca_panel"
rm -f "$out/1kGP_hg38_pca_panel.log"
cp "$KGP_SUPERPOP" "$out/1kg_superpop.tsv"
log "PCA panel: $(zstd -dc "$out/1kGP_hg38_pca_panel.pvar.zst" 2>/dev/null | grep -vc '^#' || echo '?') variants"

# ---- genome provenance ------------------------------------------------------------------
samtools dict "$ATLAS_GENOME_FASTA" \
  | awk -F'\t' 'BEGIN { print "name\tlength\tmd5" }
                $1 == "@SQ" { for (i = 2; i <= NF; i++) { split($i, a, ":"); v[a[1]] = substr($i, 4) }
                              print v["SN"] "\t" v["LN"] "\t" v["M5"] }' > "$out/atlas_genome_contigs.tsv"
log "atlas genome: $(($(wc -l < "$out/atlas_genome_contigs.tsv") - 1)) contigs"

# ---- checksums and record description ----------------------------------------------------
( cd "$out"
  for f in *; do
    case "$f" in files.tsv|SHA256SUMS|README.md|build_log.txt) ;; *) [ -f "$f" ] && echo "$f" ;; esac
  done | LC_ALL=C sort > "$tmp/list"
  : > SHA256SUMS; printf 'name\tsha256\tbytes\n' > files.tsv
  while read -r f; do
    s=$(sha256sum "$f" | cut -d' ' -f1); b=$(stat -c %s "$f")
    printf '%s  %s\n' "$s" "$f" >> SHA256SUMS; printf '%s\t%s\t%s\n' "$f" "$s" "$b" >> files.tsv
  done < "$tmp/list" )
cat > "$out/README.md" <<EOF
# Public reference data for placenta-xqtl-federated

Reference files derived from the placental xQTL atlas (MD Anderson Cancer Center,
Bhattacharya Lab) and embedded in the placenta-xqtl-federated image. No individual-level data.

- HPLRv2 long-read placental transcriptome: SQANTI3-annotated GTF, normalized GTF, classification table
- txrevise event annotations (6 GFF3) built from HPLRv2
- MuSiC single-cell reference (placenta, Campbell et al. 2023, GEO GSE182381)
- 1000 Genomes high-coverage (hg38) PCA panel: biallelic autosomal SNPs, MAF >= $PANEL_MAF,
  LD-pruned (--indep-pairwise $PANEL_PRUNE), with superpopulation labels
$( [ "$INCLUDE_REDIPORTAL" = yes ] && echo "- REDIportal v3 hg38 A-to-I editing sites (BED6)" )
- atlas_genome_contigs.tsv: contig names, lengths and MD5s of the genome FASTA used by the atlas

Checksums: SHA256SUMS. Provenance and licenses: docs/REFERENCES.md in
https://github.com/sbresnahan/placenta-xqtl-federated
EOF
log "done: $(wc -l < "$out/SHA256SUMS") files in $out"
log "next: upload the contents of $out to a public Zenodo record, then"
log "      scripts/pin_reference_record.py --record <id> --files $out/files.tsv"
