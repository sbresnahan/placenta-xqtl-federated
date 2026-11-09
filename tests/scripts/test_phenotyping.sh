#!/usr/bin/env bash
# M4 acceptance test: FASTQ -> phenotype BEDs on the chr22 SNUH fixture (8 public placenta samples,
# PRJNA820329, reads cut to chr22; docs/PHENOTYPING.md#test-fixture).
#   1. -profile test_phenotyping,<engine> --mode prepare_references builds the site indices
#      (STAR, Salmon expression + 4 txrevise sets, IRFinder-S, Picard references) for chr22;
#   2. --mode phenotyping_only produces every modality BED, the QC table and the normalized
#      matrices; check_phenotypes.py checks shapes, TSS positions, groups and the frozen /
#      novel splicing clusters;
#   3. the validator refuses phenotyping without site_refs.
# rna_editing is checked when tests/fixtures/refs_chr22/editing/rediportal_hg38.bed.gz exists, and cell
# proportions when tests/fixtures/refs_chr22/music/placenta_music_reference.rds exists (both are committed).
#
#   tests/scripts/test_phenotyping.sh [WORKDIR] [ENGINE]     ENGINE: dev (default; tools on PATH) | docker
# The chr22 genome is cut from XQTL_GENOME (default: the image's GRCh38 FASTA) unless XQTL_CHR22_FA
# points to a ready chr22 FASTA. Extra Nextflow options can be passed in NXF_OPTS_EXTRA.
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
work=${1:-$(mktemp -d)}
engine=${2:-dev}
mkdir -p "$work"; work=$(cd "$work" && pwd)
cd "$work"
pass=0 fail=0
ok()  { echo "PASS  $*"; pass=$((pass + 1)); }
bad() { echo "FAIL  $*"; fail=$((fail + 1)); }
out=results_test_phenotyping
refs=$work/site_refs
fx=$repo/tests/fixtures
# shellcheck disable=SC2086
nxf() { nextflow run "$repo/main.nf" -profile "test_phenotyping,$engine" ${NXF_OPTS_EXTRA:-} "$@"; }

# ---- chr22 genome -------------------------------------------------------------------------
if [ -n "${XQTL_CHR22_FA:-}" ]; then
  fa=$XQTL_CHR22_FA
else
  fa=$work/chr22.fa
  if [ ! -s "$fa" ]; then
    genome=${XQTL_GENOME:-/opt/xqtl/refs/genome/GRCh38_no_alt_analysis_set.fa.gz}
    if [ -s "$genome.fai" ]; then samtools faidx "$genome" chr22 > "$fa"
    else gzip -cdf "$genome" | awk '/^>/{p = ($1 == ">chr22")} p' > "$fa"; fi
  fi
fi
[ -s "$fa" ] || { echo "no chr22 FASTA (set XQTL_GENOME or XQTL_CHR22_FA)"; exit 1; }

# ---- 1. site indices ----------------------------------------------------------------------
if nxf --mode prepare_references --genome_fasta "$fa" --site_refs "$refs" -resume > prepare.log 2>&1; then
  ok "--mode prepare_references completed"
else
  bad "--mode prepare_references failed (see $work/prepare.log)"; tail -40 prepare.log; exit 1
fi
n_idx=$(find "$refs/salmon" -mindepth 2 -maxdepth 2 -name info.json | wc -l)
if grep -qP '^read_length\t101$' "$refs/MANIFEST.tsv" && [ "$n_idx" -eq 5 ] && [ -s "$refs/star/SA" ] \
   && [ -s "$refs/irfinder/IRFinder/ref-cover.bed" ] && [ -s "$refs/annotation/refFlat.txt" ]; then
  ok "site_refs: STAR, IRFinder-S, Picard references and 5 Salmon indices; MANIFEST records read length 101"
else
  bad "site_refs contents ($n_idx Salmon indices)"
fi

# ---- 2. phenotyping -----------------------------------------------------------------------
if nxf --mode phenotyping_only --site_refs "$refs" -resume > phenotyping.log 2>&1; then
  ok "--mode phenotyping_only completed"
else
  bad "--mode phenotyping_only failed (see $work/phenotyping.log)"; tail -40 phenotyping.log; exit 1
fi
mods=expression,isoform_expression,isoform_usage,splicing,intron_retention,alt_tss,alt_polya,stability
[ -s "$fx/refs_chr22/editing/rediportal_hg38.bed.gz" ] && mods=$mods,rna_editing
if python3 "$repo/tests/scripts/check_phenotypes.py" "$out" --gtf "$fx/refs_chr22/annotation/gencode.v45.chr22.PANTRY.gtf.gz" \
     --metadata "$fx/phenotyping/metadata.tsv" --restricted "$fx/phenotyping/xqtl-restricted-phenotyping-test.tar.gz" \
     --modalities "$mods" > check.log; then
  ok "phenotype BEDs and normalized matrices ($mods)"
else
  bad "phenotype checks (see $work/check.log)"
fi
cat check.log
nq=$(($(wc -l < "$out/phenotypes/qc_metrics.tsv") - 1))
if [ "$nq" -eq 8 ] && head -1 "$out/phenotypes/qc_metrics.tsv" | grep -q 'SubjectBias.GC'; then
  ok "qc_metrics.tsv: 8 samples with Picard metrics and GC / gene-length bias"
else
  bad "qc_metrics.tsv ($nq rows)"
fi
if [ -s "$fx/refs_chr22/music/placenta_music_reference.rds" ]; then
  if [ "$(($(wc -l < "$out/phenotypes/cell_proportions.tsv") - 1))" -eq 8 ]; then ok "MuSiC cell proportions for 8 samples"
  else bad "cell_proportions.tsv"; fi
fi

# ---- 3. validator ---------------------------------------------------------------------------
if nxf --mode phenotyping_only --site_refs "$work/no_such_dir" > noref.log 2>&1; then
  bad "phenotyping with a missing site_refs directory passed validation"
elif grep -q 'site_refs' noref.log; then
  ok "missing site_refs is reported by the validator"
else
  bad "missing site_refs failed, but not in validation (see noref.log)"
fi

echo "test_phenotyping.sh: $pass passed, $fail failed (work dir $work)"
[ "$fail" -eq 0 ]
