#!/usr/bin/env bash
# M6 acceptance test: --mode central on two independent test cohorts, end to end.
#   1. cohort 1 = tests/fixtures/site (TESTCOHORT, parquet export) -- reused from a finished test_site.sh
#      work dir when one is given, otherwise run here;
#   2. cohort 2 = tests/fixtures/site2 (TESTCOHORT2, 120 other 1000 Genomes samples, the same planted effects,
#      other seeds) exported as tsv.gz, so the central ingest reads both formats; its planted effects are
#      checked with check_site_truth.py (nominal effect at every planted variant; feature q reported);
#   3. a synthetic GWAS (tests/fixtures/central/make_gwas.py) with one colocalizing and one distinct signal,
#      and a BED annotation of the planted variants (+-1 kb) for the PIP enrichment;
#   4. --mode central over both studies (study 2 given as a path relative to studies.tsv);
#   5. check_central.py: ingest, meta, features, fine-mapping, SuSiE-coloc + ColocBoost, G x C, GxE tiers 1-2,
#      tier-1 scope labels, TWAS inputs, summary;
#   6. check_central_ivw.R: FE / DL-RE IVW, Q, I^2 and tau^2 recomputed in R from the raw file A rows of both
#      exports; FUSION .wgt.RDat files load and have the expected shapes.
#
#   tests/scripts/test_central.sh [WORKDIR] [ENGINE] [SITE1_WORKDIR]     ENGINE: dev (default) | docker
# Extra Nextflow options can be passed in NXF_OPTS_EXTRA (e.g. "--image ghcr.io/...:dev").
# RSCRIPT selects the Rscript used for check_central_ivw.R (default: Rscript on PATH, needs data.table).
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
work=${1:-$(mktemp -d)}
engine=${2:-dev}
site1=${3:-}
mkdir -p "$work"; work=$(cd "$work" && pwd)
pass=0 fail=0
ok()  { echo "PASS  $*"; pass=$((pass + 1)); }
bad() { echo "FAIL  $*"; fail=$((fail + 1)); }
# shellcheck disable=SC2086
nxf() { nextflow run "$repo/main.nf" -profile "test,$engine" ${NXF_OPTS_EXTRA:-} "$@"; }
fx2=$repo/tests/fixtures/site2

# ---- 1. cohort 1 ------------------------------------------------------------------------------------
if [ -n "$site1" ] && [ -s "$site1/results_test/export/TESTCOHORT/MANIFEST.tsv" ]; then
  exp1=$(cd "$site1/results_test/export/TESTCOHORT" && pwd)
  ok "cohort 1: reusing the export of $site1"
else
  mkdir -p "$work/cohort1"; cd "$work/cohort1"
  if nxf --mode site -resume > site.log 2>&1; then ok "cohort 1: --mode site completed"; else bad "cohort 1: --mode site failed (see $work/cohort1/site.log)"; tail -40 site.log; exit 1; fi
  exp1=$work/cohort1/results_test/export/TESTCOHORT
fi

# ---- 2. cohort 2 (tsv.gz export) --------------------------------------------------------------------
mkdir -p "$work/cohort2"; cd "$work/cohort2"
cat > cohort2.yml <<EOF
cohort_id: TESTCOHORT2
outdir: results_test
samplesheet: $fx2/samplesheet.csv
metadata: $fx2/metadata.tsv
phenotypes_dir: $fx2/phenotypes
genotypes:
  raw: null
  raw_build: GRCh38
  source: wgs
  imputed_dir: $fx2/genotypes
  imputed_pattern: 'chr{chr}.dose.vcf.gz'
export:
  format: tsv.gz
EOF
if nxf -params-file cohort2.yml --mode site -resume > site.log 2>&1; then ok "cohort 2: --mode site completed"; else bad "cohort 2: --mode site failed (see $work/cohort2/site.log)"; tail -40 site.log; exit 1; fi
exp2=$work/cohort2/results_test/export/TESTCOHORT2
# per-stratum feature q at n = 60 is reported (WARN) but not required here: test_site.sh requires it for
# cohort 1, and the central feature test across the four contributions is checked below
if python3 "$repo/tests/scripts/check_site_truth.py" "$exp2" "$fx2/TRUTH.tsv" EUR EAS --feature-q-warn > truth.log; then
  ok "cohort 2: planted effects at the planted variants in EUR and EAS"
else
  bad "cohort 2: truth checks"
fi
cat truth.log
if [ -n "$(find "$exp2/EUR/A_nominal" -name '*.tsv.gz' -print -quit)" ] && [ -z "$(find "$exp2" -name '*.parquet' -print -quit)" ]; then
  ok "cohort 2 exported as tsv.gz"
else
  bad "cohort 2 export format"
fi

# ---- 3. central inputs ------------------------------------------------------------------------------
mkdir -p "$work/central/gwas" "$work/central/annot"; cd "$work/central"
# the distinct (H3) GWAS variant must be in the fine-mapping variant set shared by all four contributions
python3 "$repo/tests/fixtures/central/make_gwas.py" --export-dir "$exp1/EUR" --also-export-dir "$exp1/EAS" \
  --also-export-dir "$exp2/EUR" --also-export-dir "$exp2/EAS" --out gwas/gwas_EUR.tsv.gz > gwas.log
printf 'trait_id\tpath\tancestry\tn\ttrait_type\tprop_cases\nSYNTH_T1\tgwas/gwas_EUR.tsv.gz\tEUR\t20000\tquantitative\t\n' > gwas_catalog.tsv
awk -F'\t' 'NR > 1 {split($4, v, ":"); printf "%s\t%d\t%d\n", v[1], v[2] - 1000, v[2] + 1000}' "$fx2/TRUTH.tsv" \
  | sort -k1,1 -k2,2n -u > annot/planted_regions.bed
printf 'name\tbed\nplanted_regions\tannot/planted_regions.bed\n' > annotations.tsv
rel2=$(python3 -c 'import os, sys; print(os.path.relpath(sys.argv[1], sys.argv[2]))' "$exp2" "$work/central")
printf 'study_id\texport_dir\nS1\t%s\nS2\t%s\n' "$exp1" "$rel2" > studies.tsv
cat > central.yml <<EOF
cohort_id: MDACC_central
outdir: results_central
central:
  studies: $work/central/studies.tsv
  gwas_catalog: $work/central/gwas_catalog.tsv
  annotations: $work/central/annotations.tsv
  finemap:
    shard_size: 4
  coloc:
    shard_size: 4
EOF

# ---- 4. central run ---------------------------------------------------------------------------------
if nxf -params-file central.yml --mode central -resume > central.log 2>&1; then ok "--mode central completed"; else bad "--mode central failed (see $work/central/central.log)"; tail -40 central.log; exit 1; fi
cen=$work/central/results_central/central

# ---- 5. planted signals and output contracts --------------------------------------------------------
if python3 "$repo/tests/scripts/check_central.py" check "$cen" "$fx2/TRUTH.tsv" > check.log; then
  ok "central: planted signals, scopes, TWAS inputs and summary"
else
  bad "central checks"
fi
cat check.log

# ---- 6. independent IVW recalculation (R) -----------------------------------------------------------
cat gwas.log
python3 "$repo/tests/scripts/check_central.py" extract "$cen" ivw_check --truth "$fx2/TRUTH.tsv" > extract.log
if "${RSCRIPT:-Rscript}" "$repo/tests/scripts/check_central_ivw.R" ivw_check "$cen" > ivw.log 2>&1; then
  ok "central: IVW recomputed in R; FUSION weight files valid"
else
  bad "central: IVW recalculation or FUSION weight files"
fi
cat extract.log ivw.log

echo "test_central.sh: $pass passed, $fail failed (work dir $work)"
[ "$fail" -eq 0 ]
