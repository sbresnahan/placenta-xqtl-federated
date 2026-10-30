#!/usr/bin/env bash
# M3 acceptance test: --mode site on the chr22 1000 Genomes fixture, end to end.
#   1. -profile test,<engine> --mode site runs genotype QC -> ancestry -> mapping -> LD -> export;
#   2. the planted cis effects in tests/fixtures/site/TRUTH.tsv are recovered in EUR and EAS
#      (q <= 0.05 in B; correct sign and p <= 1e-3 at the planted variant in A);
#   3. AFR (n = 10 < ancestry.min_n) is reported in I_strata but not mapped or exported;
#   4. a sample ID planted in a site output makes --mode export fail and leaves no export behind;
#   5. after removing it, --mode export rebuilds the same export.
#
#   tests/scripts/test_site.sh [WORKDIR] [ENGINE]     ENGINE: dev (default; tools on PATH) | docker
# Extra Nextflow options can be passed in NXF_OPTS_EXTRA (e.g. "--image ghcr.io/...:dev").
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
work=${1:-$(mktemp -d)}
engine=${2:-dev}
mkdir -p "$work"; work=$(cd "$work" && pwd)
cd "$work"
pass=0 fail=0
ok()  { echo "PASS  $*"; pass=$((pass + 1)); }
bad() { echo "FAIL  $*"; fail=$((fail + 1)); }
out=results_test
exp=$out/export/TESTCOHORT
# shellcheck disable=SC2086
nxf() { nextflow run "$repo/main.nf" -profile "test,$engine" ${NXF_OPTS_EXTRA:-} "$@"; }

# ---- 1. site run -------------------------------------------------------------------------
if nxf --mode site -resume > site.log 2>&1; then ok "--mode site completed"; else bad "--mode site failed (see $work/site.log)"; tail -40 site.log; exit 1; fi

# ---- 2. planted signals + manifest ---------------------------------------------------------
if python3 "$repo/tests/scripts/check_site_truth.py" "$exp" "$repo/tests/fixtures/site/TRUTH.tsv" EUR EAS > truth.log; then
  ok "planted effects recovered in EUR and EAS; manifest complete"
else
  bad "truth checks"
fi
cat truth.log

# ---- 3. AFR below min_n -----------------------------------------------------------------------
if grep -qP '^TESTCOHORT\tAFR\t10\t10\tbelow_min_n$' "$exp/TESTCOHORT_cohort_I_strata.tsv" && [ ! -e "$exp/AFR" ] && [ ! -e "$out/site/AFR" ]; then
  ok "AFR (n=10) listed as below_min_n and not mapped"
else
  bad "AFR stratum handling"
fi
if grep -q 'no_genotypes' "$exp/TESTCOHORT_cohort_I_strata.tsv"; then ok "sample without genotypes counted as no_genotypes"; else bad "no_genotypes row"; fi

# ---- 4. leak injection ------------------------------------------------------------------------
f=$out/site/EUR/EUR_expression.features.tsv
cp "$f" features.bak
awk -F'\t' 'BEGIN{OFS="\t"} NR==2{$1="HG00351"} {print}' features.bak > "$f"   # a genotype ID of the cohort
if nxf --mode export > export_bad.log 2>&1; then
  bad "export with a planted sample ID passed the leak scan"
else
  if grep -q 'LEAK SCAN FAILED' export_bad.log && [ ! -e "$out/export" ]; then
    ok "planted sample ID: export failed and no export folder was left"
  else
    bad "export failed but not via the leak scan, or an export folder remained (see export_bad.log)"
  fi
fi
cp features.bak "$f"

# ---- 5. clean rebuild -----------------------------------------------------------------------
if nxf --mode export > export.log 2>&1 && [ -s "$exp/MANIFEST.tsv" ]; then
  ok "--mode export rebuilds the export ($(($(wc -l < "$exp/MANIFEST.tsv") - 1)) files)"
else
  bad "--mode export rebuild (see export.log)"
fi

echo "test_site.sh: $pass passed, $fail failed (work dir $work)"
[ "$fail" -eq 0 ]
