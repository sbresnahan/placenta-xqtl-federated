#!/usr/bin/env bash
# End-to-end test of --mode imputation_prep (tests/fixtures/imputation_prep):
#   1. GRCh37 VCF (1000 Genomes phase 3, chr22 slice, 40 samples) -> liftover -> chr22.vcf.gz whose
#      variants sit at the GRCh38 positions and alleles of the same samples' GRCh38 genotypes
#      (tests/fixtures/site/genotypes), plus the report, liftover summary and TOPMed checklist;
#   2. GRCh38 PLINK 1 bed (3,000 SNPs, REF = major allele, ~600 not the GRCh38 REF) -> REF/ALT swapped
#      to GRCh38 with genotypes intact (ALT counts equal those of the source VCF);
#   3. the GRCh37 VCF declared as GRCh38 stops with the "Wrong build?" message and publishes no VCF;
#   4. duplicate sample IDs stop the run.
#
#   tests/scripts/test_imputation_prep.sh [WORKDIR] [ENGINE]     ENGINE: dev (default) | docker
# The chr22 GRCh38 FASTA is cut from XQTL_GENOME (default: the image's FASTA) unless XQTL_CHR22_FA
# points to a ready chr22 FASTA. Extra Nextflow options can be passed in NXF_OPTS_EXTRA.
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
work=${1:-$(mktemp -d)}
engine=${2:-dev}
mkdir -p "$work"; work=$(cd "$work" && pwd)
cd "$work"
fix=$repo/tests/fixtures/imputation_prep
pass=0 fail=0
ok()  { echo "PASS  $*"; pass=$((pass + 1)); }
bad() { echo "FAIL  $*"; fail=$((fail + 1)); }
# shellcheck disable=SC2086
nxf() { nextflow run "$repo/main.nf" -profile "test,$engine" --mode imputation_prep ${NXF_OPTS_EXTRA:-} "$@"; }

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
common=(--genome_fasta "$fa" --liftover_chain "$fix/hg19ToHg38.chr22.over.chain.gz")

params() {  # params() OUT RAW BUILD
  printf 'genotypes:\n  raw: %s\n  raw_build: %s\n  source: array\n' "$2" "$3" > "$1"
}

# ---- 1. GRCh37 VCF ------------------------------------------------------------------------------
params p37.yml "$fix/kgp37_chr22.vcf.gz" GRCh37
if nxf -params-file p37.yml "${common[@]}" --outdir out37 > run37.log 2>&1; then
  ok "GRCh37 VCF: imputation_prep completed"
else
  bad "GRCh37 VCF run failed (see $work/run37.log)"; tail -30 run37.log
fi
if python3 "$repo/tests/scripts/check_imputation_prep.py" lifted out37/imputation_prep \
     "$repo/tests/fixtures/site/genotypes/chr22.dose.vcf.gz" > check37.log 2>&1; then
  ok "GRCh37 -> GRCh38: positions and alleles match the GRCh38 genotypes; report and checklist present"
else
  bad "GRCh37 liftover checks"
fi
cat check37.log

# ---- 2. GRCh38 PLINK with REF = major allele ------------------------------------------------------
params p38.yml "$fix/kgp38_chr22" GRCh38
if nxf -params-file p38.yml "${common[@]}" --outdir out38 > run38.log 2>&1; then
  ok "GRCh38 PLINK: imputation_prep completed"
else
  bad "GRCh38 PLINK run failed (see $work/run38.log)"; tail -30 run38.log
fi
if python3 "$repo/tests/scripts/check_imputation_prep.py" swapped out38/imputation_prep \
     "$repo/tests/fixtures/site/genotypes/chr22.dose.vcf.gz" > check38.log 2>&1; then
  ok "PLINK input: REF/ALT aligned to GRCh38 with genotypes intact"
else
  bad "PLINK allele alignment checks"
fi
cat check38.log

# ---- 3. wrong build ---------------------------------------------------------------------------------
params pwrong.yml "$fix/kgp37_chr22.vcf.gz" GRCh38
if nxf -params-file pwrong.yml "${common[@]}" --outdir outwrong > runwrong.log 2>&1; then
  bad "GRCh37 data declared as GRCh38 passed"
elif grep -q 'Wrong build?' runwrong.log && ! ls outwrong/imputation_prep/chr*.vcf.gz > /dev/null 2>&1; then
  ok "GRCh37 data declared as GRCh38: stopped with 'Wrong build?' and no VCF published"
else
  bad "wrong-build run failed for another reason or published VCFs (see runwrong.log)"
fi

# ---- 4. duplicate sample IDs ------------------------------------------------------------------------
python3 - "$fix/kgp37_chr22.vcf.gz" dup.vcf <<'EOF'
import gzip, sys
with gzip.open(sys.argv[1], "rt") as fh, open(sys.argv[2], "w") as out:
    for i, line in enumerate(fh):
        if line.startswith("#CHROM"):
            f = line.rstrip("\n").split("\t")
            f[10] = f[9]
            line = "\t".join(f) + "\n"
        out.write(line)
        if i > 400:
            break
EOF
params pdup.yml "$work/dup.vcf" GRCh37
if nxf -params-file pdup.yml "${common[@]}" --outdir outdup > rundup.log 2>&1; then
  bad "duplicate sample IDs passed"
elif grep -q 'duplicate sample IDs' rundup.log; then
  ok "duplicate sample IDs stop the run"
else
  bad "duplicate-ID run failed for another reason (see rundup.log)"
fi

echo "test_imputation_prep.sh: $pass passed, $fail failed (work dir $work)"
[ "$fail" -eq 0 ]
