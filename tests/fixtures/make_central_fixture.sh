#!/usr/bin/env bash
# make_central_fixture.sh -- regenerate tests/fixtures/site2, the second cohort of the central
# acceptance test (tests/scripts/test_central.sh). Developer tool; the outputs are committed.
#
#   tests/fixtures/make_central_fixture.sh /path/1kGP_chr22.vcf.gz /path/samples_ped_population.txt
#
# Same inputs as make_site_fixture.sh. The cohort is 120 other unrelated 1kGP samples (60 EUR,
# 60 EAS; none in tests/fixtures/site or in the PCA panel) on the chr22 variants of the first
# cohort, written as a TOPMed-style chr22.dose.vcf.gz with synthetic R2 and dosages. Phenotypes are
# SYNTHETIC with the planted effects of tests/fixtures/site/TRUTH.tsv (other random seeds).
# Needs: bcftools, bgzip, tabix, python3 with numpy + pandas.
set -euo pipefail
vcf=${1:?usage: make_central_fixture.sh KGP_CHR22_VCF PED}
ped=${2:?usage: make_central_fixture.sh KGP_CHR22_VCF PED}
here=$(cd "$(dirname "$0")" && pwd)
out="$here/site2"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$out/genotypes" "$out/phenotypes"

python3 "$here/make_site_fixture.py" select --second --ped "$ped" --tmp "$tmp"

# the first cohort's variants, for the new samples
bcftools query -f '%ID\n' "$here/site/genotypes/chr22.dose.vcf.gz" > "$tmp/keep.snplist"
bcftools view -S "$tmp/cohort_samples.txt" -m2 -M2 -v snps "$vcf" -Ou \
  | bcftools annotate -x INFO,^FORMAT/GT --set-id '%CHROM:%POS:%REF:%ALT' -Ou \
  | bcftools view -i "ID=@$tmp/keep.snplist" -Oz -o "$tmp/keep.vcf.gz"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' "$tmp/keep.vcf.gz" > "$tmp/gt.tsv"
bcftools query -l "$tmp/keep.vcf.gz" > "$tmp/gt_samples.txt"

python3 "$here/make_site_fixture.py" build --second --truth "$here/site/TRUTH.tsv" --tmp "$tmp" --out "$out"
bgzip -f "$out/genotypes/chr22.dose.vcf"
echo "second cohort written to $out"
du -sh "$out"
