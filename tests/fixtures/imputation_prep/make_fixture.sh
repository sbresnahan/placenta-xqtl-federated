#!/usr/bin/env bash
# Rebuild tests/fixtures/imputation_prep (used by tests/scripts/test_imputation_prep.sh).
#   kgp37_chr22.vcf.gz            1000 Genomes phase 3 (GRCh37, contig "22"), 22:20000000-20600000, the first 40
#                                 site-fixture samples (site-fixture order) present in 1000 Genomes, MAF >= 0.02 (SNPs and indels)
#   hg19ToHg38.chr22.over.chain.gz  the chr22 -> chr22 chains of the UCSC hg19ToHg38 chain
#   kgp38_chr22.{bed,bim,fam}     3,000 biallelic SNPs of the site fixture (GRCh38) written with REF = major
#                                 allele, so ~600 sites need a REF/ALT swap back to GRCh38
#
#   tests/fixtures/imputation_prep/make_fixture.sh CHAIN    CHAIN: hg19ToHg38.over.chain.gz (refs_dir/liftover)
# Needs bcftools, tabix, plink2 and network access to the 1000 Genomes FTP (tabix remote fetch, ~1 MB).
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
chain=${1:?usage: make_fixture.sh hg19ToHg38.over.chain.gz}
site=$here/../site/genotypes/chr22.dose.vcf.gz
kgp=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
cd "$tmp"

# ---- GRCh37 VCF ----
bcftools query -l "$site" > site.ids
bcftools view -h "$kgp" | awk -F'\t' '/^#CHROM/{for (i = 10; i <= NF; i++) print $i}' > kgp.ids
awk 'NR == FNR {k[$1] = 1; next} ($1 in k)' kgp.ids site.ids | awk 'NR <= 40' > s40.txt   # site order
bcftools view -r 22:20000000-20600000 -S s40.txt "$kgp" -Ou \
  | bcftools view -i 'MAF>=0.02' -Oz -o "$here/kgp37_chr22.vcf.gz"
rm -f ./*.tbi ALL.chr22.*.tbi

# ---- chr22 chain ----
gzip -cdf "$chain" | awk '/^chain/ {p = ($3 == "chr22" && $8 == "chr22")} p' | gzip -n -9 \
  > "$here/hg19ToHg38.chr22.over.chain.gz"

# ---- GRCh38 PLINK with major-allele REF ----
plink2 --vcf "$site" --snps-only just-acgt --max-alleles 2 --thin-count 3000 --seed 7 \
  --maj-ref force --make-bed --out kgp38p > /dev/null
for e in bed bim fam; do cp "kgp38p.$e" "$here/kgp38_chr22.$e"; done
echo "fixture rebuilt in $here"
