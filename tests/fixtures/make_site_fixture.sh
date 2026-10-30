#!/usr/bin/env bash
# make_site_fixture.sh -- regenerate the site-run test fixture (tests/fixtures/site and the
# ancestry panel in tests/fixtures/refs_mini). Developer tool; the outputs are committed.
#
#   tests/fixtures/make_site_fixture.sh /path/1kGP_chr22.vcf.gz /path/samples_ped_population.txt
#
# Genotypes are real public data: 1000 Genomes high-coverage phased panel (IGSR,
# 20220422_3202_phased_SNV_INDEL_SV, chr22), free to use. The cohort is 130 unrelated
# 1kGP samples (60 EUR, 60 EAS, 10 AFR) written as a TOPMed-style chr22.dose.vcf.gz with
# synthetic R2 and dosages; the PCA panel is 160 other 1kGP samples (40 each of AFR, EAS, EUR, SAS; AMR
# is omitted because single-chromosome PCs cannot separate it from EUR).
# Everything else (phenotypes, QC metrics, cell proportions, sex, GA, ppBMI, gdm) is a
# SYNTHETIC SOFTWARE FIXTURE with planted cis effects (tests/fixtures/site/TRUTH.tsv).
# Needs: bcftools, bgzip, tabix, plink2, python3 with numpy + pandas.
set -euo pipefail
vcf=${1:?usage: make_site_fixture.sh KGP_CHR22_VCF PED}
ped=${2:?usage: make_site_fixture.sh KGP_CHR22_VCF PED}
here=$(cd "$(dirname "$0")" && pwd)
out="$here/site"
panel="$here/refs_mini/ancestry"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$out/genotypes" "$out/phenotypes" "$panel"

python3 "$here/make_site_fixture.py" select --ped "$ped" --tmp "$tmp"

# biallelic SNPs, MAF >= 0.01 over the 330 selected samples
bcftools view -S "$tmp/all_samples.txt" -m2 -M2 -v snps "$vcf" -Ou \
  | bcftools view -q 0.01:minor -Ou \
  | bcftools annotate -x INFO,^FORMAT/GT --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o "$tmp/sub.vcf.gz"
tabix -p vcf "$tmp/sub.vcf.gz"

# dense mapping region (chr22:19.0-21.5 Mb) + sparse genome-wide-like background (MAF >= 0.05)
plink2 --vcf "$tmp/sub.vcf.gz" --chr chr22 --from-bp 19000000 --to-bp 21500000 --bp-space 600 \
  --write-snplist --out "$tmp/dense" > /dev/null
plink2 --vcf "$tmp/sub.vcf.gz" --maf 0.05 --bp-space 1500 --write-snplist --out "$tmp/sparse" > /dev/null
sort -u "$tmp/dense.snplist" "$tmp/sparse.snplist" > "$tmp/keep.snplist"
bcftools view -i "ID=@$tmp/keep.snplist" "$tmp/sub.vcf.gz" -Oz -o "$tmp/keep.vcf.gz"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' "$tmp/keep.vcf.gz" > "$tmp/gt.tsv"
bcftools query -l "$tmp/keep.vcf.gz" > "$tmp/gt_samples.txt"

python3 "$here/make_site_fixture.py" build --tmp "$tmp" --out "$out"
bgzip -f "$out/genotypes/chr22.dose.vcf"

# ancestry panel, built exactly like scripts/make_public_reference_record.sh (scaled down)
bcftools view -S "$tmp/panel_samples.txt" "$tmp/keep.vcf.gz" -Oz -o "$tmp/panel.vcf.gz"
plink2 --vcf "$tmp/panel.vcf.gz" --make-pgen --out "$tmp/kgp" > /dev/null
plink2 --pfile "$tmp/kgp" --autosome --snps-only just-acgt --max-alleles 2 \
  --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1000 --rm-dup exclude-all \
  --maf 0.05 --indep-pairwise 1000 100 0.2 --out "$tmp/prune" > /dev/null
plink2 --pfile "$tmp/kgp" --autosome --snps-only just-acgt --max-alleles 2 \
  --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1000 --rm-dup exclude-all \
  --extract "$tmp/prune.prune.in" --sort-vars --make-pgen vzs pvar-cols=-xheader,-info \
  --out "$panel/1kGP_hg38_pca_panel" > /dev/null
rm -f "$panel/1kGP_hg38_pca_panel.log"
cp "$tmp/panel_superpop.tsv" "$panel/1kg_superpop.tsv"

# refs_mini manifest (paths relative to refs_mini)
(
  cd "$here/refs_mini"
  printf 'path\tsha256\tbytes\tdescription\tsource\tlicense\n'
  for f in README.txt ancestry/1kGP_hg38_pca_panel.pgen ancestry/1kGP_hg38_pca_panel.pvar.zst \
           ancestry/1kGP_hg38_pca_panel.psam ancestry/1kg_superpop.tsv; do
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$f" "$(sha256sum "$f" | cut -d' ' -f1)" "$(stat -c %s "$f")" \
      "test fixture" "1000 Genomes high-coverage chr22 subset (IGSR)" "IGSR: free, no restrictions"
  done
) > "$tmp/MANIFEST.tsv"
sed -i '2s/test fixture\t1000 Genomes high-coverage chr22 subset (IGSR)\tIGSR: free, no restrictions/smoke-test fixture\tsynthetic\tCC0/' "$tmp/MANIFEST.tsv"
mv "$tmp/MANIFEST.tsv" "$here/refs_mini/MANIFEST.tsv"
echo "site fixture written to $out; panel to $panel"
du -sh "$out" "$panel"
