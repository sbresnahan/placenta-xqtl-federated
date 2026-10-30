#!/usr/bin/env bash
# xqtl_geno_qc.sh -- post-imputation QC for one chromosome (TOPMed *.dose.vcf.gz).
#
#   xqtl_geno_qc.sh --vcf chr22.dose.vcf.gz --samples genotype_ids.txt --chrom chr22 \
#                   --rsq-min 0.4 --out chr22 [--threads N]
#
# Keeps the cohort's genotyped RNA-seq samples (genotype_id list), biallelic variants with
# imputation INFO/R2 >= --rsq-min, sets IDs to chr:pos:ref:alt, and converts dosages (FORMAT/DS)
# to a PLINK 2 fileset. MAF/MAC/HWE filters are applied later, per ancestry stratum.
# Outputs: OUT.pgen/.pvar/.psam, OUT.info.tsv (variant_id, rsq, typed), OUT.counts.tsv
set -euo pipefail
vcf="" samples="" chrom="" rsq=0.4 out="" threads=1
while [ $# -gt 0 ]; do
  case "$1" in
    --vcf) vcf=$2; shift 2 ;;
    --samples) samples=$2; shift 2 ;;
    --chrom) chrom=$2; shift 2 ;;
    --rsq-min) rsq=$2; shift 2 ;;
    --out) out=$2; shift 2 ;;
    --threads) threads=$2; shift 2 ;;
    *) echo "xqtl_geno_qc.sh: unknown option $1" >&2; exit 2 ;;
  esac
done
[ -n "$vcf" ] && [ -n "$samples" ] && [ -n "$out" ] || { echo "usage: see header" >&2; exit 2; }

bcftools query -l "$vcf" | sort -u > vcf_samples.txt
sort -u "$samples" > want.txt
comm -12 vcf_samples.txt want.txt > keep.txt
n_keep=$(wc -l < keep.txt)
if [ "$n_keep" -eq 0 ]; then
  echo "xqtl_geno_qc.sh: none of the metadata genotype_id values are samples in $vcf" >&2
  echo "  first VCF samples: $(head -3 vcf_samples.txt | tr '\n' ' ')" >&2
  echo "  first genotype_id values: $(head -3 want.txt | tr '\n' ' ')" >&2
  exit 1
fi
if ! bcftools view -h "$vcf" | grep -q '^##INFO=<ID=R2,'; then
  echo "xqtl_geno_qc.sh: $vcf has no INFO/R2 field (expected TOPMed/Michigan imputation output)" >&2
  exit 1
fi
if ! bcftools view -h "$vcf" | grep -q '^##FORMAT=<ID=DS,'; then
  echo "xqtl_geno_qc.sh: $vcf has no FORMAT/DS (dosage) field" >&2
  exit 1
fi

n_in=$(bcftools view -H "$vcf" | wc -l)
bcftools view --threads "$threads" -S keep.txt -m2 -M2 -i "INFO/R2>=${rsq}" "$vcf" -Ou \
  | bcftools annotate --threads "$threads" --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o qc.vcf.gz
bcftools query -f '%ID\t%INFO/R2\t%INFO/TYPED\n' qc.vcf.gz \
  | awk 'BEGIN{OFS="\t"; print "variant_id","rsq","typed"} {print $1, $2, ($3=="1"?1:0)}' > "$out.info.tsv"
n_out=$(($(wc -l < "$out.info.tsv") - 1))

if [ "$n_out" -gt 0 ]; then
  # Dosages are kept for mapping; hard calls are the nearest genotype (never missing because a
  # dosage is non-integer), so --geno/--mind later only see genuinely missing genotypes.
  plink2 --vcf qc.vcf.gz dosage=DS --hard-call-threshold 0.4999999 --output-chr chrM \
    --rm-dup exclude-all --make-pgen \
    --threads "$threads" --out "$out" > /dev/null
else
  echo "xqtl_geno_qc.sh: no variants left on $chrom after QC" >&2
fi
rm -f qc.vcf.gz
printf 'chrom\tvariants_in\tvariants_rsq_biallelic\tsamples_vcf\tsamples_kept\n%s\t%s\t%s\t%s\t%s\n' \
  "$chrom" "$n_in" "$n_out" "$(wc -l < vcf_samples.txt)" "$n_keep" > "$out.counts.tsv"
cp vcf_samples.txt "$out.vcf_samples.txt"
