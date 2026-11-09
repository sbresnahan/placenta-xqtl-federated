#!/usr/bin/env bash
# xqtl_impprep.sh -- --mode imputation_prep steps (atlas 01_genotype_imputation/filter_vcf_for_imputation).
#
#   xqtl_impprep.sh convert KIND RAW BUILD THREADS
#       KIND vcf (VCF/BCF file RAW) | bed (PLINK 1 prefix RAW) | pgen (PLINK 2 prefix RAW).
#       Writes converted.vcf.gz(+.csi): chromosomes renamed to chr1..chr22/chrX, other contigs
#       dropped, sorted; chroms.tsv (chrom, records); sample_missing.smiss (autosomes, plink2);
#       chrX_gt.tsv + chrX_samples.txt (non-PAR chrX SNPs with MAF >= 0.05, for the sex check).
#       Duplicate sample IDs stop the run.
#   xqtl_impprep.sh filter IN.vcf.gz CHROM GENOME GENO_MISS MAF HWE THREADS
#       Biallelic SNPs; INFO tags recomputed (+fill-tags); exclude F_MISSING > GENO_MISS,
#       MAF < MAF, HWE < HWE; REF/ALT swapped to the GRCh38 reference where the ALT allele matches
#       it (plink2 --ref-from-fa: TOPMed r3 no longer fixes allele order); sites whose REF still
#       does not match are dropped (bcftools norm --check-ref x); sorted.
#       Writes CHROM.vcf.gz(+.tbi) and CHROM.counts.tsv (chrom, step, n).
set -euo pipefail
cmd=${1:?usage: xqtl_impprep.sh convert|filter ...}
shift

count() { bcftools index -n "$1"; }

if [ "$cmd" = convert ]; then
  kind=$1 raw=$2 build=$3 threads=${4:-1}
  # ---- duplicate sample IDs (checked on the raw header: htslib refuses duplicates later) ----
  case $kind in
    vcf) if [[ $raw == *.bcf ]]; then ids=$(bcftools query -l "$raw"); else ids=$( { gzip -cdf "$raw" || true; } | awk -F'\t' '/^#CHROM/{for (i = 10; i <= NF; i++) print $i; exit}'); fi ;;
    bed) ids=$(awk '{print $2}' "$raw.fam") ;;
    pgen) ids=$(awk 'NR==1 && /^#IID/ {c=1; next} NR==1 && /^#/ {c=2; next} NR==1 {c=2} {print $c}' "$raw.psam") ;;
    *) echo "unknown input kind $kind" >&2; exit 2 ;;
  esac
  dups=$(printf '%s\n' "$ids" | sort | uniq -d | awk 'NR <= 5' | tr '\n' ' ')
  if [ -n "$dups" ]; then
    echo "ERROR: duplicate sample IDs in the genotype input (e.g. $dups). Make IDs unique before imputation." >&2
    exit 3
  fi
  printf '%s\n' "$ids" | grep -c . > n_samples.txt || true
  # ---- to VCF ----
  if [ "$kind" = vcf ]; then
    src=$raw
  else
    flag=--bfile; [ "$kind" = pgen ] && flag=--pfile
    plink2 $flag "$raw" --chr 1-22,X,XY --output-chr chrM --export vcf bgz id-paste=iid \
      --threads "$threads" --out plink_export > /dev/null
    src=plink_export.vcf.gz
  fi
  : > rename.txt
  for c in $(seq 1 22) X Y MT; do printf '%s\tchr%s\n' "$c" "$c" >> rename.txt; done
  printf '23\tchrX\n25\tchrX\nXY\tchrX\nchrXY\tchrX\nchr23\tchrX\n' >> rename.txt
  keep=$(printf 'chr%s,' $(seq 1 22) X); keep=${keep%,}
  bcftools annotate --rename-chrs rename.txt "$src" -Ou \
    | bcftools view -t "$keep" -Ou \
    | bcftools sort -m 2G -T sort_tmp -Oz -o converted.vcf.gz
  bcftools index -f converted.vcf.gz
  bcftools index -s converted.vcf.gz | awk 'BEGIN{OFS="\t"; print "chrom", "records"} {print $1, $3}' > chroms.tsv
  if ! grep -qP '^chr([1-9]|1[0-9]|2[0-2])\t' chroms.tsv; then
    echo "ERROR: no autosomal variants after renaming chromosomes (contigs found: $(bcftools index -s "$src" 2>/dev/null | awk 'NR <= 5 {print $1}' | tr '\n' ' '))" >&2
    exit 4
  fi
  # ---- per-sample missingness on the autosomes ----
  plink2 --vcf converted.vcf.gz --chr 1-22 --vcf-half-call m --allow-extra-chr --missing sample-only \
    --threads "$threads" --out sample_missing > /dev/null
  # ---- chrX genotypes for the sex check (non-PAR, common SNPs) ----
  bcftools query -l converted.vcf.gz > chrX_samples.txt
  if grep -qP '^chrX\t' chroms.tsv; then
    if [ "$build" = GRCh37 ]; then par="chrX:60001-2699520,chrX:154931044-155260560"
    else par="chrX:10001-2781479,chrX:155701383-156030895"; fi
    bcftools view -r chrX -t "^$par" -v snps -m2 -M2 converted.vcf.gz -Ou \
      | bcftools +fill-tags -Ou -- -t MAF \
      | bcftools view -i 'INFO/MAF>=0.05' -Ou \
      | bcftools query -f '[%GT\t]\n' > chrX_gt.tsv
  else
    : > chrX_gt.tsv
  fi

elif [ "$cmd" = filter ]; then
  in=$1 chrom=$2 genome=$3 gmiss=$4 maf=$5 hwe=$6 threads=${7:-1}
  n0=$(bcftools view -H -r "$chrom" "$in" | wc -l)
  bcftools view -r "$chrom" -v snps -m2 -M2 "$in" -Ou \
    | bcftools +fill-tags -Ob -o s1.bcf -- -t AF,MAF,F_MISSING,HWE
  bcftools index -f s1.bcf
  bcftools view -e "INFO/F_MISSING > $gmiss || INFO/MAF < $maf" s1.bcf -Ob -o s2.bcf && bcftools index -f s2.bcf
  bcftools view -e "INFO/HWE < $hwe" s2.bcf -Ob -o s3.bcf && bcftools index -f s3.bcf
  n1=$(count s1.bcf) n2=$(count s2.bcf) n3=$(count s3.bcf)
  swapped=0
  if [ "$n3" -gt 0 ]; then
    plink2 --bcf s3.bcf --ref-from-fa force --fa "$genome" --output-chr chrM --vcf-half-call m \
      --allow-extra-chr --export vcf bgz --threads "$threads" --out s4 > /dev/null
    swapped=$(grep -oP 'ref-from-fa force: \K[0-9]+(?= variants? changed)' s4.log || echo 0)
    bcftools norm --check-ref x -f "$genome" s4.vcf.gz -Ou 2> norm.log \
      | bcftools sort -m 1G -T sort_tmp -Oz -o "$chrom.vcf.gz"
  else
    bcftools view s3.bcf -Oz -o "$chrom.vcf.gz"
  fi
  tabix -f -p vcf "$chrom.vcf.gz"
  n4=$(bcftools view -H "$chrom.vcf.gz" | wc -l)
  {
    printf 'chrom\tstep\tn\n'
    printf '%s\tinput\t%s\n' "$chrom" "$n0"
    printf '%s\tbiallelic_snps\t%s\n' "$chrom" "$n1"
    printf '%s\tpass_missing_maf\t%s\n' "$chrom" "$n2"
    printf '%s\tpass_hwe\t%s\n' "$chrom" "$n3"
    printf '%s\tref_alt_swapped\t%s\n' "$chrom" "$swapped"
    printf '%s\tref_match\t%s\n' "$chrom" "$n4"
  } > "$chrom.counts.tsv"
  rm -f s1.bcf* s2.bcf* s3.bcf* s4.vcf.gz
else
  echo "unknown command $cmd" >&2; exit 2
fi
