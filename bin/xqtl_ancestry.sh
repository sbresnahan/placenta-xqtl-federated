#!/usr/bin/env bash
# xqtl_ancestry.sh -- genetic ancestry of the cohort against the embedded 1000 Genomes panel.
#
#   xqtl_ancestry.sh --pfile cohort --panel REFS/ancestry/1kGP_hg38_pca_panel \
#                    --superpop REFS/ancestry/1kg_superpop.tsv --metadata metadata.tsv \
#                    --min-n 50 --prob-min 0.8 --groups AFR,AMR,EAS,EUR,SAS --out cohort [--threads N]
#
# Atlas logic (01_genotype_imputation/PC_AiR.R): cohort and panel are intersected on
# chr:pos:ref:alt, merged, filtered (MAF 0.05, missingness 0.05) and LD-pruned
# (--indep-pairwise 1000 100 0.1); KING-robust kinship and PC-AiR (1000 Genomes samples as the
# unrelated set) give the PCs; each sample is assigned to the nearest superpopulation centroid
# on PCs 1-10. Added for the federated pipeline: an assignment confidence (posterior of the
# assigned label from a QDA classifier on the panel PCs); samples below --prob-min are "unassigned".
# Outputs (prefix OUT): _ancestry.tsv, _strata.tsv, _pca.png, _pcair_pcs.tsv,
# stratum_<ANC>.samples.tsv for every ancestry group with n >= --min-n.
set -euo pipefail
pfile="" panel="" superpop="" meta="" min_n=50 prob_min=0.8 groups="AFR,AMR,EAS,EUR,SAS" out="" threads=1
while [ $# -gt 0 ]; do
  case "$1" in
    --pfile) pfile=$2; shift 2 ;;
    --panel) panel=$2; shift 2 ;;
    --superpop) superpop=$2; shift 2 ;;
    --metadata) meta=$2; shift 2 ;;
    --min-n) min_n=$2; shift 2 ;;
    --prob-min) prob_min=$2; shift 2 ;;
    --groups) groups=$2; shift 2 ;;
    --out) out=$2; shift 2 ;;
    --threads) threads=$2; shift 2 ;;
    *) echo "xqtl_ancestry.sh: unknown option $1" >&2; exit 2 ;;
  esac
done
for v in pfile panel superpop meta out; do [ -n "${!v}" ] || { echo "xqtl_ancestry.sh: --$v is required" >&2; exit 2; }; done
pz=""; [ -f "$panel.pvar.zst" ] && pz=" vzs"

# shellcheck disable=SC2016  # plink2 --set-all-var-ids template, not shell variables
ids='@:#:$r:$a'
plink2 --pfile "$pfile" --autosome --snps-only just-acgt --max-alleles 2 --output-chr chrM \
  --set-all-var-ids "$ids" --new-id-max-allele-len 1000 --rm-dup exclude-all \
  --maf 0.05 --geno 0.05 --make-pgen --threads "$threads" --out anc_cohort > /dev/null
# shellcheck disable=SC2086
plink2 --pfile "$panel"$pz --autosome --snps-only just-acgt --max-alleles 2 --output-chr chrM \
  --set-all-var-ids "$ids" --new-id-max-allele-len 1000 --rm-dup exclude-all \
  --make-pgen --threads "$threads" --out anc_panel > /dev/null
comm -12 <(grep -v '^#' anc_cohort.pvar | cut -f3 | sort) <(grep -v '^#' anc_panel.pvar | cut -f3 | sort) > common.txt
n_common=$(wc -l < common.txt)
echo "cohort/panel shared variants: $n_common"
[ "$n_common" -ge 100 ] || { echo "xqtl_ancestry.sh: only $n_common variants shared with the 1000 Genomes panel (need >= 100); are the imputed VCFs on GRCh38?" >&2; exit 1; }
overlap=$(comm -12 <(grep -v '^#' anc_cohort.psam | cut -f1 | sort) <(grep -v '^#' anc_panel.psam | cut -f1 | sort) | wc -l)
[ "$overlap" -eq 0 ] || { echo "xqtl_ancestry.sh: $overlap cohort sample IDs are also 1000 Genomes panel IDs; rename them" >&2; exit 1; }

plink2 --pfile anc_cohort --extract common.txt --make-pgen --threads "$threads" --out anc_c2 > /dev/null
plink2 --pfile anc_panel --extract common.txt --make-pgen --threads "$threads" --out anc_p2 > /dev/null
# plink2 cannot yet merge disjoint sample sets on shared variants; go through bcftools.
for f in anc_c2 anc_p2; do
  plink2 --pfile "$f" --export vcf bgz id-paste=iid --output-chr chrM --threads "$threads" --out "$f" > /dev/null
  bcftools index -f "$f.vcf.gz"
done
bcftools merge -m none -Oz -o anc_merged.vcf.gz anc_c2.vcf.gz anc_p2.vcf.gz
plink2 --vcf anc_merged.vcf.gz --output-chr chrM --make-pgen --threads "$threads" --out anc_merged > /dev/null
plink2 --pfile anc_merged --maf 0.05 --geno 0.05 --indep-pairwise 1000 100 0.1 --threads "$threads" --out anc_prune > /dev/null
plink2 --pfile anc_merged --extract anc_prune.prune.in --make-bed --threads "$threads" --out anc_pruned > /dev/null
echo "pruned variants for PC-AiR: $(wc -l < anc_pruned.bim)"

xqtl_ancestry.R --bfile anc_pruned --superpop "$superpop" --metadata "$meta" \
  --cohort-psam anc_cohort.psam --min-n "$min_n" --prob-min "$prob_min" --groups "$groups" \
  --threads "$threads" --out "$out"
rm -f anc_*.pgen anc_*.bed anc_*.vcf.gz anc_*.vcf.gz.csi
