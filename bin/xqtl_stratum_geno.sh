#!/usr/bin/env bash
# xqtl_stratum_geno.sh -- genotype QC, genotype PCs and PC-outlier removal within one ancestry stratum.
#
#   xqtl_stratum_geno.sh --pfile cohort --info cohort.info.tsv --samples stratum_EUR.samples.tsv \
#       --anc EUR [--maf 0.01 --mac 5 --hwe 1e-6 --geno 0.05 --mind 0.05 --npcs 5 --sd 6 --threads N]
#
# Atlas logic: variant filters within the stratum; cohort-only PCA (--indep-pairwise 200 50 0.2,
# --pca 20; 22_genotype_pca.sh); PC1-PC5 are the covariates; samples beyond 6 SD on any of them
# are removed and the MAC floor is re-applied to the remaining samples (24_outlier_exclusion.py).
# PCs are not recomputed after outlier removal (as in the atlas).
#
# Outputs (prefix ANC):
#   ANC.{pgen,pvar,psam}        final genotypes (dosages) for mapping
#   ANC.geno_pcs.tsv            genotype_id + PC1..PCn (sample-level; never exported)
#   ANC.geno_outliers.tsv       removed samples (sample-level; never exported)
#   ANC.samples.tsv             rnaseq_id, genotype_id after outlier removal
#   ANC.variant_qc.tsv          per-variant summary (export table D)
#   ANC.geno_qc_counts.tsv      step, n_variants, n_samples
set -euo pipefail
pfile="" info="" samples="" anc="" maf=0.01 mac=5 hwe=1e-6 geno=0.05 mind=0.05 npcs=5 sd=6 threads=1
while [ $# -gt 0 ]; do
  case "$1" in
    --pfile) pfile=$2; shift 2 ;;
    --info) info=$2; shift 2 ;;
    --samples) samples=$2; shift 2 ;;
    --anc) anc=$2; shift 2 ;;
    --maf) maf=$2; shift 2 ;;
    --mac) mac=$2; shift 2 ;;
    --hwe) hwe=$2; shift 2 ;;
    --geno) geno=$2; shift 2 ;;
    --mind) mind=$2; shift 2 ;;
    --npcs) npcs=$2; shift 2 ;;
    --sd) sd=$2; shift 2 ;;
    --threads) threads=$2; shift 2 ;;
    *) echo "xqtl_stratum_geno.sh: unknown option $1" >&2; exit 2 ;;
  esac
done
for v in pfile info samples anc; do [ -n "${!v}" ] || { echo "xqtl_stratum_geno.sh: --$v is required" >&2; exit 2; }; done

counts="$anc.geno_qc_counts.tsv"
printf 'step\tn_variants\tn_samples\n' > "$counts"
nv() { grep -vc '^#' "$1.pvar"; }
ns() { grep -vc '^#' "$1.psam"; }

awk -F'\t' 'NR>1{print $2}' "$samples" | sort -u > keep.txt
plink2 --pfile "$pfile" --keep keep.txt --output-chr chrM --make-pgen --threads "$threads" --out s0 > /dev/null
printf 'stratum_input\t%s\t%s\n' "$(nv s0)" "$(ns s0)" >> "$counts"

plink2 --pfile s0 --mind "$mind" --geno "$geno" --maf "$maf" --mac "$mac" --hwe "$hwe" \
  --output-chr chrM --make-pgen --threads "$threads" --out s1 > /dev/null
printf 'variant_filters\t%s\t%s\n' "$(nv s1)" "$(ns s1)" >> "$counts"

n1=$(ns s1)
npca=20; [ "$n1" -le 21 ] && npca=$((n1 - 1))
plink2 --pfile s1 --indep-pairwise 200 50 0.2 --threads "$threads" --out pca_prune > /dev/null
plink2 --pfile s1 --extract pca_prune.prune.in --pca "$npca" --threads "$threads" --out pca > /dev/null

# PC outliers (|PC - mean| > sd * SD on PC1..npcs)
awk -v k="$npcs" -v t="$sd" -F'\t' '
  NR==1 { for (i = 1; i <= NF; i++) if ($i ~ /^#?IID$/) id = i; first = id + 1; next }
  { n++; ids[n] = $id; for (j = 0; j < k; j++) { v[n, j] = $(first + j); s[j] += $(first + j); ss[j] += $(first + j)^2 } }
  END {
    for (j = 0; j < k; j++) { m[j] = s[j] / n; sdv[j] = sqrt((ss[j] - n * m[j]^2) / (n - 1)) }
    print "genotype_id\treason"
    for (i = 1; i <= n; i++) for (j = 0; j < k; j++)
      if (sdv[j] > 0 && (v[i, j] - m[j] > t * sdv[j] || m[j] - v[i, j] > t * sdv[j])) { printf "%s\tPC%d beyond %s SD\n", ids[i], j + 1, t; break }
  }' pca.eigenvec > "$anc.geno_outliers.tsv"
n_out=$(($(wc -l < "$anc.geno_outliers.tsv") - 1))
echo "$anc: $n_out genotype PC outlier(s)"

awk -v k="$npcs" -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i ~ /^#?IID$/) id=i; printf "genotype_id"; for(j=1;j<=k;j++) printf "\tPC%d", j; print ""; next}
  {printf "%s", $id; for(j=1;j<=k;j++) printf "\t%s", $(id+j); print ""}' pca.eigenvec > pcs_all.tsv

if [ "$n_out" -gt 0 ]; then
  tail -n +2 "$anc.geno_outliers.tsv" | cut -f1 > remove.txt
  plink2 --pfile s1 --remove remove.txt --mac "$mac" --maf "$maf" --output-chr chrM --make-pgen pvar-cols=-info --threads "$threads" --out "$anc" > /dev/null
  awk -F'\t' 'NR==FNR{r[$1]=1; next} !($1 in r)' remove.txt pcs_all.tsv > "$anc.geno_pcs.tsv"
else
  plink2 --pfile s1 --output-chr chrM --make-pgen pvar-cols=-info --threads "$threads" --out "$anc" > /dev/null
  mv pcs_all.tsv "$anc.geno_pcs.tsv"
fi
printf 'pc_outliers_removed\t%s\t%s\n' "$(nv "$anc")" "$(ns "$anc")" >> "$counts"

grep -v '^#' "$anc.psam" | cut -f1 | sort > final_ids.txt
awk -F'\t' 'NR==FNR{k[$1]=1; next} FNR==1 || ($2 in k)' final_ids.txt "$samples" > "$anc.samples.tsv"

# variant summary (table D): dosage allele frequency, HWE p (hard calls), imputation R2, typed flag
plink2 --pfile "$anc" --freq cols=+pos --hardy cols=+pos --output-chr chrM --threads "$threads" --out vq > /dev/null
python3 - "$anc" "$info" <<'PY'
import sys
import pandas as pd
anc, info = sys.argv[1], sys.argv[2]
f = pd.read_csv("vq.afreq", sep="\t").rename(columns={"#CHROM": "chrom", "POS": "pos", "ID": "variant_id",
                                                     "REF": "ref", "ALT": "alt", "ALT_FREQS": "af", "OBS_CT": "obs_ct"})
h = pd.read_csv("vq.hardy", sep="\t")[["ID", "P"]].rename(columns={"ID": "variant_id", "P": "hwe_p"})
i = pd.read_csv(info, sep="\t", dtype={"variant_id": str})
d = f[["variant_id", "chrom", "pos", "ref", "alt", "af", "obs_ct"]].merge(h, on="variant_id", how="left")
d = d.merge(i, on="variant_id", how="left")
d["maf"] = d["af"].where(d["af"] <= 0.5, 1 - d["af"])
d["n"] = (d.pop("obs_ct") // 2).astype(int)
d[["variant_id", "chrom", "pos", "ref", "alt", "af", "maf", "n", "hwe_p", "rsq", "typed"]].to_csv(
    f"{anc}.variant_qc.tsv", sep="\t", index=False, float_format="%.6g")
PY
rm -f s0.* s1.* keep.txt remove.txt final_ids.txt vq.* pca_prune.*
cat "$counts"
