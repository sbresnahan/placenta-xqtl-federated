#!/usr/bin/env bash
# xqtl_geno_merge.sh -- merge per-chromosome PLINK 2 filesets from xqtl_geno_qc.sh.
#
#   xqtl_geno_merge.sh --out cohort [--threads N] chr1 chr2 ... chr22      (filesets by prefix)
#
# Outputs: OUT.pgen/.pvar/.psam, OUT.info.tsv (concatenated), OUT.qc_counts.tsv
set -euo pipefail
out="" threads=1 prefixes=()
while [ $# -gt 0 ]; do
  case "$1" in
    --out) out=$2; shift 2 ;;
    --threads) threads=$2; shift 2 ;;
    *) prefixes+=("$1"); shift ;;
  esac
done
[ -n "$out" ] && [ "${#prefixes[@]}" -gt 0 ] || { echo "usage: xqtl_geno_merge.sh --out OUT PREFIX..." >&2; exit 2; }

have=()
for p in "${prefixes[@]}"; do [ -s "$p.pgen" ] && have+=("$p"); done
[ "${#have[@]}" -gt 0 ] || { echo "xqtl_geno_merge.sh: no chromosome has variants left after QC" >&2; exit 1; }

# all chromosomes must hold the same samples in the same order
ref=$(cut -f1 "${have[0]}.psam" | md5sum)
for p in "${have[@]}"; do
  [ "$(cut -f1 "$p.psam" | md5sum)" = "$ref" ] || { echo "xqtl_geno_merge.sh: $p has a different sample set from ${have[0]}" >&2; exit 1; }
done

if [ "${#have[@]}" -eq 1 ]; then
  plink2 --pfile "${have[0]}" --output-chr chrM --make-pgen --threads "$threads" --out "$out" > /dev/null
else
  printf '%s\n' "${have[@]}" > merge_list.txt
  plink2 --pmerge-list merge_list.txt pfile --output-chr chrM --sort-vars --make-pgen \
    --threads "$threads" --out "$out" > /dev/null
fi
{ echo -e "variant_id\trsq\ttyped"; for p in "${prefixes[@]}"; do [ -f "$p.info.tsv" ] && tail -n +2 "$p.info.tsv"; done; } > "$out.info.tsv"
{ head -1 "${prefixes[0]}.counts.tsv"; for p in "${prefixes[@]}"; do [ -f "$p.counts.tsv" ] && tail -n +2 "$p.counts.tsv"; done; } > "$out.qc_counts.tsv"
# sample IDs present in the imputed VCFs (input to the export leak scan)
cat "${prefixes[@]/%/.vcf_samples.txt}" 2>/dev/null | sort -u > "$out.vcf_samples.txt"
echo "merged ${#have[@]} chromosome(s): $(grep -vc '^#' "$out.pvar") variants, $(grep -vc '^#' "$out.psam") samples"
