#!/usr/bin/env bash
# MDACC only: build the restricted file xqtl-restricted-<version>.tar.gz from atlas outputs.
#
#   scripts/make_restricted_bundle.sh --env scripts/packaging/seadragon.env --out DIR [--dry-run]
#
# Contents (see bin/xqtl_restricted.py): the frozen leafCutter cluster set (atlas pooled
# clusters, merged across ancestries) and the MDACC signal list (atlas lead + conditionally
# independent signals, 9 modalities, all ancestries, canonical IDs). Junction coordinates and
# variant/feature IDs only: no per-sample data.
#
# The file is NOT published. Send it by institutional secure email (or as an encrypted
# archive, password by a separate channel) to approved sites only.
# Needs python3 with pandas + pyarrow (the xqtl environment or the image).
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
envf="" out="" dry=0
while [ $# -gt 0 ]; do
  case "$1" in
    --env) envf=$2; shift 2 ;;
    --out) out=$2; shift 2 ;;
    --dry-run) dry=1; shift ;;
    -h|--help) sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown argument $1" >&2; exit 2 ;;
  esac
done
[ -n "$envf" ] && [ -n "$out" ] || { sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//' >&2; exit 2; }
# shellcheck disable=SC1090
source "$envf"
: "${RESTRICTED_VERSION:?}" "${ANCESTRIES:?}" "${ATLAS_SPLICING_COUNTS:?}" "${ATLAS_QTL_RESULTS:?}"
: "${QVALUE_MAX:=0.05}" "${INCLUDE_UNGROUPED:=no}"
log() { printf '[restricted] %s\n' "$*"; }

log "version $RESTRICTED_VERSION; ancestries: $ANCESTRIES; q <= $QVALUE_MAX; ungrouped scans: $INCLUDE_UNGROUPED"
missing=0
for anc in $ANCESTRIES; do
  p=${ATLAS_SPLICING_COUNTS//\{ANC\}/$anc}
  if [ -s "$p" ]; then log "splicing $anc: $(du -h "$p" | cut -f1) $p"; else log "splicing $anc: MISSING $p"; missing=$((missing + 1)); fi
  for mod in expression isoform_expression isoforms splicing intron_retention alt_TSS alt_polyA RNA_editing stability; do
    found=""
    for suf in cisqtl cisqtl_independent; do
      [ -s "$ATLAS_QTL_RESULTS/${anc}_${mod}_${suf}.parquet" ] && found="$found $suf"
    done
    log "  $anc $mod:${found:- MISSING}"
    [ -z "$found" ] && missing=$((missing + 1))
  done
done
out_file="$out/xqtl-restricted-$RESTRICTED_VERSION.tar.gz"
if [ "$dry" -eq 1 ]; then
  log "dry run: would write $out_file (+ .sha256)"
  [ "$missing" -eq 0 ] && log "dry run: READY" || log "dry run: NOT READY ($missing missing input(s))"
  exit 0
fi
[ "$missing" -eq 0 ] || { log "ERROR: $missing input(s) missing"; exit 1; }

mkdir -p "$out"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
ung=(); [ "$INCLUDE_UNGROUPED" = yes ] && ung=(--include-ungrouped)
python3 "$here/scripts/lib/atlas_restricted_tables.py" \
  --ancestries "$ANCESTRIES" --splicing-counts "$ATLAS_SPLICING_COUNTS" \
  --qtl-results "$ATLAS_QTL_RESULTS" --qvalue-max "$QVALUE_MAX" "${ung[@]}" \
  --out-clusters "$tmp/frozen_clusters.tsv.gz" --out-signals "$tmp/mdacc_signals.tsv.gz"
python3 "$here/bin/xqtl_restricted.py" make --clusters "$tmp/frozen_clusters.tsv.gz" \
  --signals "$tmp/mdacc_signals.tsv.gz" --version "$RESTRICTED_VERSION" --out "$out_file"
python3 "$here/bin/xqtl_restricted.py" verify "$out_file" --expected-version "$RESTRICTED_VERSION"
(cd "$out" && sha256sum "$(basename "$out_file")" > "$(basename "$out_file").sha256")
log "wrote $out_file"
log "send it only to approved sites, by secure email; keep a copy with the release records"
