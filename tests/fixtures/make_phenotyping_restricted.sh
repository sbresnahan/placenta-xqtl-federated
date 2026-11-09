#!/usr/bin/env bash
# Build tests/fixtures/phenotyping/xqtl-restricted-phenotyping-test.tar.gz, the restricted file used
# by -profile test_phenotyping.
#
# The real restricted file holds the MDACC frozen leafCutter clusters, which cannot be published. For the
# test we cluster the 8 fixture samples' junctions with the pipeline's own leafCutter settings and
# "freeze" every other cluster (renumbered clu_<1000+k>_<strand> so frozen and local IDs differ).
# A phenotyping run then sees half its clusters matched to the frozen set and half flagged novel,
# which exercises both branches of AGG_SPLICING. Signals are copied from the site test file.
#
# Usage: make_phenotyping_restricted.sh JUNC_DIR
#   JUNC_DIR  regtools junction files <sample>.junc for the fixture samples (STAR 2-pass BAMs,
#             regtools junctions extract -a 8 -m 50 -M 500000 -s XS)
# Needs python3 with pandas on PATH. Writes into tests/fixtures/phenotyping/.
set -euo pipefail

JUNC=$(cd "${1:?usage: make_phenotyping_restricted.sh JUNC_DIR}" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
BIN="$HERE/../../bin"
OUT="$HERE/phenotyping"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

ls "$JUNC"/*.junc > "$TMP/juncfiles.txt"
mkdir -p "$TMP/lc"
python3 "$BIN/leafcutter_cluster_regtools_py3.py" -j "$TMP/juncfiles.txt" -r "$TMP/lc" -o leafcutter \
    --maxintronlen 100000 --minclureads 30 --mincluratio 0.001 > "$TMP/lc.log" 2>&1

# refined lines: "chrom:strand start:end:count start:end:count ..."; keep odd-numbered clusters
python3 - "$TMP/lc/leafcutter_refined" "$TMP/frozen_clusters.tsv" <<'EOF'
import sys
rows, k = [], 0
for i, line in enumerate(open(sys.argv[1]), start=1):
    if i % 2 == 0:
        continue
    k += 1
    loc, *introns = line.split()
    chrom, strand = loc.rsplit(":", 1)
    for itr in introns:
        s, e, _ = itr.split(":")
        rows.append((chrom, int(s), int(e), strand, f"clu_{1000 + k}_{strand}"))
rows.sort(key=lambda r: (r[0], r[1], r[2]))
with open(sys.argv[2], "w") as fh:
    fh.write("chrom\tstart\tend\tstrand\tcluster_id\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
print(f"frozen test clusters: {k} ({len(rows)} introns)", file=sys.stderr)
EOF
gzip -9 -n "$TMP/frozen_clusters.tsv"

python3 "$BIN/xqtl_restricted.py" extract "$HERE/xqtl-restricted-test.tar.gz" "$TMP/site"
python3 "$BIN/xqtl_restricted.py" make --clusters "$TMP/frozen_clusters.tsv.gz" \
    --signals "$TMP/site/signals/mdacc_signals.tsv.gz" --version 1 \
    --out "$OUT/xqtl-restricted-phenotyping-test.tar.gz"
python3 "$BIN/xqtl_restricted.py" verify "$OUT/xqtl-restricted-phenotyping-test.tar.gz"
