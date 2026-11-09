#!/usr/bin/env bash
# make_restricted_test.sh -- rebuild tests/fixtures/xqtl-restricted-test.tar.gz (version 1), the
# restricted file used by -profile test. Developer tool; the archive is committed.
#
# The frozen splicing clusters are kept from the current archive. The MDACC signal list is
# written here; each pair is chosen to hit one G x C / G x E tier-2 flag in the site fixture
# (tests/fixtures/site, planted effects in TRUTH.tsv):
#   ENSGFIX0000011       planted G x C lead, rank 1          -> source both (site lead == MDACC lead)
#   ENSGFIX0000031 (+)   planted G x E splicing lead, rank 1 -> source both
#   ENSGFIX0000035       null gene, common variant, rank 1   -> source mdacc, tested
#   ENSGFIX0000006       rank 1, not the site lead variant   -> site and mdacc pairs flagged separately
#   ENSGFIX0000011       second variant, rank 2              -> not a G x C pair; phenotype in tier 2
#   ENSGFIX0000003       variant not in the genotypes        -> untested, variant_absent
#   ENSG00000100000      phenotype not in the stratum        -> untested, phenotype_absent
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
bin="$here/../../bin"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

python3 "$bin/xqtl_restricted.py" extract "$here/xqtl-restricted-test.tar.gz" "$tmp/old" > /dev/null
printf '%s\n' \
  $'phenotype_id\tvariant_id\tmodality\tancestry\tsignal_rank' \
  $'ENSGFIX0000011\tchr22:19730097:A:G\texpression\tEUR\t1' \
  $'ENSGFIX0000031__chr22:20762038:20762738:clu_31_+\tchr22:20636375:C:T\tsplicing\tEAS\t1' \
  $'ENSGFIX0000035\tchr22:20918446:T:C\texpression\tEAS\t1' \
  $'ENSGFIX0000006\tchr22:19505121:G:A\texpression\tEAS\t1' \
  $'ENSGFIX0000011\tchr22:19734243:G:A\texpression\tEUR\t2' \
  $'ENSGFIX0000003\tchr22:20000001:A:C\texpression\tEUR\t1' \
  $'ENSG00000100000\tchr22:17000500:A:G\texpression\tEUR\t1' \
  | gzip -9 -n > "$tmp/mdacc_signals.tsv.gz"

python3 "$bin/xqtl_restricted.py" make --clusters "$tmp/old/splicing/frozen_clusters.tsv.gz" \
    --signals "$tmp/mdacc_signals.tsv.gz" --version 1 --out "$tmp/xqtl-restricted-test.tar.gz"
python3 "$bin/xqtl_restricted.py" verify "$tmp/xqtl-restricted-test.tar.gz"
mv "$tmp/xqtl-restricted-test.tar.gz" "$here/xqtl-restricted-test.tar.gz"
