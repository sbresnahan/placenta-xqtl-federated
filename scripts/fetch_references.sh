#!/usr/bin/env bash
# Fetch and assemble the public reference data embedded in the image (/opt/xqtl/refs).
#
#   scripts/fetch_references.sh download --cache DIR [--sources TSV] [--allow-unpinned] [--dry-run]
#   scripts/fetch_references.sh build    --cache DIR --out DIR [--sources TSV] [--threads N]
#                                        [--expected TSV --expected-mode strict|report]
#   scripts/fetch_references.sh all      --cache DIR --out DIR ...   (download, then build)
#
# download  needs curl + sha256sum only (runs on the CI runner or any host). Every file is
#           checked against the sha256 pinned in assets/refs_sources.tsv; a mismatch is fatal.
# build     needs bgzip + samtools + python3 (runs inside the image build, so the outputs are
#           made by the image's pinned tools and are byte-identical between builds). Writes
#           OUT/MANIFEST.tsv (path, sha256, bytes, description, source, license) and
#           OUT/README.txt. With --expected, compares MANIFEST.tsv with the release manifest
#           (assets/refs_manifest.tsv): strict = fail on any difference (release builds),
#           report = print differences only (development builds).
#
# assets/refs_sources.tsv columns:
#   id  required  url  sha256  bytes  action  dest  description  source  license
#   required: yes | optional    (optional rows are skipped when unpinned or unavailable)
#   action:   copy | bgzip_fasta (gunzip -> bgzip -> samtools faidx)
#   {ZENODO_RECORD} in a url is replaced by the record id on the '# ZENODO_RECORD=' line.
#   sha256 TO_BE_PINNED = not yet published (see scripts/pin_reference_record.py).
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

cmd=${1:-}; [ $# -gt 0 ] && shift
sources="$here/assets/refs_sources.tsv" cache="" out="" threads=4
expected="" expected_mode=report allow_unpinned=0 dry_run=0
while [ $# -gt 0 ]; do
  case "$1" in
    --sources) sources=$2; shift 2 ;;
    --cache) cache=$2; shift 2 ;;
    --out) out=$2; shift 2 ;;
    --threads) threads=$2; shift 2 ;;
    --expected) expected=$2; shift 2 ;;
    --expected-mode) expected_mode=$2; shift 2 ;;
    --allow-unpinned) allow_unpinned=1; shift ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "fetch_references: unknown argument $1" >&2; exit 2 ;;
  esac
done
case "$cmd" in download|build|all) ;; *) sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//' >&2; exit 2 ;; esac
[ -n "$cache" ] || { echo "fetch_references: --cache is required" >&2; exit 2; }
[ -f "$sources" ] || { echo "fetch_references: sources table $sources not found" >&2; exit 2; }
case "$expected_mode" in strict|report) ;; *) echo "--expected-mode must be strict or report" >&2; exit 2 ;; esac

record=$(sed -n 's/^# *ZENODO_RECORD=\([^[:space:]]*\).*/\1/p' "$sources" | awk 'NR==1')
record=${record:-TO_BE_PINNED}

log() { printf '[fetch_references] %s\n' "$*"; }
die() { printf '[fetch_references] ERROR: %s\n' "$*" >&2; exit 1; }

# rows: id required url sha256 bytes action dest description source license
rows() { grep -v '^#' "$sources" | awk -F'\t' 'NR > 1 && NF >= 10'; }
cache_name() { local u=${1%%\?*}; basename "$u"; }

do_download() {
  mkdir -p "$cache"
  local unpinned=() n=0
  # shellcheck disable=SC2034
  while IFS=$'\t' read -r id required url sha bytes action dest desc src lic; do
    url=${url//\{ZENODO_RECORD\}/$record}
    if [ "$sha" = "TO_BE_PINNED" ] || [[ "$url" == *TO_BE_PINNED* ]]; then
      if [ "$required" = "yes" ]; then unpinned+=("$id"); fi
      log "skip  $id (not yet pinned${required:+, $required})"
      continue
    fi
    local f; f="$cache/$(cache_name "$url")"
    if [ -s "$f" ] && echo "$sha  $f" | sha256sum -c --status - 2>/dev/null; then
      log "ok    $id (cached)"; n=$((n + 1)); continue
    fi
    if [ "$dry_run" -eq 1 ]; then log "would download $id <- $url ($bytes bytes)"; continue; fi
    log "get   $id <- $url"
    if ! curl -fL --retry 5 --retry-delay 10 --connect-timeout 60 -sS -o "$f.part" "$url"; then
      rm -f "$f.part"
      [ "$required" = "yes" ] && die "download failed: $id ($url)"
      log "skip  $id (optional, download failed)"; continue
    fi
    echo "$sha  $f.part" | sha256sum -c --status - \
      || { rm -f "$f.part"; die "checksum mismatch for $id ($url): the source changed or the download is corrupt"; }
    mv "$f.part" "$f"; n=$((n + 1))
  done < <(rows)
  if [ "${#unpinned[@]}" -gt 0 ]; then
    if [ "$allow_unpinned" -eq 1 ]; then
      log "WARNING: required references not yet pinned (development build only): ${unpinned[*]}"
    else
      die "required references are not pinned yet: ${unpinned[*]}.
  These come from the public MDACC Zenodo reference record. Build it with
  scripts/make_public_reference_record.sh, upload it to Zenodo, then run
  scripts/pin_reference_record.py --record <id> --checksums <record>/SHA256SUMS
  (or pass --allow-unpinned for a development image without them)."
    fi
  fi
  log "download complete: $n file(s) in $cache"
}

do_build() {
  [ -n "$out" ] || die "--out is required for build"
  for t in bgzip samtools python3 sha256sum; do command -v "$t" >/dev/null || die "$t not found on PATH"; done
  mkdir -p "$out"
  local built=0
  # shellcheck disable=SC2034  # desc/src/lic are read for column alignment; used by the manifest step
  while IFS=$'\t' read -r id required url sha bytes action dest desc src lic; do
    url=${url//\{ZENODO_RECORD\}/$record}
    local f; f="$cache/$(cache_name "$url")"
    if [ ! -s "$f" ]; then
      [ "$required" = "yes" ] && [ "$allow_unpinned" -eq 0 ] && [ "$sha" != "TO_BE_PINNED" ] \
        && die "$id: $f missing; run 'download' first"
      log "skip  $id (not downloaded)"; continue
    fi
    echo "$sha  $f" | sha256sum -c --status - || die "$id: cached file $f fails its checksum"
    mkdir -p "$out/$(dirname "$dest")"
    case "$action" in
      copy)
        cp "$f" "$out/$dest"; chmod 0644 "$out/$dest" ;;
      bgzip_fasta)
        log "bgzip $id -> $dest"
        gzip -dc "$f" | bgzip -@ "$threads" -c > "$out/$dest"
        samtools faidx "$out/$dest"
        chmod 0644 "$out/$dest" "$out/$dest".* ;;
      *) die "$id: unknown action $action" ;;
    esac
    built=$((built + 1))
    log "built $id -> $dest"
  done < <(rows)
  [ "$built" -gt 0 ] || die "nothing was built"

  # Provenance check: does the public genome equal the one the atlas was built on?
  # (genome/atlas_genome_contigs.tsv is written by make_public_reference_record.sh)
  local fa; fa=$(rows | awk -F'\t' '$6 == "bgzip_fasta" && !n++ {print $7}')  # no early exit: SIGPIPE under pipefail
  if [ -n "$fa" ] && [ -s "$out/genome/atlas_genome_contigs.tsv" ]; then
    local tmpd; tmpd=$(mktemp -d)
    samtools dict "$out/$fa" | awk -F'\t' '$1=="@SQ"{for(i=2;i<=NF;i++){split($i,a,":"); v[a[1]]=substr($i,4)} print v["SN"]"\t"v["LN"]"\t"v["M5"]}' \
      | sort > "$tmpd/public_contigs.tsv"
    grep -v '^#' "$out/genome/atlas_genome_contigs.tsv" | awk -F'\t' 'NR>1' | sort > "$tmpd/atlas_contigs.tsv"
    if cmp -s "$tmpd/public_contigs.tsv" "$tmpd/atlas_contigs.tsv"; then
      log "genome: identical to the atlas genome (all contig names, lengths and MD5s match)"
    else
      log "WARNING: genome differs from the atlas genome in these contigs (name, length, md5):"
      diff "$tmpd/atlas_contigs.tsv" "$tmpd/public_contigs.tsv" | grep '^[<>]' | sed -n '1,20p' || true
    fi
    rm -rf "$tmpd"
  fi

  python3 - "$sources" "$out" "$record" <<'PY'
import csv, hashlib, os, sys
sources, out, record = sys.argv[1:4]
meta = {}
with open(sources) as fh:
    lines = [l for l in fh if not l.startswith("#")]
for r in csv.DictReader(lines, delimiter="\t"):
    meta[r["dest"]] = r
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()
rows = []
for d, _, fs in os.walk(out):
    for f in fs:
        p = os.path.join(d, f)
        rel = os.path.relpath(p, out)
        if rel in ("MANIFEST.tsv", "README.txt"):
            continue
        m = meta.get(rel) or next((v for k, v in meta.items() if rel.startswith(k + ".")), None)
        desc = m["description"] if m else "derived file"
        if m and rel != m["dest"]:
            desc = f"index of {m['dest']}"
        rows.append([rel, sha(p), str(os.path.getsize(p)), desc,
                     m["source"] if m else "", m["license"] if m else ""])
rows.sort()
with open(os.path.join(out, "MANIFEST.tsv"), "w") as fh:
    fh.write("path\tsha256\tbytes\tdescription\tsource\tlicense\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
with open(os.path.join(out, "README.txt"), "w") as fh:
    fh.write("Public reference data for placenta-xqtl-federated.\n"
             "Provenance, versions and licenses: docs/REFERENCES.md. Checksums: MANIFEST.tsv.\n"
             f"MDACC public reference record: Zenodo {record}\n")
print(f"[fetch_references] MANIFEST.tsv: {len(rows)} files")
PY

  if [ -n "$expected" ]; then
    if [ ! -s "$expected" ] || [ "$(grep -vc '^#' "$expected")" -le 1 ]; then
      log "expected manifest $expected has no entries yet (first build): commit the generated"
      log "MANIFEST.tsv as assets/refs_manifest.tsv before tagging a release"
      if [ "$expected_mode" = strict ]; then die "release builds need a pinned assets/refs_manifest.tsv"; fi
    elif diff <(grep -v '^#' "$expected") "$out/MANIFEST.tsv" > "$cache/.manifest.diff" 2>&1; then
      log "MANIFEST.tsv matches the release manifest $expected"
    else
      log "MANIFEST.tsv differs from $expected:"; head -40 "$cache/.manifest.diff"
      if [ "$expected_mode" = strict ]; then die "embedded references do not match the release manifest"; fi
    fi
  fi
}

case "$cmd" in
  download) do_download ;;
  build) do_build ;;
  all) do_download; [ "$dry_run" -eq 1 ] || do_build ;;
esac
