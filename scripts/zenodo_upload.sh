#!/usr/bin/env bash
# zenodo_upload.sh -- archive release files (the .sif and its checksum) on Zenodo.
# Used by .github/workflows/image.yml on a release tag; can be run by hand.
#
#   ZENODO_TOKEN=... scripts/zenodo_upload.sh --version 1.0.0 [--concept ID] [--publish]
#                     [--sandbox] [--dry-run] FILE [FILE ...]
#     --concept ID   an existing deposition (any version) of the image record: a new version
#                    of it is created, so every release shares one concept DOI.
#                    Without it a new record is created.
#     --publish      publish the record (irreversible). Default: leave it as a draft to review
#                    and publish on zenodo.org.
#     --sandbox      use sandbox.zenodo.org (testing)
# Prints the deposition URL (and the DOI when published).
# Needs: curl, python3.
set -euo pipefail

api=https://zenodo.org/api; version=""; concept=""; publish=no; dry=no; files=()
while [ $# -gt 0 ]; do
  case "$1" in
    --version) version=$2; shift 2 ;;
    --concept) concept=$2; shift 2 ;;
    --publish) publish=yes; shift ;;
    --sandbox) api=https://sandbox.zenodo.org/api; shift ;;
    --dry-run) dry=yes; shift ;;
    -h|--help) sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    -*) echo "zenodo_upload: unknown option $1" >&2; exit 2 ;;
    *) files+=("$1"); shift ;;
  esac
done
[ -n "$version" ] || { echo "zenodo_upload: --version is required" >&2; exit 2; }
[ "${#files[@]}" -gt 0 ] || { echo "zenodo_upload: no files given" >&2; exit 2; }
for f in "${files[@]}"; do [ -f "$f" ] || { echo "zenodo_upload: $f not found" >&2; exit 2; }; done

metadata=$(python3 - "$version" <<'PY'
import json, sys
v = sys.argv[1]
print(json.dumps({"metadata": {
    "title": f"placenta-xqtl-federated {v}: offline Apptainer image",
    "upload_type": "software",
    "version": v,
    "description": (
        "Self-contained, offline Apptainer/Singularity image (.sif) of the placenta-xqtl-federated "
        "Nextflow pipeline: software environment, public reference data, Java runtime, Nextflow and "
        "pipeline code. Source and documentation: https://github.com/sbresnahan/placenta-xqtl-federated. "
        "The restricted file (frozen splicing clusters and the MDACC signal list) is NOT included; "
        "approved sites request it from Sean T. Bresnahan (stbresnahan@mdanderson.org)."),
    "creators": [{"name": "Bresnahan, Sean T.", "affiliation": "The University of Texas MD Anderson Cancer Center"}],
    "license": "GPL-3.0",
    "access_right": "open",
    "related_identifiers": [{"identifier": "https://github.com/sbresnahan/placenta-xqtl-federated",
                             "relation": "isSupplementTo", "scheme": "url"}],
}}))
PY
)

if [ "$dry" = yes ]; then
  echo "API: $api"; echo "concept: ${concept:-<new record>}"; echo "publish: $publish"
  echo "files:"; for f in "${files[@]}"; do echo "  $f ($(stat -c %s "$f") bytes)"; done
  echo "metadata: $metadata"; exit 0
fi
[ -n "${ZENODO_TOKEN:-}" ] || { echo "zenodo_upload: ZENODO_TOKEN is not set" >&2; exit 2; }
auth=(-H "Authorization: Bearer $ZENODO_TOKEN")
jget() {  # jget key [key ...]: print a nested value of the JSON on stdin
  python3 -c 'import json,sys
d=json.load(sys.stdin)
for k in sys.argv[1:]: d=d[k]
print(d)' "$@"
}

if [ -n "$concept" ]; then
  # new version of an existing record; the response links to the new draft
  draft=$(curl -fsS "${auth[@]}" -X POST "$api/deposit/depositions/$concept/actions/newversion" \
          | jget links latest_draft)
  dep=$(curl -fsS "${auth[@]}" "$draft")
  # a new version starts with the previous version's files: remove them
  echo "$dep" | python3 -c 'import json,sys; [print(f["links"]["self"]) for f in json.load(sys.stdin).get("files", [])]' \
    | while read -r link; do curl -fsS "${auth[@]}" -X DELETE "$link" > /dev/null; done
else
  dep=$(curl -fsS "${auth[@]}" -H "Content-Type: application/json" -X POST "$api/deposit/depositions" -d '{}')
fi
id=$(echo "$dep" | jget id)
bucket=$(echo "$dep" | jget links bucket)
html=$(echo "$dep" | jget links html)

for f in "${files[@]}"; do
  echo "uploading $(basename "$f") ($(stat -c %s "$f") bytes)"
  curl -fsS "${auth[@]}" --upload-file "$f" "$bucket/$(basename "$f")" > /dev/null
done
curl -fsS "${auth[@]}" -H "Content-Type: application/json" -X PUT "$api/deposit/depositions/$id" \
     -d "$metadata" > /dev/null

if [ "$publish" = yes ]; then
  doi=$(curl -fsS "${auth[@]}" -X POST "$api/deposit/depositions/$id/actions/publish" | jget doi)
  echo "published: $html  DOI $doi"
else
  echo "draft ready for review (not published): $html"
fi
