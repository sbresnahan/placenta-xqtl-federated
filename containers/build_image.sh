#!/usr/bin/env bash
# build_image.sh -- build the self-contained offline image (used locally and by
# .github/workflows/image.yml).
#
#   containers/build_image.sh [options]
#     --tag NAME           image tag            (default xqtl-federated:<version from nextflow.config>)
#     --version V          XQTL_VERSION baked into the image (default: nextflow.config manifest version)
#     --refs-cache DIR     downloaded public references   (default ./refs_cache)
#     --variant cpu|cuda   tensorQTL torch build           (default cpu; cuda uses the cu132 wheels)
#     --release            strict mode: every reference row must be pinned and the processed
#                          references must match assets/refs_manifest.tsv exactly
#     --skip-download      use --refs-cache as is (no fetch_references.sh download)
#     --skip-book          use book/_book as is (no render)
#     --dry-run            print the commands only
#   Any argument after `--` is passed to docker build (e.g. -- --push --platform linux/amd64).
#
# Steps: (1) download + sha256-verify public references into the cache, (2) render the
# documentation book if absent, (3) docker build with the two named build contexts.
# Needs: docker with BuildKit (buildx), curl, sha256sum; Rscript + bookdown + pandoc to render the book.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)

die() { echo "build_image.sh: $*" >&2; exit 1; }

version=$(sed -n "s/^[[:space:]]*version[[:space:]]*=[[:space:]]*'\([^']*\)'.*/\1/p" "$repo/nextflow.config" | awk 'NR==1')
tag=""; cache="$repo/refs_cache"; variant=cpu; release=no
skip_download=no; skip_book=no; dry=no; extra=()
while [ $# -gt 0 ]; do
  case "$1" in
    --tag) tag=$2; shift 2 ;;
    --version) version=$2; shift 2 ;;
    --refs-cache) cache=$2; shift 2 ;;
    --variant) variant=$2; shift 2 ;;
    --release) release=yes; shift ;;
    --skip-download) skip_download=yes; shift ;;
    --skip-book) skip_book=yes; shift ;;
    --dry-run) dry=yes; shift ;;
    --) shift; extra=("$@"); break ;;
    -h|--help) sed -n '2,19p' "$0"; exit 0 ;;
    *) die "unknown option $1 (see --help)" ;;
  esac
done
[ -n "$version" ] || die "could not read the version from nextflow.config; pass --version"

case "$variant" in
  cpu)  torch_index=https://download.pytorch.org/whl/cpu; suffix="" ;;
  cuda) torch_index=https://download.pytorch.org/whl/cu132; suffix="-cuda" ;;
  *) die "--variant must be cpu or cuda" ;;
esac
[ -n "$tag" ] || tag="xqtl-federated:${version}${suffix}"

run() { echo "+ $*"; if [ "$dry" = no ]; then "$@"; fi; }

# ---- 1. public references --------------------------------------------------------------
if [ "$skip_download" = no ]; then
  dl=("$repo/scripts/fetch_references.sh" download --cache "$cache")
  [ "$release" = yes ] || dl+=(--allow-unpinned)
  run "${dl[@]}"
fi
if [ "$dry" = no ] && [ ! -d "$cache" ]; then die "reference cache $cache does not exist"; fi

# ---- 2. documentation --------------------------------------------------------------------
if [ "$skip_book" = no ] && [ ! -f "$repo/book/_book/index.html" ]; then
  command -v Rscript > /dev/null || die "book/_book is missing and Rscript is not available to render it"
  run Rscript "$repo/book/render_book.R"
fi
if [ "$dry" = no ] && [ ! -f "$repo/book/_book/index.html" ]; then die "book/_book/index.html missing"; fi

# ---- 3. image ------------------------------------------------------------------------------
refs_expected=report
[ "$release" = yes ] && refs_expected=strict
build=(docker buildx build --load
  -f "$repo/containers/Dockerfile"
  --build-context "refs_cache=$cache"
  --build-context "docs_html=$repo/book/_book"
  --build-arg "XQTL_VERSION=$version"
  --build-arg "TORCH_INDEX=$torch_index"
  --build-arg "REFS_EXPECTED=$refs_expected"
  -t "$tag")
# an explicit --push replaces --load
for a in "${extra[@]}"; do
  if [ "$a" = "--push" ]; then build=("${build[@]/--load/}"); fi
done
build=("${build[@]}" "${extra[@]}" "$repo")
# drop empty elements left by the substitution above
clean=(); for a in "${build[@]}"; do [ -n "$a" ] && clean+=("$a"); done
run "${clean[@]}"
echo "built $tag (version $version, variant $variant, references $refs_expected)"
if [ "$dry" = no ] && [ "${#extra[@]}" -eq 0 ]; then
  echo "verify:  docker run --rm $tag versions && docker run --rm $tag refs-check"
fi
