#!/usr/bin/env bash
# Fallback publisher: render the book locally and force-push it to the orphan
# `bookdown` branch (same result as .github/workflows/docs.yml).
#   book/publish_pages.sh [remote]      (default remote: origin)
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
remote=${1:-origin}
cd "$repo"
Rscript book/render_book.R
sha=$(git rev-parse --short HEAD)
tmp=$(mktemp -d)
cp -r book/_book/. "$tmp/"
touch "$tmp/.nojekyll"
cd "$tmp"
git init -q -b bookdown
git add -A
git commit -q -m "docs: $sha"
git push -f "$(git -C "$repo" remote get-url "$remote")" bookdown:bookdown
rm -rf "$tmp"
echo "published book for $sha to $remote/bookdown"
