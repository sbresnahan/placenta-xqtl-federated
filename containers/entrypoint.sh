#!/usr/bin/env bash
# Image entry point (`apptainer run <sif> ...`, `docker run <image> ...`).
#   bootstrap DEST [options]   export launcher + pipeline + docs to the host; writes run_xqtl.sh
#   versions                   tool versions in this image (TSV)
#   refs-check [--quick]       verify the embedded references against this release's checksums
#   docs                       where the offline documentation is
#   help                       this message
# Anything else is executed as a command, so Nextflow tasks (`/bin/bash .command.run`) and
# `apptainer exec` work as usual.
set -euo pipefail
root=${XQTL_ROOT:-/opt/xqtl}
ver=$(cat "$root/XQTL_VERSION" 2>/dev/null || echo unknown)

usage() {
  cat <<EOF
placenta-xqtl-federated image, version $ver

Usage: apptainer run xqtl-federated_${ver}.sif <command>

Commands
  bootstrap DEST     export the pipeline, Java runtime, Nextflow and documentation to DEST
                     on the host and write DEST/run_xqtl.sh (see 'bootstrap --help')
  versions           software versions in this image
  refs-check         verify the embedded reference data (add --quick for sizes only)
  docs               location of the offline documentation
  help               this message

Then run the pipeline from the host:  DEST/run_xqtl.sh -profile <site> ... --mode <mode>
Documentation: https://sbresnahan.github.io/placenta-xqtl-federated/
EOF
}

cmd=${1:-help}
case "$cmd" in
  bootstrap)
    shift; exec bash "$root/pipeline/containers/bootstrap.sh" "$@" ;;
  versions)
    cat "$root/VERSIONS.tsv" ;;
  refs-check)
    shift
    exec python3 "$root/pipeline/bin/xqtl_refs_check.py" --refs-dir "$root/refs" \
        --expected "$root/pipeline/assets/refs_manifest.tsv" "$@" ;;
  docs)
    echo "Offline documentation: $root/docs_html/index.html (exported to DEST/docs by bootstrap)" ;;
  help|-h|--help)
    usage ;;
  *)
    exec "$@" ;;
esac
