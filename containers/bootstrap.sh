#!/usr/bin/env bash
# bootstrap: export what the Nextflow head job needs from the image to a host directory.
#
#   apptainer run xqtl-federated_<ver>.sif bootstrap DEST [--sif PATH] [--force]
#
# Writes
#   DEST/run_xqtl.sh      launcher: offline Nextflow + this image (NXF_OFFLINE, JAVA_HOME, NXF_HOME, --sif)
#   DEST/pipeline/        pipeline code (main.nf, conf/, bin/, ...), identical to the image copy
#   DEST/jre/             Java 17 runtime for the Nextflow head job
#   DEST/nextflow/        Nextflow launcher (self-contained -dist build) and NXF_HOME
#   DEST/docs/            offline documentation (open DEST/docs/index.html)
#   DEST/VERSIONS.tsv, DEST/XQTL_VERSION, DEST/refs_MANIFEST.tsv
# The reference data stay inside the image; tasks read them at /opt/xqtl/refs.
#
# DEST must be on a host filesystem that is visible inside the container (Apptainer binds
# $HOME and the current directory by default; otherwise add  -B /path).
# XQTL_ROOT (default /opt/xqtl) is overridable for testing against a directory mock.
set -euo pipefail
root=${XQTL_ROOT:-/opt/xqtl}

usage() {
  sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//'
}

dest="" sif="" force=0
while [ $# -gt 0 ]; do
  case "$1" in
    --sif) sif=${2:?--sif needs a path}; shift 2 ;;
    --sif=*) sif=${1#--sif=}; shift ;;
    --force) force=1; shift ;;
    -h|--help) usage; exit 0 ;;
    -*) echo "bootstrap: unknown option $1" >&2; usage >&2; exit 2 ;;
    *) [ -z "$dest" ] || { echo "bootstrap: only one DEST allowed" >&2; exit 2; }; dest=$1; shift ;;
  esac
done
[ -n "$dest" ] || { usage >&2; exit 2; }

for d in pipeline jre nextflow docs_html; do
  [ -d "$root/$d" ] || { echo "bootstrap: $root/$d not found; is this the xqtl-federated image?" >&2; exit 1; }
done
ver=$(cat "$root/XQTL_VERSION" 2>/dev/null || echo unknown)

# The image file on the host: Apptainer/Singularity export its path inside the container.
if [ -z "$sif" ]; then
  sif=${APPTAINER_CONTAINER:-${SINGULARITY_CONTAINER:-}}
fi
if [ -z "$sif" ]; then
  echo "bootstrap: cannot tell where the .sif file is (not running under Apptainer/Singularity)." >&2
  echo "           Re-run with  --sif /path/to/xqtl-federated_${ver}.sif" >&2
  exit 1
fi
case "$sif" in /*) ;; *) sif="$(pwd)/$sif" ;; esac

mkdir -p "$dest"
dest=$(cd "$dest" && pwd)
if [ -f "$dest/XQTL_VERSION" ] && [ "$force" -ne 1 ]; then
  echo "bootstrap: $dest already holds version $(cat "$dest/XQTL_VERSION")." >&2
  echo "           Use a new directory per version, or --force to replace the exported files" >&2
  echo "           (your own files, work directories and results are not touched)." >&2
  exit 1
fi
if [ "$force" -ne 1 ] && [ -n "$(ls -A "$dest" 2>/dev/null)" ]; then
  echo "bootstrap: $dest is not empty; choose an empty directory or use --force." >&2
  exit 1
fi

echo "Exporting xqtl-federated $ver to $dest"
for pair in pipeline:pipeline jre:jre nextflow:nextflow docs_html:docs; do
  src=${pair%%:*} tgt=${pair##*:}
  rm -rf "${dest:?}/$tgt"
  cp -R "$root/$src" "$dest/$tgt"
  chmod -R u+w "$dest/$tgt"          # files in a SIF are read-only
done
mkdir -p "$dest/nextflow/home"
cp "$root/VERSIONS.tsv" "$dest/VERSIONS.tsv" 2>/dev/null || true
cp "$root/refs/MANIFEST.tsv" "$dest/refs_MANIFEST.tsv" 2>/dev/null || true
echo "$ver" > "$dest/XQTL_VERSION"

cat > "$dest/run_xqtl.sh" <<EOF
#!/usr/bin/env bash
# Launcher written by 'bootstrap' (xqtl-federated $ver). Runs Nextflow offline from this
# directory, with every task inside the image below.
#   ./run_xqtl.sh -profile <site> [-params-file cohort.yml] --mode <mode> [-resume]
# Environment overrides: XQTL_SIF (image path), NXF_OPTS (head JVM heap), NXF_HOME.
set -euo pipefail
XQTL_HOME="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
XQTL_SIF="\${XQTL_SIF:-$sif}"
if [ ! -r "\$XQTL_SIF" ]; then
  echo "run_xqtl.sh: image not found or not readable: \$XQTL_SIF" >&2
  echo "             set XQTL_SIF=/path/to/xqtl-federated_${ver}.sif, or re-run bootstrap" >&2
  exit 1
fi
export JAVA_HOME="\$XQTL_HOME/jre"
export NXF_HOME="\${NXF_HOME:-\$XQTL_HOME/nextflow/home}"
export NXF_OFFLINE=true
export NXF_DISABLE_CHECK_LATEST=true
export NXF_OPTS="\${NXF_OPTS:--Xms1g -Xmx4g}"
export PATH="\$JAVA_HOME/bin:\$PATH"
unset NXF_VER NXF_JAVA_HOME JAVA_CMD
exec "\$XQTL_HOME/nextflow/nextflow" run "\$XQTL_HOME/pipeline/main.nf" --sif "\$XQTL_SIF" "\$@"
EOF
chmod 0755 "$dest/run_xqtl.sh"

cat <<EOF

Done. Next steps (see the Setup chapter in $dest/docs/index.html):
  1. Site check (about 10 minutes):
       $dest/run_xqtl.sh -profile <site>,smoke
  2. Place the restricted file xqtl-restricted-<version>.tar.gz and set 'restricted:' in cohort.yml.
     Not received? Request it from Sean T. Bresnahan (stbresnahan@mdanderson.org), stating your
     cohort, the pipeline version ($ver) and the approved data-sharing agreement.
  3. Build the alignment indices (once):
       $dest/run_xqtl.sh -profile <site> -params-file cohort.yml --mode prepare_references
EOF
