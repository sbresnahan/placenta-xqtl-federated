#!/usr/bin/env bash
# Test containers/bootstrap.sh + the generated run_xqtl.sh against a directory mock of the
# image (/opt/xqtl layout), without a container runtime:
#   * the pinned Java runtime and Nextflow -dist launcher (URLs and sha256 read from
#     containers/Dockerfile) are downloaded and verified, exactly as in the image build;
#   * bootstrap exports them with the pipeline and docs, and writes run_xqtl.sh;
#   * run_xqtl.sh runs --mode smoke with only the exported Java (no java on PATH) and
#     NXF_OFFLINE, using the local tools (-profile local,dev).
#
#   tests/scripts/test_bootstrap.sh [WORKDIR]       (default: a new temporary directory)
# Needs: curl, the development environment on PATH (scripts/dev_env.sh), book/_book rendered
# (or it uses a placeholder page).
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
work=${1:-$(mktemp -d)}
mkdir -p "$work"; work=$(cd "$work" && pwd)
pass=0 fail=0
ok()   { echo "PASS  $*"; pass=$((pass + 1)); }
bad()  { echo "FAIL  $*"; fail=$((fail + 1)); }
arg()  { sed -n "s/^ARG $1=//p" "$repo/containers/Dockerfile" | head -1; }

# ---- 1. build the mock image tree -------------------------------------------------------
dl="$work/downloads"; mkdir -p "$dl"
jre_url=$(arg JRE_URL) jre_sha=$(arg JRE_SHA256) nxf_ver=$(arg NXF_VERSION) nxf_sha=$(arg NXF_SHA256)
[ -s "$dl/jre.tar.gz" ] || curl -fsSL -o "$dl/jre.tar.gz" "$jre_url"
echo "$jre_sha  $dl/jre.tar.gz" | sha256sum -c --quiet - && ok "JRE pin verifies ($jre_sha)" || bad "JRE checksum"
[ -s "$dl/nextflow" ] || curl -fsSL -o "$dl/nextflow" \
  "https://github.com/nextflow-io/nextflow/releases/download/v${nxf_ver}/nextflow-${nxf_ver}-dist"
echo "$nxf_sha  $dl/nextflow" | sha256sum -c --quiet - && ok "Nextflow $nxf_ver -dist pin verifies" || bad "Nextflow checksum"

root="$work/mock/opt/xqtl"; rm -rf "$work/mock"; mkdir -p "$root/jre" "$root/nextflow/home" "$root/refs"
tar -xzf "$dl/jre.tar.gz" --strip-components=1 -C "$root/jre"
install -m 0755 "$dl/nextflow" "$root/nextflow/nextflow"
# pipeline code as it would be copied into the image (respecting .dockerignore's main exclusions)
mkdir -p "$root/pipeline"
(cd "$repo" && git ls-files -co --exclude-standard | grep -v -e '^book/_book/' -e '^book/_build/' \
   | tar -cf - -T -) | tar -xf - -C "$root/pipeline"
if [ -d "$repo/book/_book" ]; then cp -R "$repo/book/_book" "$root/docs_html"
else mkdir -p "$root/docs_html"; echo "<html>placeholder</html>" > "$root/docs_html/index.html"; fi
cp -R "$repo/tests/fixtures/refs_mini/." "$root/refs/"
"$repo/bin/xqtl_tool_versions.sh" > "$root/VERSIONS.tsv" 2>/dev/null || true
echo "0.0.0-mock" > "$root/XQTL_VERSION"
chmod -R a-w "$root"                     # files inside a SIF are read-only
touch "$work/mock/xqtl-federated_0.0.0-mock.sif"
sif="$work/mock/xqtl-federated_0.0.0-mock.sif"
export XQTL_ROOT="$root"

# ---- 2. entrypoint + bootstrap behaviour --------------------------------------------------
ep="$root/pipeline/containers/entrypoint.sh"
bash "$ep" help | grep -q "version 0.0.0-mock" && ok "entrypoint help" || bad "entrypoint help"
[ "$(bash "$ep" echo passthrough)" = passthrough ] && ok "entrypoint passes other commands through" || bad "passthrough"
bash "$ep" refs-check --quick | tail -1 | grep -q '^OK' && ok "entrypoint refs-check" || bad "refs-check"

dest="$work/site/xqtl"
if (unset APPTAINER_CONTAINER SINGULARITY_CONTAINER; bash "$ep" bootstrap "$dest") 2> "$work/err.txt"; then
  bad "bootstrap without a SIF path should fail"
else grep -q -- "--sif" "$work/err.txt" && ok "bootstrap asks for --sif outside Apptainer" || bad "missing --sif message"; fi

APPTAINER_CONTAINER="$sif" bash "$ep" bootstrap "$dest" > "$work/boot.txt"
for f in run_xqtl.sh pipeline/main.nf jre/bin/java nextflow/nextflow docs/index.html VERSIONS.tsv XQTL_VERSION refs_MANIFEST.tsv; do
  [ -e "$dest/$f" ] || { bad "bootstrap output $f missing"; continue; }
done
[ -x "$dest/run_xqtl.sh" ] && ok "bootstrap wrote run_xqtl.sh and exported pipeline, JRE, Nextflow, docs" || bad "run_xqtl.sh"
grep -q "XQTL_SIF:-$sif" "$dest/run_xqtl.sh" && ok "run_xqtl.sh records the SIF path (APPTAINER_CONTAINER)" || bad "SIF path"
grep -q "NXF_OFFLINE=true" "$dest/run_xqtl.sh" && ok "run_xqtl.sh sets NXF_OFFLINE" || bad "NXF_OFFLINE"
grep -q "stbresnahan@mdanderson.org" "$work/boot.txt" && ok "bootstrap prints the restricted-file request note" || bad "email note"
touch "$dest/pipeline/should_be_replaced" "$dest/cohort.yml"
if APPTAINER_CONTAINER="$sif" bash "$ep" bootstrap "$dest" > /dev/null 2>&1; then bad "re-bootstrap without --force should fail"
else ok "re-bootstrap into a used directory refuses without --force"; fi
APPTAINER_CONTAINER="$sif" bash "$ep" bootstrap "$dest" --force > /dev/null
[ ! -e "$dest/pipeline/should_be_replaced" ] && [ -e "$dest/cohort.yml" ] \
  && ok "--force replaces exported files and keeps user files" || bad "--force behaviour"

# ---- 3. run_xqtl.sh: offline Nextflow with the exported Java only ---------------------------
devpath=$(command -v python3 | xargs dirname)
clean_path="$(dirname "$(command -v STAR 2>/dev/null || echo /nonexistent/x)"):$devpath:$(dirname "$(command -v IRFinder 2>/dev/null || echo /nonexistent/x)"):$repo/bin:/usr/bin:/bin"
if PATH="$clean_path" command -v java > /dev/null 2>&1; then
  echo "NOTE  a java is on the clean PATH ($(PATH="$clean_path" command -v java)); the run still uses JAVA_HOME"
fi
run="$work/site/run"; mkdir -p "$run"
( cd "$run" && env -i HOME="$HOME" USER="${USER:-root}" PATH="$clean_path" TMPDIR="${TMPDIR:-/tmp}" \
    "$dest/run_xqtl.sh" -profile local,dev,smoke \
      --refs_dir "$repo/tests/fixtures/refs_mini" \
      --restricted "$repo/tests/fixtures/xqtl-restricted-test.tar.gz" ) > "$work/smoke.txt" 2>&1 || true
summary="$run/xqtl_smoke/smoke/smoke_summary.tsv"
if grep -qP '^VERDICT\tGO\t' "$summary" 2>/dev/null; then ok "run_xqtl.sh --mode smoke: VERDICT GO"
else bad "smoke run (no VERDICT GO in $summary)"; grep -E 'FAIL|VERDICT' "$summary" 2>/dev/null; tail -30 "$work/smoke.txt"; fi
grep -q "17.0.21" "$run/.nextflow.log" 2>/dev/null && ok "head job ran on the exported Java 17.0.21" || bad "Java version not found in .nextflow.log"
[ -d "$dest/nextflow/home" ] && ok "NXF_HOME is inside the bootstrap directory" || bad "NXF_HOME"
if XQTL_SIF=/nonexistent.sif "$dest/run_xqtl.sh" -profile local 2> "$work/err2.txt"; then bad "missing SIF should fail"
else grep -q "image not found" "$work/err2.txt" && ok "run_xqtl.sh reports a missing image" || bad "missing-image message"; fi

chmod -R u+w "$work/mock"
echo "----"; echo "bootstrap test: $pass passed, $fail failed (work dir $work)"
[ "$fail" -eq 0 ]
