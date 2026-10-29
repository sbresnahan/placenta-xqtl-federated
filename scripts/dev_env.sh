#!/usr/bin/env bash
# Build a local environment equivalent to the image's software layer, for development
# and tests on machines without a container runtime. Not used at sites.
#   scripts/dev_env.sh [PREFIX]      (default ./.xqtl-env)
# Then:  export PATH="$PREFIX/align/bin:$PREFIX/xqtl/bin:$PREFIX/irfinder/bin:$PWD/bin:$PATH"
set -euo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
prefix=${1:-$here/.xqtl-env}
mm=$(command -v micromamba || command -v mamba || command -v conda)
"$mm" create -y -p "$prefix/xqtl"  -f "$here/containers/env/xqtl.yml"
"$mm" create -y -p "$prefix/align" -f "$here/containers/env/align.yml"
"$prefix/xqtl/bin/python" -m pip install --no-cache-dir -r "$here/containers/requirements-pip.txt" \
    --extra-index-url https://download.pytorch.org/whl/cpu
"$prefix/xqtl/bin/Rscript" "$here/containers/install_r_extra.R"
# IRFinder-S from source (needs a C++ compiler <= gcc 13 on the host, plus boost/zlib from conda)
# (own environment: Salmon 1.10.2 pins an older boost, so boost must not go into align/)
"$mm" create -y -p "$prefix/irfinder" -c conda-forge "libboost-devel=1.84" zlib make
CPATH="$prefix/irfinder/include" BOOST_LIBDIR="$prefix/irfinder/lib" "$here/containers/build_irfinder.sh" "$prefix/irfinder"
echo "export PATH=\"$prefix/align/bin:$prefix/xqtl/bin:$prefix/irfinder/bin:$here/bin:\$PATH\""
