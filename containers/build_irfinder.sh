#!/usr/bin/env bash
# Build IRFinder-S v2.0.1 (intron retention) from source into PREFIX/IRFinder and link
# PREFIX/bin/IRFinder. The upstream release asset is itself a Singularity image, which
# cannot run inside our container, so we compile the three C/C++ helpers instead.
# Two one-line portability patches for current compilers (no behavioural change):
#   * CoverageBlocks.cpp uses uint64_t without <cstdint>
#   * winflat_with_beta.c uses strcmp without <string.h>
# Set BOOST_LIBDIR when boost is outside the system library path.
# Needs: git, make, g++/gcc (<= 13 tested), boost-iostreams and zlib headers/libs.
#   containers/build_irfinder.sh /opt/xqtl            (image)
#   containers/build_irfinder.sh .xqtl-env/align      (dev environment)
set -euo pipefail
prefix=$(realpath -m "${1:?usage: build_irfinder.sh PREFIX}")
tag=v2.0.1
commit=14f053aaa65c2a56df57cb866b241542710898d2
src=$(mktemp -d)
git clone -q --branch "$tag" https://github.com/RitchieLabIGH/IRFinder.git "$src/IRFinder"
cd "$src/IRFinder"
[ "$(git rev-parse HEAD)" = "$commit" ] || { echo "IRFinder $tag is not commit $commit" >&2; exit 1; }
sed -i '1i #include <cstdint>' src/irfinder/src/ReadBlock/CoverageBlocks.cpp
sed -i '1i #include <string.h>' src/winflat/winflat_with_beta.c
jobs=$(nproc 2>/dev/null || echo 2)
(cd src/trim && make clean >/dev/null && make >/dev/null)
(cd src/winflat && make clean >/dev/null && make >/dev/null)
libs="-lboost_iostreams"   # embed an rpath when boost is not in a system location (dev envs)
[ -n "${BOOST_LIBDIR:-}" ] && libs="-L$BOOST_LIBDIR -lboost_iostreams -Wl,-rpath,$BOOST_LIBDIR"
(cd src/irfinder/Release && make clean >/dev/null && make -j"$jobs" LIBS="$libs" >/dev/null)
mkdir -p bin/util "$prefix/bin"
cp src/trim/trim src/winflat/winflat src/irfinder/Release/irfinder bin/util/
chmod -R a+rx bin
rm -rf "$prefix/IRFinder" && mkdir -p "$prefix/IRFinder"
cp -r bin REF LICENSE CHANGELOG.md "$prefix/IRFinder/"
ln -sf "$prefix/IRFinder/bin/IRFinder" "$prefix/bin/IRFinder"
rm -rf "$src"
"$prefix/bin/IRFinder" --help >/dev/null 2>&1 || true
echo "IRFinder-S $tag installed in $prefix/IRFinder"
