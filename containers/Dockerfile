# syntax=docker/dockerfile:1.7
# =============================================================================
# placenta-xqtl-federated: one self-contained, offline image
#   software (atlas pins) + public references + Java runtime + Nextflow + pipeline code
#
# Build from the repository root (containers/build_image.sh does all of this):
#   scripts/fetch_references.sh download --cache /path/refs_cache   # public refs, sha256-verified
#   Rscript book/render_book.R                                      # documentation -> book/_book
#   docker build -f containers/Dockerfile \
#       --build-context refs_cache=/path/refs_cache --build-context docs_html=book/_book \
#       -t xqtl-federated:<ver> .
# The two named build contexts keep the downloaded references and the rendered book out of
# the pipeline-code layer.
# CUDA variant (GPU tensorQTL):
#   docker build ... --build-arg TORCH_INDEX=https://download.pytorch.org/whl/cu132 -t ...:<ver>-cuda
#
# Layers are ordered by how often they change (OS -> environments -> tools -> references
# -> documentation -> pipeline code), so a code-only release re-uploads only the top layers.
# Nothing in the image needs network access at run time.
# =============================================================================

# Base image: Debian 12 (bookworm-slim), pinned by digest in both FROM lines below.

# -----------------------------------------------------------------------------
# Stage 1: IRFinder-S, compiled with the system toolchain (Debian gcc 12 + boost 1.74)
# -----------------------------------------------------------------------------
FROM debian:bookworm-slim@sha256:99c0b25c46d0214d5a2091092e7908907dcaecd6a962cde0ac700acee222b6f3 AS irfinder
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates git g++ make \
        libboost-iostreams-dev zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*
COPY containers/build_irfinder.sh /build/build_irfinder.sh
RUN bash /build/build_irfinder.sh /opt/xqtl

# -----------------------------------------------------------------------------
# Stage 2: operating system + conda environments (the large, stable layers)
# -----------------------------------------------------------------------------
FROM debian:bookworm-slim@sha256:99c0b25c46d0214d5a2091092e7908907dcaecd6a962cde0ac700acee222b6f3 AS envs
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
ARG MICROMAMBA_URL=https://github.com/mamba-org/micromamba-releases/releases/download/2.3.2-0/micromamba-linux-64
ARG MICROMAMBA_SHA256=ffc3cb8d52d4d6b354bdbb979c407719c485392b74e462cbd50811aa88e58f85
ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PYTHONNOUSERSITE=1 \
    PYTHONDONTWRITEBYTECODE=1
# procps: `ps` for Nextflow task metrics; libboost-iostreams/zlib: IRFinder-S runtime
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl procps bzip2 pigz less \
        libboost-iostreams1.74.0 zlib1g libgomp1 \
 && rm -rf /var/lib/apt/lists/* \
 && curl -fsSL -o /usr/local/bin/micromamba "${MICROMAMBA_URL}" \
 && echo "${MICROMAMBA_SHA256}  /usr/local/bin/micromamba" | sha256sum -c - \
 && chmod 0755 /usr/local/bin/micromamba

# Analysis environment: Python 3.11 + R 4.3.1 / Bioconductor 3.18 + genotype tools
COPY containers/env/xqtl.yml /build/env/xqtl.yml
RUN micromamba create -y -p /opt/conda/envs/xqtl -f /build/env/xqtl.yml \
 && micromamba clean -a -y

# Alignment environment (separate: Salmon 1.10.2 needs an ICU that conflicts with R 4.3.1)
COPY containers/env/align.yml /build/env/align.yml
RUN micromamba create -y -p /opt/conda/envs/align -f /build/env/align.yml \
 && micromamba clean -a -y

# tensorQTL stack (Phase-0 validated pins). TORCH_INDEX selects CPU or CUDA torch wheels.
ARG TORCH_INDEX=https://download.pytorch.org/whl/cpu
COPY containers/requirements-pip.txt /build/requirements-pip.txt
RUN /opt/conda/envs/xqtl/bin/python -m pip install --no-cache-dir \
        -r /build/requirements-pip.txt --extra-index-url "${TORCH_INDEX}"

# R packages not on conda (pinned commits / CRAN version). Rhcpp and colocboost compile C++
# with the environment's conda compiler (x86_64-conda-linux-gnu-c++), so the environment's
# bin directory must be on PATH for this step.
COPY containers/install_r_extra.R /build/install_r_extra.R
RUN PATH="/opt/conda/envs/xqtl/bin:${PATH}" Rscript /build/install_r_extra.R \
 && rm -rf /tmp/Rtmp* /root/.cache

# Environment lock files (exact package builds), kept in the image for provenance
RUN mkdir -p /opt/xqtl/env_lock \
 && micromamba env export -p /opt/conda/envs/xqtl --explicit > /opt/xqtl/env_lock/xqtl.conda-explicit.txt \
 && micromamba env export -p /opt/conda/envs/align --explicit > /opt/xqtl/env_lock/align.conda-explicit.txt \
 && /opt/conda/envs/xqtl/bin/python -m pip freeze --all > /opt/xqtl/env_lock/xqtl.pip-freeze.txt \
 && /opt/conda/envs/xqtl/bin/Rscript -e 'ip <- installed.packages(fields = "RemoteSha"); write.table(ip[, c("Package", "Version", "RemoteSha")], sep = "\t", quote = FALSE, row.names = FALSE)' \
        > /opt/xqtl/env_lock/xqtl.R-packages.tsv

# -----------------------------------------------------------------------------
# Stage 3: public reference data, assembled with the image's own pinned tools so the
# result is byte-identical between builds (refs/MANIFEST.tsv checksums are stable)
# -----------------------------------------------------------------------------
FROM envs AS refs
ARG REFS_EXPECTED=report
COPY assets/refs_sources.tsv assets/refs_manifest.tsv /build/refs/
COPY scripts/fetch_references.sh /build/refs/fetch_references.sh
# refs_cache and docs_html are named build contexts (--build-context), not stages
# hadolint ignore=DL3022
COPY --from=refs_cache . /build/refs_cache/
RUN PATH="/opt/conda/envs/xqtl/bin:${PATH}" bash /build/refs/fetch_references.sh build \
        --sources /build/refs/refs_sources.tsv --cache /build/refs_cache --out /opt/xqtl/refs \
        --expected /build/refs/refs_manifest.tsv --expected-mode "${REFS_EXPECTED}" --threads 4

# -----------------------------------------------------------------------------
# Stage 4: final image
# -----------------------------------------------------------------------------
FROM envs AS final
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Intron retention (IRFinder-S v2.0.1)
COPY --from=irfinder /opt/xqtl/IRFinder /opt/xqtl/IRFinder

# Java runtime + Nextflow launcher. These are exported to the host by `bootstrap`, so the
# Nextflow head job needs nothing but this image.
ARG JRE_URL=https://github.com/adoptium/temurin17-binaries/releases/download/jdk-17.0.21%2B9/OpenJDK17U-jre_x64_linux_hotspot_17.0.21_9.tar.gz
ARG JRE_SHA256=9f41bf88c7900c7097444a1bdeb84146ae692bc08842490ab098063c74bf9a8a
ARG NXF_VERSION=26.10.4
ARG NXF_SHA256=2cbda6e7fdbbb1e7a413973eff549653f566ba6bb7bdfbb58eb78d0f831fdec4
RUN curl -fsSL -o /tmp/jre.tar.gz "${JRE_URL}" \
 && echo "${JRE_SHA256}  /tmp/jre.tar.gz" | sha256sum -c - \
 && mkdir -p /opt/xqtl/jre && tar -xzf /tmp/jre.tar.gz --strip-components=1 -C /opt/xqtl/jre \
 && rm /tmp/jre.tar.gz \
 && mkdir -p /opt/xqtl/nextflow/home \
 && curl -fsSL -o /opt/xqtl/nextflow/nextflow \
        "https://github.com/nextflow-io/nextflow/releases/download/v${NXF_VERSION}/nextflow-${NXF_VERSION}-dist" \
 && echo "${NXF_SHA256}  /opt/xqtl/nextflow/nextflow" | sha256sum -c - \
 && chmod 0755 /opt/xqtl/nextflow/nextflow \
 && JAVA_HOME=/opt/xqtl/jre NXF_HOME=/tmp/nxf_home NXF_OFFLINE=true /opt/xqtl/nextflow/nextflow -version \
 && rm -rf /tmp/nxf_home

# Public references (genome, annotation, txrevise, MuSiC, 1kGP panel, editing sites, chain)
COPY --from=refs /opt/xqtl/refs /opt/xqtl/refs

# Offline copy of the documentation site
# hadolint ignore=DL3022
COPY --from=docs_html . /opt/xqtl/docs_html/

# Pipeline code (smallest, most frequently changing layer)
COPY . /opt/xqtl/pipeline/
ENV PATH=/opt/xqtl/pipeline/bin:/opt/conda/envs/align/bin:/opt/conda/envs/xqtl/bin:/opt/xqtl/IRFinder/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    XQTL_ROOT=/opt/xqtl \
    R_PROFILE_USER=/dev/null \
    R_ENVIRON_USER=/dev/null \
    MPLCONFIGDIR=/tmp
RUN chmod -R a+rX /opt/xqtl \
 && chmod 0755 /opt/xqtl/pipeline/containers/entrypoint.sh /opt/xqtl/pipeline/containers/bootstrap.sh \
 && /opt/xqtl/pipeline/bin/xqtl_tool_versions.sh > /opt/xqtl/VERSIONS.tsv \
 && grep -q 'IRFinder-S' /opt/xqtl/VERSIONS.tsv \
 && python3 /opt/xqtl/pipeline/bin/xqtl_refs_check.py --refs-dir /opt/xqtl/refs --quick

ARG XQTL_VERSION=1.0.0-dev
LABEL org.opencontainers.image.title="placenta-xqtl-federated" \
      org.opencontainers.image.description="Federated placental xQTL pipeline: software, public references, Nextflow and pipeline code (offline)" \
      org.opencontainers.image.source="https://github.com/sbresnahan/placenta-xqtl-federated" \
      org.opencontainers.image.licenses="GPL-3.0-only" \
      org.opencontainers.image.version="${XQTL_VERSION}"
RUN echo "${XQTL_VERSION}" > /opt/xqtl/XQTL_VERSION

WORKDIR /tmp
ENTRYPOINT ["/opt/xqtl/pipeline/containers/entrypoint.sh"]
CMD ["help"]
