# Changelog

All notable changes to this pipeline are recorded here. Versions follow [semantic versioning](https://semver.org/): a change in the export format (files A–J) or in the restricted-file version is a major version change.

## 1.0.0-dev (unreleased)

### Added
- Repository skeleton: `main.nf` with run modes, `nextflow.config`, resource classes (`conf/base.config`), executor profiles for LSF, SLURM, PBS Pro, Torque, SGE and local runs, Apptainer/Singularity/Docker profiles, the built-in `seadragon` profile and site templates.
- `--mode smoke`: site validation (12 checks with a GO / NO-GO verdict), folded in from the Phase 0 seadragon kit.
- Input validation (`bin/xqtl_validate.py`) against `assets/schema/cohort.schema.json`, run at the start of every mode, with `--validate_only`.
- Restricted-file tooling (`bin/xqtl_restricted.py`): byte-reproducible archive, version and checksum verification.
- Embedded-reference manifest check (`bin/xqtl_refs_check.py`) and tool version report (`bin/xqtl_tool_versions.sh`).
- Software environment definition (`containers/env/`, `containers/requirements-pip.txt`, `containers/install_r_extra.R`, `containers/build_irfinder.sh`) and `scripts/dev_env.sh` for local development.
- Documentation: collaboration brief, HPC configuration guide, seadragon validation; bookdown site builder (`book/`) and GitHub Pages workflow (`.github/workflows/docs.yml`).
- Offline image recipe (`containers/Dockerfile`): software environments, IRFinder-S, public references, Java 17 runtime, Nextflow 26.10.4, pipeline code and the rendered documentation in one image; entry point with `bootstrap`, `versions`, `refs-check` and `docs` commands.
- `bootstrap` (`containers/bootstrap.sh`): exports the pipeline, Java runtime, Nextflow and documentation to the host and writes `run_xqtl.sh` (offline Nextflow with the image path, `JAVA_HOME` and `NXF_HOME` set).
- Public reference data: sources pinned by sha256 (`assets/refs_sources.tsv`), download and deterministic processing (`scripts/fetch_references.sh`), release checksum list (`assets/refs_manifest.tsv`).
- MDACC packaging: `scripts/make_public_reference_record.sh` (Zenodo reference record, including the pre-pruned 1000 Genomes PCA panel), `scripts/pin_reference_record.py`, `scripts/make_restricted_bundle.sh` (frozen splicing clusters merged across ancestries, MDACC signal list), settings in `scripts/packaging/seadragon.env`.
- Image build and release: `containers/build_image.sh`; workflows `ci.yml` (unit tests, profile resolution, `nextflow lint`, shellcheck, hadolint, actionlint, image self-checks) and `image.yml` (development image on `main`; on a version tag: CPU and CUDA images, `.sif` pushed to GHCR as an OCI artifact, Zenodo archive via `scripts/zenodo_upload.sh`, release assets).
- Documentation: setup guide and reference-data provenance page.

### Notes
- The CUDA image variant uses PyTorch `cu132` wheels (no older CUDA build exists for the pinned PyTorch release); GPU nodes need a driver supporting CUDA 13.2.

### Not yet available
- `prepare_references`, `imputation_prep`, `site`, `phenotyping_only`, `export` and `central` modes stop with a message until their milestones are complete.
