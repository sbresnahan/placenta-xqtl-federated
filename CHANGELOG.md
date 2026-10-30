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
- `--mode site` (genotypes to export, with ready phenotype inputs from `phenotypes_dir`):
  - post-imputation QC and merge (`xqtl_geno_qc.sh`, `xqtl_geno_merge.sh`);
  - ancestry strata from PC-AiR against the 1000 Genomes panel, with QDA-posterior confidence (`xqtl_ancestry.sh`, `xqtl_ancestry.R`);
  - per-stratum genotype QC and PCs (`xqtl_stratum_geno.sh`);
  - normalization (`xqtl_normalize.R`), HCP factors (`xqtl_hcp.R`), covariates (`xqtl_covariates.py`) and k selection (`xqtl_select_k.py`);
  - tensorQTL permutation, chunked conditionally independent and nominal passes (`xqtl_tensorqtl.py`, `xqtl_merge_parquet.py`);
  - the combined-modality arm (`xqtl_combine_modalities.py`) and in-sample LD (`xqtl_ld.py`).
- Export files A–F, I and J with a mandatory leak scan (`xqtl_export.py`); `--mode export` rebuilds the export from a finished site run. Specification: `docs/EXPORT_SPEC.md`.
- Test fixture `tests/fixtures/site/`: 130 1000 Genomes samples (chr22) with synthetic phenotypes and planted cis effects (`make_site_fixture.sh`). `tests/scripts/test_site.sh` checks signal recovery, the AFR below-minimum stratum and the leak scan end to end, and runs in CI with the development image.

### Fixed
- Image build: `containers/install_r_extra.R` now runs with the environment's `bin` on `PATH` so Rhcpp and colocboost find the conda C++ compiler; it also stops early if the bioconda `GenomeInfoDbData` post-link download failed.
- CI: base image pinned by digest directly in the `FROM` lines (hadolint DL3006); `actions/checkout` and `actions/setup-java` v5.
- Leak scan: text columns were skipped under pandas >= 3 (new `str` dtype), so sample IDs in string columns of parquet/TSV exports could pass; the scanner now checks every non-numeric, non-boolean, non-datetime column.
- CI unit tests: dependencies pinned to the image versions (`tests/requirements-ci.txt`), plus a second run with the latest pandas/pyarrow/numpy.
- `image.yml`: references are downloaded and cached in their own steps (cached even when the build fails), the book is rendered in its own step, and a failed `docker build` posts its last log lines as a public annotation and uploads `build.log`.

### Notes
- The CUDA image variant uses PyTorch `cu132` wheels (no older CUDA build exists for the pinned PyTorch release); GPU nodes need a driver supporting CUDA 13.2.

### Changed
- `image.yml` no longer runs on pushes or tags until M6 is complete; manual runs still work.
- `conf/test.config` uses the site fixture and re-applies `resourceLimits` (2 CPU, 6 GB), because `base.config` evaluates the limits before profiles set `params.max_*`.

### Not yet available
- `prepare_references`, `imputation_prep`, `phenotyping_only` and `central` modes stop with a message until their milestones are complete. `--mode site` currently needs `phenotypes_dir` (FASTQ phenotyping arrives with M4), and G×C/GxE (files G and H) arrive with M5.
