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

### Not yet available
- `prepare_references`, `imputation_prep`, `site`, `phenotyping_only`, `export` and `central` modes stop with a message until their milestones are complete.
