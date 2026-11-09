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
- Documentation: setup guide and reference-data provenance page; configuration reference (every `cohort.yml` setting and command-line parameter), end-to-end walkthrough, troubleshooting guide, and the central-analysis page with the limits of the shared files. A unit test checks that every schema setting is documented.
- `--mode site` (genotypes to export, with ready phenotype inputs from `phenotypes_dir`):
  - post-imputation QC and merge (`xqtl_geno_qc.sh`, `xqtl_geno_merge.sh`);
  - ancestry strata from PC-AiR against the 1000 Genomes panel, with QDA-posterior confidence (`xqtl_ancestry.sh`, `xqtl_ancestry.R`);
  - per-stratum genotype QC and PCs (`xqtl_stratum_geno.sh`);
  - normalization (`xqtl_normalize.R`), HCP factors (`xqtl_hcp.R`), covariates (`xqtl_covariates.py`) and k selection (`xqtl_select_k.py`);
  - tensorQTL permutation, chunked conditionally independent and nominal passes (`xqtl_tensorqtl.py`, `xqtl_merge_parquet.py`);
  - the combined-modality arm (`xqtl_combine_modalities.py`) and in-sample LD (`xqtl_ld.py`).
- Export files A–F, I and J with a mandatory leak scan (`xqtl_export.py`); `--mode export` rebuilds the export from a finished site run. Specification: `docs/EXPORT_SPEC.md`.
- `--mode prepare_references`: site indices in `site_refs` (genome with Picard dictionary, annotation derivatives, STAR, Salmon expression and txrevise indices, IRFinder-S reference, REDIportal sites, `MANIFEST.tsv`).
- Phenotyping from FASTQs (`modules/local/phenotyping.nf`, `workflows/phenotyping.nf`), now run by `--mode site` when `phenotypes_dir` is not set: STAR, Salmon, regtools, featureCounts, IRFinder-S, RNA-editing levels and Picard QC per sample; the nine modality BEDs, count BEDs, `qc_metrics.tsv` and MuSiC cell-type proportions per cohort. Splicing uses the frozen clusters from the restricted file plus site-specific novel clusters. See `docs/PHENOTYPING.md`.
- `--mode phenotyping_only`: phenotype BEDs and normalized matrices without genotypes; nothing is exported.
- File E has a `novel` column (splicing junctions in site-specific clusters); `isoform_expression.counts.bed.gz` is accepted in `phenotypes_dir`.
- Validator checks that `site_refs` is complete when phenotyping runs and warns if it was built for a different read length; `xqtl_restricted.py verify` checks the frozen cluster ID format.
- Test fixture `tests/fixtures/phenotyping/` and `tests/fixtures/refs_chr22/`: 8 public SNUH placenta RNA-seq runs (PRJNA820329) cut to chr22, with GENCODE v45 chr22 and txrevise events as a stand-in annotation, plus chr22 cuts of the REDIportal sites and the MuSiC single-cell reference (`tests/fixtures/extract_atlas_refs_chr22.sh`). `tests/scripts/compare_atlas_phenotypes.R` compares a phenotyping run with atlas phenotypes for the same samples (results in `docs/PHENOTYPING.md`). `tests/scripts/test_phenotyping.sh` runs `prepare_references` and `phenotyping_only` end to end.
- Test fixture `tests/fixtures/site/`: 130 1000 Genomes samples (chr22) with synthetic phenotypes and planted cis effects (`make_site_fixture.sh`). `tests/scripts/test_site.sh` checks signal recovery, the AFR below-minimum stratum and the leak scan end to end, and runs in CI with the development image.
- Interaction analyses in `--mode site` (`modules/local/interactions.nf`; [Interaction analyses](docs/INTERACTIONS.md)):
  - signal test sets per stratum × modality (`xqtl_signals.py`): site leads plus MDACC rank-1 pairs for G×C, site-lead plus MDACC phenotypes for GxE tier 2, with `source` and `untested_reason` flags; MAJIQ segments in the MDACC list are mapped to IRFinder-S introns;
  - G×C at lead pairs (`xqtl_gxc.py`, port of atlas `_fit_gxc`), with maternal-contamination refit and per-pair primary cell type;
  - GxE (`xqtl_gxe.py`, `xqtl_gxe_core.py`, port of atlas `07_gxe_mapping`): exposure preparation, tier-1 adaptive-permutation scan of every phenotype, tier-2 nominal scan of all cis variants of the test-set phenotypes (`gxe.tier2: pairs` as a fallback), SNP × covariate sensitivity refits (atlas 55), and a status row per component;
  - `gxc` settings `fdr`, `min_n`, `maternal_threshold`; `gxe` settings `sensitivity`, `tier1_permutation_blocks`, `dosages`, `tier2`, `seed`, `fdr`. `gxe.transmitted` was removed (atlas 56 is not ported). The validator checks enabled sensitivity covariates.
- Export files G (G×C) and H (GxE); `I_strata` gains `gxc` and `gxe` run-state columns; `--mode export` rebuilds G and H.
- `--mode imputation_prep` (`workflows/imputation_prep.nf`, `modules/local/imputation_prep.nf`, `xqtl_impprep.sh`, `xqtl_impprep.py`; atlas `01_genotype_imputation`): VCF/BCF, PLINK 1 or PLINK 2 genotypes from arrays, WES or WGS to per-chromosome GRCh38 VCFs for the TOPMed Imputation Server. Chromosome names harmonized; GRCh37 input lifted with Picard LiftoverVcf; biallelic SNPs filtered on missingness, MAF and HWE (`imputation_prep` settings); REF/ALT aligned to GRCh38 without strand flipping; a reference-match rate below `imputation_prep.ref_match_min` stops the run ("Wrong build?") and publishes no VCF. Duplicate sample IDs stop the run; sample missingness and a chrX sex check against the metadata are reported. Outputs `imputation_prep_report.tsv`, `imputation_prep_samples.tsv`, `liftover_summary.tsv` and `TOPMED_CHECKLIST.md`. The upload stays manual. New parameter `liftover_chain`; the validator checks the raw files, FASTA and chain.
- Test fixture `tests/fixtures/imputation_prep/` (1000 Genomes phase 3 GRCh37 chr22 slice of 40 site-fixture samples, chr22 hg19→hg38 chains, GRCh38 PLINK set with REF = major allele; rebuilt by `make_fixture.sh`); `tests/scripts/test_imputation_prep.sh` runs the GRCh37, GRCh38 PLINK, wrong-build and duplicate-ID cases, also in CI.
- Test fixture: planted G × Cytotrophoblast and G × GA interactions, a `maternal_age` column, and MDACC test signals that exercise every test-set flag; `tests/scripts/check_site_interactions.py` checks their recovery in `test_site.sh`.
- `--mode central` (MDACC; `workflows/central.nf`, `modules/local/central.nf`, `bin/xqtl_central.py`): meta-analysis of the exports listed in `central.studies`.
  - Ingest: MANIFEST sha256 check, export-format version check (same major pipeline version and restricted-file version), file schemas, one contribution per study × mapped ancestry, allele harmonization to one orientation per site (swap, strand flip; palindromic sites dropped only when orientations differ and MAF > `central.palindromic_maf_max`), frequency-outlier flags.
  - Nominal fixed- and DerSimonian-Laird random-effects IVW (Q, I², τ²) of file A per modality × chromosome, within each ancestry and across all contributions; novel splicing clusters are excluded.
  - Feature level: ACAT of per-contribution permutation p-values (file B), Storey q-values (BH fallback), meta lead variant.
  - Multi-ancestry fine-mapping with SuShiE in summary-statistics mode (`bin/xqtl_finemap.py`, own `sushie` environment in the image) on within-ancestry meta z-scores and N-weighted in-sample LD (file C): PIPs, credible sets, cross-ancestry effect correlations and prediction weights.
  - Colocalization with GWAS summary statistics listed in `central.gwas_catalog` (`bin/xqtl_coloc.R`): SuSiE-coloc (coloc.abf fallback when a trait has no credible set) and ColocBoost, using the ancestry-matched xQTL LD for both traits.
  - G×C pair × cell-type IVW with central primary cell type; GxE tier-1 synthesis with scope labels (atlas 54) and tier-2 variant-level IVW after rescaling interaction effects to the pooled exposure SD.
  - TWAS input contract: OTTERS / mvOTTERS stage-I files per ancestry and FUSION `.wgt.RDat` weights from SuShiE (`bin/xqtl_fusion_wgt.R`).
  - `central_summary.tsv` and PIP enrichment (Fisher and PIP-weighted) for TSS-distance classes and annotation BEDs in `central.annotations`.
  - New `central` block in the cohort schema; the validator checks the studies table, export folders, GWAS catalog and annotation tables. Unit tests compare the IVW with a hand calculation and cover ACAT, harmonization, exposure rescaling and validation.
- Central acceptance test `tests/scripts/test_central.sh` (also in CI, skipped while the development image lacks SuShiE): `--mode site` for two cohorts, a synthetic GWAS (`tests/fixtures/central/make_gwas.py`), then `--mode central`. Second test cohort `tests/fixtures/site2/` (`make_central_fixture.sh`; 120 other 1000 Genomes samples, 60 EUR and 60 EAS, with the planted effects of the site fixture and other seeds, exported as `tsv.gz`); `make_site_fixture.py` gains `--second`. `check_central.py` checks the planted signals at every central stage, the colocalization posteriors, the tier-1 scope labels, the TWAS inputs and the summary; `check_central_ivw.R` recomputes the IVW in R from the raw file A rows and loads the FUSION weights. `make_gwas.py --also-export-dir` keeps the distinct (H3) GWAS variant in the variant set shared by all contributions; `check_site_truth.py --feature-q-warn` reports per-stratum feature q-values without requiring them.

### Fixed
- MDACC restricted-file builder: frozen cluster IDs now carry the strand (`clu_<n>_<strand>`), and atlas splicing signal IDs in the underscore form (`GENE__chr_start_end_clu_N_strand`) are converted to the frozen clusters.
- Image build: `containers/install_r_extra.R` now runs with the environment's `bin` on `PATH` so Rhcpp and colocboost find the conda C++ compiler; it also stops early if the bioconda `GenomeInfoDbData` post-link download failed.
- CI: base image pinned by digest directly in the `FROM` lines (hadolint DL3006); `actions/checkout` and `actions/setup-java` v5.
- Export file C: the `file` column of `*_C_ld.windows.tsv` now names the exported LD files (`<cohort>_<ANC>_C_ld_...`); central ingest also matches older exports by window coordinates.
- Colocalization: when one trait had no credible set, the `coloc.abf` fallback wrote empty (NaN) posteriors because the `.abf` suffix of the summary names was not stripped; `PP.H0`–`PP.H4` are now filled.
- Leak scan: text columns were skipped under pandas >= 3 (new `str` dtype), so sample IDs in string columns of parquet/TSV exports could pass; the scanner now checks every non-numeric, non-boolean, non-datetime column.
- CI unit tests: dependencies pinned to the image versions (`tests/requirements-ci.txt`), plus a second run with the latest pandas/pyarrow/numpy.
- `image.yml`: references are downloaded and cached in their own steps (cached even when the build fails), the book is rendered in its own step, and a failed `docker build` posts its last log lines as a public annotation and uploads `build.log`.

### Notes
- The CUDA image variant uses PyTorch `cu132` wheels (no older CUDA build exists for the pinned PyTorch release); GPU nodes need a driver supporting CUDA 13.2.

### Changed
- `image.yml` builds again on pushes to `main` that change the image inputs and on version tags (paused during development of the analysis modes; manual runs always worked).
- `conf/test.config` uses the site fixture and re-applies `resourceLimits` (2 CPU, 6 GB), because `base.config` evaluates the limits before profiles set `params.max_*`.

