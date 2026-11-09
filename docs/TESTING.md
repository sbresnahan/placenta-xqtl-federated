# Testing and validation

This page is for maintainers and developers. It describes the automated tests and their test data, the checks run on MD Anderson's seadragon cluster, and how the image is built and released. Sites do not need any of it: they validate their own cluster with the smoke test ([Validating a new site](HPC.md#validating-a-new-site)).

All test data are either public or synthetic. The test cohorts use real 1000 Genomes genotypes (chr22) with **synthetic phenotypes carrying planted effects**: they check that the software recovers known signals, not anything about placental biology. The only real placental RNA-seq in the tests is 8 public runs of BioProject PRJNA820329, cut to chr22.

## Development environment

```bash
scripts/dev_env.sh /path/to/env        # conda/micromamba environment matching the image
export PATH=/path/to/env/bin:/path/to/env/align/bin:/path/to/env/xqtl/bin:/path/to/env/irfinder/bin:$PWD/bin:$PATH
python -m pytest tests/unit            # unit tests
scripts/check_configs.sh               # every profile resolves
nextflow run main.nf -profile local,dev,smoke --refs_dir tests/fixtures/refs_mini \
    --restricted tests/fixtures/xqtl-restricted-test.tar.gz
Rscript book/render_book.R             # documentation site into book/_book
tests/scripts/test_bootstrap.sh        # bootstrap + run_xqtl.sh against a mock image tree
```

`-profile dev` (`conf/containers/none.config`) runs every tool from `PATH` without a container. It is for development only, not for production runs. `-profile test` (`conf/test.config`) points every input at the test data below and caps resources at 2 CPUs and 6 GB; combine it with an engine, e.g. `-profile test,docker` or `-profile test,dev`.

## Unit tests

`tests/unit/` (pytest) covers input validation, the restricted file, reference download and checks, the MDACC packaging scripts, splicing clusters, mapping and export helpers (including the leak scan), G×C and GxE statistics, imputation preparation and the central analysis. The central tests compare the IVW meta-analysis with a hand calculation and cover ACAT, allele harmonization and exposure rescaling. `test_docs.py` checks that every `cohort.yml` setting in the schema is documented in [Configuration](CONFIGURATION.md).

## End-to-end tests

Each script runs the pipeline on test data and then checks the outputs. Usage: `tests/scripts/<test>.sh [WORKDIR] [ENGINE]`, with `ENGINE` = `dev` (tools on `PATH`) or `docker`.

### Site run: `test_site.sh`

Test cohort `tests/fixtures/site/` (built by `make_site_fixture.sh` and `make_site_fixture.py`): 130 1000 Genomes samples on chr22 with synthetic phenotypes and planted cis effects (`TRUTH.tsv`). The test runs `--mode site` with `-profile test` and checks that:

- the planted cis effects are recovered in EUR and EAS (q ≤ 0.05 in file B; correct sign and p ≤ 10⁻³ at the planted variant in file A; `check_site_truth.py`);
- AFR (n = 10, below `ancestry.min_n`) is listed in `I_strata` but not mapped or exported;
- a sample ID planted in a site output makes `--mode export` fail and leaves no export behind, and `--mode export` rebuilds the same export once it is removed.

**Interactions.** The site cohort also carries planted interactions:

- a G × Cytotrophoblast interaction on `ENSGFIX0000011`;
- a G × GA interaction on `ENSGFIX0000021` (planted variant chr22:20290948:A:G);
- a G × GA splicing plant on clu_31.

The restricted test file (`tests/fixtures/xqtl-restricted-test.tar.gz`, built by `make_restricted_test.sh`) adds MDACC pairs that exercise every test-set flag: a phenotype absent from the site, a variant absent from the site, a different MDACC lead for a site eGene, shared pairs, and an MDACC-only pair. With GA and the maternal-age sensitivity enabled, `check_site_interactions.py` checks, in EUR and EAS:

- the primary cell type of `ENSGFIX0000011`;
- the tier-1 q-value and sign for `ENSGFIX0000021`;
- the tier-2 top variant for `ENSGFIX0000021`;
- the test-set sources and untested reasons;
- H_status;
- the sensitivity refit.

It asserts only the planted pairs: with 60 samples per stratum, another pair can reach q ≤ 0.05 by chance. GxE recovery is checked on the expression plant only; normalization of splicing ratios attenuates the splicing plant too much to detect at n = 60, so that plant is checked only for its main effect and test-set flags.

### Central analysis: `test_central.sh`

`test_central.sh` runs `--mode site` for two independent test cohorts, then `--mode central`.

- **Cohort 1** is the site cohort above (parquet export).
- **Cohort 2** (`tests/fixtures/site2/`, built by `make_central_fixture.sh`) has 120 other 1000 Genomes samples (60 EUR, 60 EAS), the same planted effects with other random seeds, and is exported as `tsv.gz`. The cohorts do not share samples, and each stratum has at least `ancestry.min_n` = 50 samples.
- **GWAS.** `tests/fixtures/central/make_gwas.py` writes a synthetic GWAS from one export's LD window: one causal variant shared with a planted cis effect (colocalizing) and one nearby independent causal variant (a distinct signal, H3).
- The test adds a BED of the planted variants (±1 kb) as an enrichment annotation and gives study 2 as a path relative to `studies.tsv`.

`check_central.py` checks ingest, the planted signals at meta, feature, fine-mapping, colocalization, G×C and GxE level, the tier-1 scope labels, the TWAS input files and the summary. `check_central_ivw.R` recomputes the fixed- and random-effects IVW (Q, I², τ²) in R from the raw file A rows of both exports and loads every FUSION weight file.

### Phenotyping: `test_phenotyping.sh`

`test_phenotyping.sh` runs `--mode prepare_references` and `--mode phenotyping_only` (`-profile test_phenotyping`) and checks the phenotype BEDs (`check_phenotypes.py`).

- **Reads.** 8 public placenta RNA-seq runs of BioProject PRJNA820329 (Seoul National University Hospital), as read pairs with both mates on chr22 (`tests/fixtures/fastq/`, `tests/fixtures/phenotyping/`). `extract_snuh_chr22.sh` cuts them from STAR alignments and `make_fastq_fixture.sh` writes the fixture. Base qualities are set to a constant (Phred 30), as in SRA Lite, so base-quality filters act differently from real data.
- **Annotation.** GENCODE v45 chr22 (`gencode.v45.annotation.gtf.gz`, sha256-pinned in `make_refs_fixture.sh`), written in the HPLRv2 PANTRY GTF form and processed with Pantry's txrevise scripts, is a public stand-in for HPLRv2 (`tests/fixtures/refs_chr22/`). Real runs use the HPLRv2 annotation in the image.
- **Editing sites and cell types.** `refs_chr22/editing/` and `music/` are chr22 cuts of the public reference record (`extract_refs_chr22.sh`): 413,357 REDIportal sites on chr22, and a mini single-cell reference of 254 chr22 genes and 3,857 cells (at most 30 cells per cell type and donor; 29 cell types, 5 donors). With so few genes, the test's cell proportions only check that the step runs; they are not realistic estimates.
- **Splicing.** The test's restricted file freezes every other leafCutter cluster of these samples (`make_phenotyping_restricted.sh`), so both frozen and novel clusters are exercised.
- The chr22 genome is cut from the image's GRCh38 FASTA. On it, IRFinder-S writes an empty `intergenic.ROI.bed`, because it keeps intergenic regions only on chromosome names of up to two characters. This does not affect IR ratios.

### Imputation preparation: `test_imputation_prep.sh`

Test data `tests/fixtures/imputation_prep/` (rebuilt by `make_fixture.sh`): a 1000 Genomes phase 3 GRCh37 chr22 slice of 40 site-cohort samples, chr22 hg19 → hg38 chains, and a GRCh38 PLINK set with REF = major allele. The test checks:

1. GRCh37 VCF → liftover → `chr22.vcf.gz` at the GRCh38 positions and alleles of the same samples' GRCh38 genotypes, plus the report, liftover summary and TOPMed checklist;
2. GRCh38 PLINK input → REF/ALT swapped to GRCh38 with genotypes intact;
3. the GRCh37 VCF declared as GRCh38 stops with "Wrong build?" and publishes no VCF;
4. duplicate sample IDs stop the run.

### Smoke test and bootstrap

`tests/fixtures/refs_mini/` is a small stand-in reference set (labelled "TEST FIXTURE: not real reference data") for running `--mode smoke` outside the image. `test_bootstrap.sh` runs `containers/bootstrap.sh` and the generated `run_xqtl.sh` against a directory mock of the image, with the pinned Java runtime and Nextflow downloaded and verified as in the image build.

## Continuous integration

`.github/workflows/ci.yml` runs on every push and pull request:

- **unit tests**, with dependencies at the image versions (`tests/requirements-ci.txt`) and again with the latest pandas, pyarrow and numpy;
- **configs and linters:** every profile resolves (`scripts/check_configs.sh`), `nextflow lint`, executable bits, shellcheck, hadolint and actionlint;
- **image checks** in the development image `ghcr.io/sbresnahan/placenta-xqtl-federated:dev`: `versions`, `refs-check`, the smoke test, and the site, central, phenotyping and imputation-preparation tests above. These steps are skipped, with a notice, while no development image is published.

`.github/workflows/image.yml` builds the development image on pushes to `main` that change the image inputs, and the release images on version tags (see below). `.github/workflows/docs.yml` renders the documentation site.

## Cluster checks on seadragon

The scheduler and container behaviour the pipeline relies on was checked on seadragon on 2026-10-29, with a stand-alone set of checks equivalent to `--mode smoke`, in both supported head modes. Both runs gave **GO**. Smoke-test identifiers are in brackets.

| Item | Value |
|---|---|
| Scheduler | IBM Spectrum LSF |
| Container engine | Apptainer 1.5.1, invoked as `singularity` |
| Nextflow | 26.10.4 |
| Head, run 1 | Login node; tasks as LSF jobs (job array 79979021) |
| Head, run 2 | Inside LSF job 79979246 (queue `medium`, 1 core, 4 GB, host cdragon233) |
| Software | R 4.3.1 (+ Bioconductor `qvalue`), tensorQTL 1.0.10 / torch 2.15.0 (CPU), plink2 |

Results of run 2 (head inside an LSF job); run 1 gave the same results:

| Check [smoke test] | Status | Observation |
|---|---|---|
| Container start, R [S01] | PASS | R 4.3.1 ran inside the container on compute node cdragon208 |
| plink2 [S05] | PASS | `--freq`, `--hardy` and `--r2-unphased` on a dummy dataset |
| Resources [S06] | PASS | 4 CPUs and 8 GB granted (`rusage[mem=8192MB]`, memlimit 8 G); `nproc=4`: CPUs confined to the request |
| Fan-out [S07] | PASS | 5 independent jobs on 4 different nodes |
| Job array [S08] | PASS | One LSF array job (79979266[1–5]) with `LSB_JOBINDEX` 1–5 |
| Bind read [S09] | PASS | A file on `/rsrch5` was readable inside the container |
| Memory enforcement and retry [S10] | PASS | Attempt 1 (2 GB request, 3 GB used) killed with exit 130; attempt 2 with 4 GB succeeded |
| tensorQTL CPU [S11] | PASS | tensorQTL 1.0.10, 4 threads; planted feature `pval_beta` = 1.97 × 10⁻¹⁴ |
| tensorQTL GPU [S11] | SKIP | Not requested |
| Storey q-value [S12] | PASS | `compute_qvalues.R` bridge, qvalue 2.34.0 on R 4.3.1; planted feature q = 7.87 × 10⁻¹³ |
| Resume | PASS | 18 of 18 tasks reported as cached on `-resume` |
| Run exit | PASS | Nextflow exit status 0 |

How these results shape the `seadragon` profile and the image:

- **Queues.** The site's submission filter enforces walltime ranges (`short` ≤ 3 h; `medium` more than 3 h and under 24 h; `long` 24–120 h), so the profile picks the queue from each task's walltime (`queue_rules`).
- **Memory.** `LSF_UNIT_FOR_LIMITS=GB` and `LSB_MEMLIMIT_ENFORCE=y`: memory is requested per job and enforced, and over-limit jobs exit with status 130, which is in the pipeline's retry list.
- **CPU confinement.** `LSB_RESOURCE_ENFORCE="cpu gpu"`.
- **No internet on compute nodes.** Without `NXF_OFFLINE`, Nextflow waits about 130 s per launch for `www.nextflow.io`; `run_xqtl.sh` sets `NXF_OFFLINE=true`, and the image contains everything needed at run time.
- **Registries.** From login nodes `ghcr.io` and `quay.io` worked at about 0.6 MB/s; Docker Hub, `gcr.io` and Amazon ECR were blocked or unreliable. The image is published on `ghcr.io` and Zenodo.
- **No image builds on seadragon.** `apptainer build --fakeroot` with apt fails (exit 100); images are built by GitHub Actions.
- **`ps` inside the image.** Nextflow's task wrapper needs `ps`; the image includes `procps`.
- **q-values without rpy2.** Storey q-values are computed in R through `compute_qvalues.R`.

## Building and releasing

### Building an image

`containers/build_image.sh` downloads and verifies the references, renders this documentation if needed, and runs `docker buildx build` with the right build contexts and arguments:

```bash
containers/build_image.sh                       # development build (unpinned rows allowed)
containers/build_image.sh --variant cuda        # GPU variant
containers/build_image.sh --release             # every row pinned; references must equal assets/refs_manifest.tsv
```

### Releasing a version

1. Pin the public reference record ([Reference data](REFERENCES.md#building-the-public-reference-record-once-per-reference-release)) and push to `main`. The `image` workflow builds `ghcr.io/sbresnahan/placenta-xqtl-federated:dev` and uploads `refs_MANIFEST.tsv` as an artifact.
2. Commit that file as `assets/refs_manifest.tsv` whenever a reference changes.
3. Set `manifest.version` in `nextflow.config`, update the changelog, and push the tag `v<version>`. The workflow then builds the CPU and `-cuda` images in strict mode, pushes them, builds the `.sif`, pushes it to `ghcr.io/sbresnahan/xqtl-federated-sif:<version>` (pulled with `apptainer pull oras://...`), archives it on Zenodo, and attaches the `.sif` checksum, reference manifest and version list to the GitHub release.

Release checklist (in order):

1. All CI checks pass on `main`.
2. The public reference record is pinned in `assets/refs_sources.tsv` (no unpinned rows) and `assets/refs_manifest.tsv` matches the `refs_MANIFEST.tsv` artifact of the last `:dev` build.
3. The restricted file for this version is built, its version equals `expected_restricted_version` in `nextflow.config`, and its sha256 is recorded with the MDACC records.
4. The `:dev` image passes the smoke test on seadragon (`-profile seadragon,smoke`, `VERDICT GO`) and a full `--mode site` run on a real cohort subset.
5. `manifest.version` in `nextflow.config` and `params.image` carry the new version; the changelog section is dated; the export-format or restricted-file version change, if any, is a major version change.
6. The documentation book renders with no broken links (`Rscript book/render_book.R`).
7. Push the tag `v<version>`; after the workflow finishes, check the GitHub release assets, the GHCR packages, the Zenodo deposition (publish it if `ZENODO_PUBLISH` is not set) and `apptainer pull oras://ghcr.io/sbresnahan/xqtl-federated-sif:<version>` from a site-like machine.
8. Email the restricted file to approved sites and announce the release with the SIF checksum.

One-time repository settings:

- **Workflow permissions**: Settings → Actions → General → Workflow permissions → Read and write.
- **Secret `ZENODO_TOKEN`**: a Zenodo personal access token with `deposit:write` and `deposit:actions`. Optional repository variables: `ZENODO_CONCEPT_ID` (an earlier deposition of the image record, so every release becomes a new version under one concept DOI) and `ZENODO_PUBLISH=true` (publish automatically; otherwise the deposition is left as a draft to review on zenodo.org).
- **Package visibility**: after the first push, set both GHCR packages (`placenta-xqtl-federated` and `xqtl-federated-sif`) to public under the repository's Packages settings, so sites can pull without logging in.
- **Documentation site**: Settings → Pages → Deploy from a branch → `bookdown` / (root). The `docs` workflow renders the book and pushes it to that branch.
