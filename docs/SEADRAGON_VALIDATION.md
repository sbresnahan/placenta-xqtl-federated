# Seadragon validation (Phase 0)

Before the pipeline was written, a stand-alone smoke kit checked that every scheduler and container behaviour the pipeline relies on works on MD Anderson's seadragon cluster. It was run on 2026-10-29 in both supported head modes. Both runs gave **GO**. The kit's checks are now built into the pipeline as `--mode smoke` (see [Validating a new site](HPC.md#validating-a-new-site)); test identifiers below are those of the Phase 0 kit, with the pipeline equivalent in brackets.

## Environment

| Item | Value |
|---|---|
| Scheduler | IBM Spectrum LSF |
| Container engine | Apptainer 1.5.1, invoked as `singularity` |
| Nextflow | 26.10.4 |
| Head, run 1 | Login node; tasks as LSF jobs (job array 79979021) |
| Head, run 2 | Inside LSF job 79979246 (queue `medium`, 1 core, 4 GB, host cdragon233) |
| Images | R 4.3.1 (+ Bioconductor `qvalue`), tensorQTL 1.0.10 / torch 2.15.0 (CPU), plink2 |

## Results (run 2, head inside an LSF job)

| Phase 0 test [pipeline] | Status | Observation |
|---|---|---|
| T1/T2 container R [S01] | PASS | R 4.3.1 ran inside the container on compute node cdragon208 |
| T3 plink2 [S05] | PASS | `--freq`, `--hardy` and `--r2-unphased` on a dummy dataset |
| T4 resources [S06] | PASS | 4 CPUs and 8 GB granted (`rusage[mem=8192MB]`, memlimit 8 G); `nproc=4`: CPUs confined to the request |
| T5 fan-out [S07] | PASS | 5 independent jobs on 4 different nodes |
| T5b job array [S08] | PASS | One LSF array job (79979266[1–5]) with `LSB_JOBINDEX` 1–5 |
| T6 bind read [S09] | PASS | A file on `/rsrch5` was readable inside the container |
| T7 memory enforcement and retry [S10] | PASS | Attempt 1 (2 GB request, 3 GB used) killed with exit 130; attempt 2 with 4 GB succeeded |
| T9 tensorQTL CPU [S11] | PASS | tensorQTL 1.0.10, 4 threads; planted feature `pval_beta` = 1.97 × 10⁻¹⁴ |
| T9b tensorQTL GPU [S11] | SKIP | Not requested |
| T9q Storey q-value [S12] | PASS | `compute_qvalues.R` bridge, qvalue 2.34.0 on R 4.3.1; planted feature q = 7.87 × 10⁻¹³ |
| T8 resume | PASS | 18 of 18 tasks reported as cached on `-resume` |
| R1 run exit | PASS | Nextflow exit status 0 |

Run 1 (head on a login node) gave the same results.

## Findings carried into the pipeline

- **Queues.** The site's submission filter enforces walltime ranges (`short` ≤ 3 h; `medium` more than 3 h and under 24 h; `long` 24–120 h). The seadragon profile therefore picks the queue from each task's walltime (`queue_rules`).
- **Memory.** `LSF_UNIT_FOR_LIMITS=GB` and `LSB_MEMLIMIT_ENFORCE=y`: memory is requested per job and enforced, and over-limit jobs exit with status 130. Exit 130 is in the pipeline's retry list.
- **CPU confinement.** `LSB_RESOURCE_ENFORCE="cpu gpu"`.
- **No internet on compute nodes.** Nextflow spent about 130 s per launch waiting for a connection to `www.nextflow.io`. `run_xqtl.sh` sets `NXF_OFFLINE=true`, and the image contains everything needed at run time.
- **Registries.** From login nodes `ghcr.io` and `quay.io` worked at about 0.6 MB/s; Docker Hub, `gcr.io` and Amazon ECR were blocked or unreliable. The pipeline image is published on `ghcr.io` and Zenodo.
- **No image builds on seadragon.** `apptainer build --fakeroot` with apt failed (exit 100). Images are built by GitHub Actions only.
- **`ps` inside the image.** Nextflow's task wrapper needs `ps`; the image includes `procps`.
- **q-values without rpy2.** The atlas file bridge (`compute_qvalues.R`) is used for Storey q-values.
