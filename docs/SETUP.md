# Setup

This page takes a site from nothing to a working installation: get the image, check it, export the launcher, and confirm that the cluster can run the pipeline. It takes about 30 minutes plus the image download.

> **Restricted file: request it before you start.** Site runs need `xqtl-restricted-<version>.tar.gz`, a small file holding the shared (frozen) splicing clusters and the MDACC signal list. It is not in the image or on any public server. MDACC sends it **by email** to sites with an approved data-sharing agreement. If you have not received it, request it from **Sean T. Bresnahan (stbresnahan@mdanderson.org)**, stating:
>
> 1. your cohort name,
> 2. the pipeline version you will run (for example `1.0.0`), and
> 3. the approved agreement (title or reference number).
>
> Keep the file inside your secure environment and set its path as `restricted:` in `cohort.yml`. The smoke test (step 5) reports whether it is present and intact.

## What you need

| Requirement | Details |
|---|---|
| Linux, x86-64 | Login/head node plus compute nodes (or a single large workstation) |
| Apptainer ≥ 1.0 or Singularity ≥ 3.8 | On every node that runs tasks. Check with `apptainer --version` or `singularity --version` |
| A batch scheduler (optional) | LSF, SLURM, PBS Pro/Torque or SGE. See [HPC configuration](HPC.md) |
| Shared filesystem | Visible from the head node and all compute nodes, for the work directory and results |
| Disk | The image file (its size is listed on the release page and the Zenodo record), plus work space for your cohort (see [HPC configuration](HPC.md#work-directory-temporary-space-and-cleanup)) |
| Internet | Only to download the image once. Nothing at run time |

Nothing else is installed on the host: no conda, Java, Nextflow, R or Python. All of them come from the image.

## 1. Get the image

The image is one file, `xqtl-federated_<version>.sif`, available from two places that hold the same bytes.

**From the GitHub container registry** (needs outbound HTTPS to `ghcr.io` from the machine you download on):

```bash
apptainer pull xqtl-federated_1.0.0.sif oras://ghcr.io/sbresnahan/xqtl-federated-sif:1.0.0
# Singularity: singularity pull xqtl-federated_1.0.0.sif oras://ghcr.io/sbresnahan/xqtl-federated-sif:1.0.0
```

**From Zenodo** (a plain HTTPS download; works through most proxies and can be done on a laptop and copied over): download `xqtl-federated_<version>.sif` and `xqtl-federated_<version>.sif.sha256` from the Zenodo record linked on the GitHub release page.

**Verify the download.** The checksum file is attached to the GitHub release and to the Zenodo record:

```bash
sha256sum -c xqtl-federated_1.0.0.sif.sha256        # must print: xqtl-federated_1.0.0.sif: OK
```

Put the `.sif` on the shared filesystem, for example next to where you will keep the launcher. One copy serves every run; compute nodes read it directly, so no image cache is needed.

A Docker image of the same build (`ghcr.io/sbresnahan/placenta-xqtl-federated:<version>`) exists for workstations and continuous integration. Most clusters do not allow Docker; use the `.sif`.

## 2. Check the image

```bash
apptainer run xqtl-federated_1.0.0.sif help          # commands available in the image
apptainer run xqtl-federated_1.0.0.sif versions      # software versions (TSV)
apptainer run xqtl-federated_1.0.0.sif refs-check    # verify the embedded reference data (sha256; ~1 min)
```

`refs-check` compares every embedded reference file with the checksums recorded for this release and ends with `OK <n> reference files verified`. A failure means the image file is damaged; download it again.

## 3. Export the launcher (`bootstrap`)

Nextflow's head job runs on the host, outside the container, so it can submit jobs to your scheduler. `bootstrap` copies what the head job needs out of the image into a directory you choose:

```bash
apptainer run xqtl-federated_1.0.0.sif bootstrap ./xqtl
```

| Written | Contents |
|---|---|
| `xqtl/run_xqtl.sh` | The launcher you run from now on |
| `xqtl/pipeline/` | Pipeline code (identical to the copy in the image) |
| `xqtl/jre/` | Java 17 runtime for the head job |
| `xqtl/nextflow/` | Nextflow launcher (self-contained build) and its home directory |
| `xqtl/docs/` | This documentation for offline reading: open `xqtl/docs/index.html` |
| `xqtl/VERSIONS.tsv`, `xqtl/XQTL_VERSION`, `xqtl/refs_MANIFEST.tsv` | Software versions, pipeline version, reference checksums |

The reference data stay inside the image; tasks read them at `/opt/xqtl/refs`.

Notes:

- The destination must be on a filesystem visible inside the container. Apptainer binds `$HOME` and the current directory by default; for other locations add `-B /path` (for example `apptainer run -B /scratch xqtl-federated_1.0.0.sif bootstrap /scratch/xqtl`).
- `bootstrap` records the path of the `.sif` in `run_xqtl.sh`. If you move the image later, set `XQTL_SIF=/new/path/xqtl-federated_1.0.0.sif` before running, or bootstrap again.
- Use **one directory per pipeline version**. `bootstrap` refuses to write into a directory that already holds an installation; `--force` replaces the exported files and leaves your own files, work directories and results untouched.
- If `bootstrap` says it cannot tell where the `.sif` is (this happens with some older Singularity versions), add `--sif /full/path/xqtl-federated_1.0.0.sif`.

## 4. What `run_xqtl.sh` does

`run_xqtl.sh` runs `nextflow run xqtl/pipeline/main.nf --sif <image> <your arguments>` with:

| Setting | Value | Why |
|---|---|---|
| `JAVA_HOME` | `xqtl/jre` | No Java needed on the host |
| `NXF_HOME` | `xqtl/nextflow/home` | Nextflow keeps its files inside the installation, not in `$HOME/.nextflow` |
| `NXF_OFFLINE` | `true` | Nextflow never tries to reach the internet (on seadragon each attempt cost about 130 s) |
| `NXF_OPTS` | `-Xms1g -Xmx4g` unless you set it | Head-job memory; see [HPC configuration](HPC.md) |
| `XQTL_SIF` | the `.sif` path | Every task runs in this image |

Any `NXF_VER`, `NXF_JAVA_HOME` or `JAVA_CMD` in your environment is ignored, so a Nextflow or Java installation already on the host cannot interfere.

Run it from the directory where you want the work directory and results, usually on scratch:

```bash
cd /scratch/$USER/xqtl_run
/path/to/xqtl/run_xqtl.sh -profile <site profile> -params-file cohort.yml --mode <mode>
```

Long runs should be launched as a batch job or inside `tmux`/`screen`; see [Running the head inside a batch job](HPC.md#running-the-head-inside-a-batch-job) and the templates in `pipeline/assets/head_job/`.

## 5. Check that your cluster works (smoke test)

```bash
/path/to/xqtl/run_xqtl.sh -profile seadragon,smoke --restricted /secure/path/xqtl-restricted-1.tar.gz
# other sites:  -profile slurm,apptainer,smoke -c conf/sites/<yoursite>.config ...
```

The smoke test takes about 10 minutes and runs 12 checks (container start-up, tools, references, restricted file, job fan-out and arrays, memory enforcement, filesystem binds, a small tensorQTL run). It writes `xqtl_smoke/smoke/smoke_summary.tsv` ending in `VERDICT GO` or `VERDICT NO-GO`. [HPC configuration](HPC.md#validating-a-new-site) explains each check and what to do about a FAIL or WARN. Send MDACC the summary file if you would like a second opinion; it holds no data.

## 6. Build the alignment indices (once per site)

STAR, Salmon and IRFinder-S indices are too large to ship in the image, so each site builds them once from the references inside it:

```bash
/path/to/xqtl/run_xqtl.sh -profile <site profile> -params-file cohort.yml --mode prepare_references
```

The indices are written to `site_refs:` in `cohort.yml` and reused by every later run. This mode becomes available with the phenotyping release (see the [changelog](../CHANGELOG.md)).

## Upgrading to a new pipeline version

1. Download and verify the new `.sif` (step 1).
2. Bootstrap it into a new directory, for example `./xqtl-1.1.0`.
3. Check the [changelog](../CHANGELOG.md): a new major version changes the export format or the restricted-file version, and MDACC will send a new restricted file.
4. Re-run the smoke test. Site indices from `prepare_references` are reused unless the changelog says otherwise.

Keep the old directory until the new version has run; results record the version that produced them.

## Offline documentation

Every image carries this documentation. After `bootstrap`, open `xqtl/docs/index.html` in a browser (copy the folder to a laptop if the cluster has no browser). `apptainer run <sif> docs` prints its location inside the image. The online copy is at <https://sbresnahan.github.io/placenta-xqtl-federated/>.
