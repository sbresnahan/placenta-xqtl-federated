# HPC configuration guide

This guide explains how the pipeline requests compute resources and how to adapt it to your cluster. Most sites only need to choose a profile, copy one site template and set five or six values. Validate every new configuration with the smoke test (see [Validating a new site](#validating-a-new-site)) before running real data.

The pipeline runs in two parts:

- **The head** is a small Java process (Nextflow) that runs on a login node, a workflow node or inside a batch job. It submits every task to the scheduler and tracks it. It needs about 1 CPU and 2–4 GB of memory for the whole run, which can take several days.
- **Tasks** are batch jobs. Each task runs inside the same pipeline image (a `.sif` file, run with Apptainer or Singularity) on a compute node.

Neither part needs internet access. The image contains all software and public reference data, and `run_xqtl.sh` (written by the image's `bootstrap` command; see `SETUP.md`) starts Nextflow in offline mode.

## Configuration layers

Resource settings are layered so that the pipeline code never contains anything site-specific:

| Layer | File | What it sets | Who edits it |
|---|---|---|---|
| 1. Process labels | `modules/`, `workflows/` | Each process declares a resource class (`label 'process_medium'`), never a queue, account or path | Pipeline developers |
| 2. Resource classes | `conf/base.config` | CPUs, memory and walltime for each class; retry policy; queue selection rule; caps | Rarely changed |
| 3. Executor profile | `conf/executors/<scheduler>.config` | Scheduler type, job submission rate, polling | Rarely changed |
| 4. Container profile | `conf/containers/<engine>.config` | Apptainer, Singularity or Docker; bind paths | Rarely changed |
| 5. Site profile | `conf/sites/<site>.config` | Queues, account codes, GPU syntax, bind paths, node limits | **Each site** |
| 6. Command line | `--max_memory 64.GB` etc. | Overrides any parameter for one run | Each run |

Later layers override earlier ones. Profiles are combined with commas, for example `-profile slurm,apptainer`, and a site file is added with `-c conf/sites/<site>.config`. Seadragon has a built-in profile that bundles layers 3–5: `-profile seadragon`.

To see the configuration that a combination resolves to, without running anything:

```bash
nextflow config -profile slurm,apptainer -c conf/sites/mysite.config .
```

### Resource classes

Defined in `conf/base.config`. Memory and time are multiplied by the attempt number, so a task that is killed for exceeding its memory or walltime is resubmitted with twice (then three times) as much.

| Label | CPUs | Memory (attempt 1) | Walltime (attempt 1) | Used for |
|---|---|---|---|---|
| `process_single` | 1 | 4 GB | 2 h | Small bookkeeping, validation, plotting |
| `process_low` | 2 | 8 GB | 4 h | Per-sample light steps (featureCounts, RegTools) |
| `process_medium` | 8 | 32 GB | 8 h | STAR/Salmon per sample, normalization, HCP fits |
| `process_high` | 16 | 64 GB | 16 h | Index builds, cohort-level aggregation |
| `process_long` | (class) | (class) | 48 h | Added to long tasks (permutation scans) |
| `process_high_memory` | (class) | 128 GB | (class) | STAR genome index build |
| `process_gpu_capable` | 8 (CPU) / 4 (GPU) | 32 GB | 12 h | tensorQTL scans; GPU only with `--use_gpu` |
| `process_gpu` | 4 | 32 GB | 4 h | Steps that require a GPU (none in round 1) |
| `process_local` | 1 | 2 GB | 1 h | Trivial steps |
| `process_array` | (class) | (class) | (class) | Added to many-small-task steps: submitted as job arrays |

To change a class for your site, override it in your site file:

```groovy
process {
    withLabel: 'process_medium' {
        cpus   = 12
        memory = { 48.GB * task.attempt }
    }
}
```

To change a single process, use its name: `withName: 'STAR_ALIGN' { cpus = 16 }`. Process names are listed in the execution report (`<outdir>/pipeline_info/report_*.html`).

### Caps (`--max_cpus`, `--max_memory`, `--max_time`)

Every request is capped at `--max_cpus`, `--max_memory` and `--max_time` (Nextflow `process.resourceLimits`). Set these to the largest node and longest queue you may use; a request above the cap is lowered to the cap instead of being rejected by the scheduler. Defaults are 16 CPUs, 128 GB and 120 h.

Check node sizes with `bhosts -l` / `lshosts` (LSF), `sinfo -o "%n %c %m"` (SLURM), `pbsnodes -a` (PBS) or `qhost` (SGE).

## Executors

| Scheduler | Profile | Notes |
|---|---|---|
| IBM Spectrum LSF | `lsf` | Memory is per job (`perJobMemLimit`); see [Memory](#memory) |
| SLURM | `slurm` | `--mem` per job by default; `--mem-per-cpu` option in the template |
| PBS Pro / OpenPBS | `pbspro` | `select=1:ncpus=N:mem=M` syntax |
| Torque / PBS | `pbs` | `nodes=1:ppn=N,mem=M` syntax |
| SGE / UGE / Son of Grid Engine | `sge` | Uses parallel environment `smp`; change it if your site uses another name |
| No scheduler (one large node) | `local` | All tasks run on the head's node; capped by `--max_cpus` / `--max_memory` |

Site templates with every value marked `EDIT` are in `conf/sites/`:

- `TEMPLATE_slurm.config`
- `TEMPLATE_lsf.config`
- `TEMPLATE_pbspro.config` (also works for Torque with `-profile pbs`)
- `TEMPLATE_sge.config`

Copy one to `conf/sites/<yoursite>.config`, edit it and run with:

```bash
./run_xqtl.sh -profile slurm,apptainer -c conf/sites/<yoursite>.config -params-file cohort.yml --mode smoke
```

### Single node (`-profile local`)

On a workstation or one large allocated node, `-profile local,apptainer` runs every task on that machine, never using more than `--max_cpus` and `--max_memory` at once. A full cohort needs many CPU-days (see [Tuning](#tuning-table)); this mode suits small cohorts, the smoke test and re-running the export.

## Queues and partitions

Many schedulers reject or misroute jobs whose walltime does not fit the queue. The pipeline chooses a queue for each task from its walltime with `--queue_rules`:

```text
--queue_rules 'short:3h,medium:23h59m,long:120h'
```

Rules are read left to right; the first queue whose limit is at least the task's walltime is used, and the last queue takes anything longer. Because walltime grows with each retry, a retried task can move to a longer queue automatically. List queues in increasing walltime order.

- `--queue NAME` sends every task to one queue and ignores the rules.
- With neither set, no queue is requested and the scheduler's default is used.
- `--gpu_queue` is used only for GPU tasks with `--use_gpu`.

Find your limits with `bqueues -l` (LSF), `scontrol show partition` (SLURM), `qstat -Qf` (PBS) or `qconf -sq <queue>` (SGE).

## Memory

The scheduler must receive the same memory number that the task will use, and in the units it expects. Getting this wrong is the most common cause of jobs being killed or never starting.

### LSF

- `perJobMemLimit = true` (the default here, `--lsf_per_job_mem true`) tells Nextflow that memory limits apply to the whole job. Without it, Nextflow divides the request by the number of CPUs, and multi-core tasks are killed far below their request on clusters configured per job.
- Nextflow reads `LSF_UNIT_FOR_LIMITS` from `$LSF_ENVDIR/lsf.conf` to choose the unit for `-M` and `rusage[mem=...]`. The head must therefore run in a shell where the LSF environment is loaded (`LSF_ENVDIR` set). Check one submitted script (`.command.run` in any task directory) and confirm that the `#BSUB -M` value is in the unit your site expects.
- If `LSB_MEMLIMIT_ENFORCE=y` (as on seadragon), a task that exceeds its request is killed with exit status 130 and retried with more memory.

To check how your site is configured:

```bash
grep -E 'LSF_UNIT_FOR_LIMITS|LSB_MEMLIMIT_ENFORCE|LSB_RESOURCE_ENFORCE|LSB_JOB_MEMLIMIT' $LSF_ENVDIR/lsf.conf
```

### SLURM

- By default the pipeline requests memory per job (`--mem`). This is correct on almost all SLURM sites.
- Some sites require `--mem-per-cpu`. The SLURM template contains a commented `clusterOptions` line that converts the request; enable it if jobs are rejected with a memory-specification error.
- Sites that set `DefMemPerCPU` without enforcing `--mem` can still kill jobs through cgroups; the smoke test (S10) shows whether memory is enforced.

### PBS and SGE

- PBS Pro: memory goes into the `select` statement (`mem=`), per job.
- Torque: `-l mem=` per job; some Torque sites need `pmem` (per process) instead. Add it with `--cluster_options`.
- SGE: memory is requested with `-l h_vmem` (per slot on many sites). If your site expects per-slot values, divide the request in your site file:

  ```groovy
  process.clusterOptions = { "-l h_vmem=${(task.memory.toMega() / task.cpus) as int}M" }
  ```

### Retries after memory or time kills

`conf/base.config` retries a task up to 3 times when it exits with a status that schedulers use for resource kills: 104, 130, 134, 137, 139, 140, 143, 247 or 271. Other failures stop the run after running tasks finish (`errorStrategy 'finish'`), so you can fix the problem and continue with `-resume`.

If a task fails with a resource error on all attempts, raise its class or the cap in your site file rather than editing the pipeline.

## CPUs and threads

Each tool receives the number of CPUs from `task.cpus` explicitly (for example `STAR --runThreadN ${task.cpus}`, `salmon -p`, `plink2 --threads`, `torch.set_num_threads`). Tools therefore never use more cores than requested, even on schedulers that do not confine CPUs.

The container receives the scheduler's job variables (`LSB_*`, `SLURM_*`, `PBS_*`, `JOB_ID`, `SGE_TASK_ID`, `NSLOTS`) through the engine's environment whitelist. Every multi-threaded process also sets `OMP_NUM_THREADS`, `MKL_NUM_THREADS` and `OPENBLAS_NUM_THREADS` to `task.cpus`, which keeps R and numpy linear algebra within the request.

On LSF sites with `LSB_RESOURCE_ENFORCE="cpu"` (seadragon), a task only sees the cores it requested; the smoke test reports this (S06 shows `nproc` equal to the request).

## Job arrays, queue size and submission rate

- **Job arrays.** Processes with many short identical tasks (per-chromosome or per-chunk steps, label `process_array`) are submitted as job arrays of up to `--job_array_size` tasks (default 50). This reduces scheduler load and queue wait. Arrays work with LSF, SLURM, PBS Pro and SGE; the `local` executor ignores them. `--job_array_size 0` turns arrays off.
- **`--queue_size`** (default 100; seadragon 200) is the maximum number of jobs Nextflow keeps queued or running at once. Lower it if your site limits jobs per user.
- **`--submit_rate_limit`** (default `'10 sec'`, meaning 10 submissions per second) throttles submission. Use `'1 sec'` or `'30 min'`-style values only if the scheduler administrators ask.
- Polling: status is checked every 30 s and the queue every 2 min, which is gentle on the scheduler.

## GPUs

GPUs are optional. All tensorQTL steps run on CPU by default, which takes longer but needs nothing special. To use GPUs:

1. Use the CUDA image (`xqtl-federated_<version>-cuda.sif`) through `--sif`.
2. Add `--use_gpu true`.
3. Set the GPU queue and request syntax in your site file:

| Scheduler | `gpu_cluster_options` example |
|---|---|
| LSF | `-gpu num=1` (seadragon) or `-gpu "num=1:mode=shared"` |
| SLURM | `--gres=gpu:1` or `--gpus=1` |
| PBS Pro | `-l select=1:ngpus=1` |
| SGE | `-l gpu=1` |

With `--use_gpu`, tasks labelled `process_gpu_capable` go to `--gpu_queue`, request 4 CPUs plus the GPU flags, and Apptainer/Singularity are started with `--nv` so the container can see the NVIDIA driver. Without it, the same tasks run on CPU with 8 CPUs. Results are numerically equivalent within floating-point tolerance. The smoke test checks GPU tensorQTL only when `--use_gpu true` is given.

## Container runtime

### Apptainer or Singularity

- **Apptainer ≥ 1.0** is recommended (`-profile <scheduler>,apptainer`).
- On sites where Apptainer is installed under the name `singularity`, or where only SingularityCE ≥ 3.8 is available, use `-profile <scheduler>,singularity`. Seadragon is in this group (Apptainer 1.5 invoked as `singularity`).
- Docker is used only for CI (`-profile test,docker`).

The image is a local file. `run_xqtl.sh` passes its path as `--sif`, so nothing is pulled at run time and no image cache is needed (`NXF_APPTAINER_CACHEDIR` and `NXF_SINGULARITY_CACHEDIR` are unused).

### Bind paths

Inside the container, only paths that are bound are visible. `autoMounts` binds the work directory and the paths of input files automatically, but data reached through symlinks or listed inside a samplesheet often live on other filesystems. List every filesystem that holds FASTQs, genotypes, the restricted file, the work directory and outputs:

```text
--bind_paths '/scratch,/project,/data'
```

The smoke test checks a path you choose with `--bind_test_path /project/somefile.txt` (S09).

### Work directory, temporary space and cleanup

- **Work directory (`-w`).** Every task runs in its own subdirectory of the work directory. Put it on a large shared filesystem visible to all compute nodes, ideally scratch: `-w /scratch/$USER/xqtl_work`. A full cohort produces about 2 TB of transient alignment files at peak.
- **Temporary files.** Tasks write temporary files inside their work directory, not to `/tmp`, because `/tmp` on compute nodes is often small. To use fast node-local disk instead, set `process.scratch = '/local/scratch'` (or `'$TMPDIR'`) in your site file; Nextflow then runs each task in that directory and copies outputs back.
- **Apptainer build/temp directories** are not used at run time, since nothing is built or pulled.
- **Cleanup.** Keep the work directory until the run has finished and the export has been reviewed, because `-resume` needs it. Then delete it, or run `nextflow clean -f -before <last run name>` from the launch directory. Intermediate BAMs are deleted by the pipeline once their counts are written.

## Running offline

Compute nodes at most sites have no internet access, and the pipeline never needs it:

- **`NXF_OFFLINE=true`** is set by `run_xqtl.sh`. Without it, Nextflow tries to contact `www.nextflow.io` at every launch; on seadragon that added about 130 s per launch before timing out.
- **No plugins.** The pipeline uses no Nextflow plugins, so nothing is downloaded at launch. Input validation is done by a Python step inside the image (`bin/xqtl_validate.py`).
- **The launcher is exported from the image.** `apptainer run <sif> bootstrap ./xqtl` copies the pipeline code, a Java runtime and the Nextflow launcher out of the image, and writes `run_xqtl.sh` with `JAVA_HOME`, `NXF_HOME` and the image path set. The host needs only Apptainer or Singularity and a scheduler.
- **References are inside the image.** Processes read the public references from `/opt/xqtl/refs` inside the container. Indices that must be built for your read length are written once by `--mode prepare_references` to `--site_refs`.

## Running the head inside a batch job

Running Nextflow on a login node is fine for short runs, but many sites kill long login-node processes. For full runs, submit the head as an ordinary batch job; it then submits tasks as further jobs. Phase 0 confirmed this works on seadragon (head job on the `medium` queue).

Templates are in `assets/head_job/` (`lsf.sh`, `slurm.sh`, `pbs.sh`, `sge.sh`). The essential settings are:

| Setting | Value | Why |
|---|---|---|
| CPUs | 1–2 | The head does little computation |
| Memory | 4–8 GB | Java heap plus headroom |
| `NXF_OPTS` | `-Xms1g -Xmx4g` | Keeps the JVM within the job's memory limit |
| Walltime | Longest queue (e.g. 120 h) | The head must outlive all tasks; resume if it is killed |
| Queue | One that allows jobs to submit jobs | Some sites forbid submission from compute nodes; ask your administrators |

Example for LSF:

```bash
#BSUB -J xqtl_head
#BSUB -q long
#BSUB -W 120:00
#BSUB -n 1
#BSUB -M 8
#BSUB -R "rusage[mem=8]"
#BSUB -o xqtl_head.%J.out
export NXF_OPTS="-Xms1g -Xmx4g"
cd "$LS_SUBCWD"
./xqtl/run_xqtl.sh -profile seadragon -params-file cohort.yml --mode site -resume
```

If the head job is killed, resubmit it with `-resume`: finished tasks are reused, and jobs that were still running are resubmitted.

## Site-specific settings

### Account, project and QOS codes

Put flags required on every job in `--cluster_options`, usually in your site file:

| Scheduler | Example |
|---|---|
| LSF | `-P myproject` or `-G mygroup` |
| SLURM | `--account=myproj --qos=normal` |
| PBS | `-A myproject` |
| SGE | `-P myproject` |

GPU tasks receive both `cluster_options` and `gpu_cluster_options`.

### Module systems

Nothing needs to be loaded inside tasks, because all software is in the image. On the head node, Apptainer or Singularity must be on `PATH` for the head and for the batch jobs. If your site provides it as a module, add `module load apptainer` to your shell profile or to the head job script, and add `process.beforeScript = 'module load apptainer'` to your site file if compute nodes do not inherit it.

## Validating a new site

Before running real data, run the smoke test with exactly the profile you will use for real runs:

```bash
./xqtl/run_xqtl.sh -profile seadragon,smoke --bind_test_path /rsrch5/home/$USER/somefile.txt \
    --restricted /path/to/xqtl-restricted-1.tar.gz
```

For other sites: `-profile slurm,apptainer,smoke -c conf/sites/mysite.config`. Then run the same command again with `-resume`; every task should be reported as cached. The test takes 5–15 minutes, depending on queue wait, and uses only synthetic data generated inside the tasks.

Results are written to `xqtl_smoke/smoke/smoke_summary.tsv`, ending with `VERDICT GO` or `NO-GO`:

| Test | What it checks | If it fails |
|---|---|---|
| S01 container R | The image starts on compute nodes; R 4.3.1 runs | Engine not on `PATH` on compute nodes; wrong `--sif` path; image not readable from nodes |
| S02 embedded tools | Every tool in the image reports a version (`tool_versions.tsv`) | Wrong or damaged image; re-verify the sha256 |
| S03 embedded references | Reference files match the image manifest | Wrong image or `--refs_dir` |
| S04 restricted file | The restricted file opens, has the expected version and checksums | WARN if not given; FAIL if damaged or the wrong version. Request it by email (see `SETUP.md`) |
| S05 plink2 | Genotype tooling runs | Image problem |
| S06 resources | 4 CPUs and 8 GB are granted; a 4 GB allocation succeeds; whether CPUs are confined (scheduler view in `sched_view.txt`) | Memory units or per-job/per-slot setting wrong (see [Memory](#memory)) |
| S07 fan-out | 5 independent jobs run, usually on several nodes | Queue or account settings |
| S08 job array | A 5-task job array runs | Arrays disabled at your site: set `--job_array_size 0` |
| S09 bind read | A file on your project filesystem is visible in the container | Add the filesystem to `--bind_paths` |
| S10 memory enforcement and retry | A task that exceeds its 2 GB request is killed and retried with 4 GB | WARN means memory is not enforced; set `--max_memory` conservatively |
| S11 tensorQTL | CPU (or GPU with `--use_gpu`) cis-mapping recovers a planted signal | Image or GPU driver problem |
| S12 Storey q-value | The R `qvalue` bridge runs on tensorQTL output | Image problem |

WARN and SKIP lines do not block a GO, but they tell you which defaults to adjust. Send the summary file to MDACC with your first export.

## Worked example: seadragon (MD Anderson)

Seadragon runs IBM LSF with Apptainer 1.5.1, invoked as `singularity`. The built-in profile `conf/sites/seadragon.config` encodes the Phase 0 results (appendix *Seadragon validation*):

```groovy
includeConfig '../executors/lsf.config'
includeConfig '../containers/singularity.config'
params {
    queue_rules         = 'short:3h,medium:23h59m,long:120h'
    gpu_queue           = 'gpu'
    gpu_cluster_options = '-gpu num=1'
    bind_paths          = '/rsrch5,/rsrch9'
    max_cpus            = 16
    max_memory          = '120.GB'
    max_time            = '120.h'
    queue_size          = 200
}
executor { perJobMemLimit = true }
```

Why each value:

- **Queues.** Seadragon's submission filter (esub) rejects walltimes outside each queue's range: `short` ≤ 3 h, `medium` more than 3 h and under 24 h, `long` 24–120 h. The rule `medium:23h59m` keeps 24-hour tasks out of `medium`.
- **Memory.** `LSF_UNIT_FOR_LIMITS=GB` and `LSB_MEMLIMIT_ENFORCE=y`: memory is per job, in GB, and enforced. Phase 0 saw the over-limit task killed with exit 130 and retried successfully with 4 GB.
- **CPUs.** `LSB_RESOURCE_ENFORCE="cpu gpu"` confines each job to its requested cores (Phase 0: `nproc=4` for a 4-core request).
- **Bind paths.** Project data live on `/rsrch5` and `/rsrch9`.
- **Network.** Compute nodes have no internet. Of the container registries, `ghcr.io` and `quay.io` were reachable from login nodes at about 0.6 MB/s; Docker Hub, `gcr.io` and Amazon ECR were blocked or unreliable. Download the SIF once (about 3–4 h at that speed) or copy it by Globus.
- **Building images.** `apptainer build --fakeroot` with package installation fails on seadragon (apt exit 100). Never build images there; use the published SIF.
- **Head.** Both a login-node head and a head inside an LSF job (`medium` queue, 1 core, 4 GB) passed Phase 0.

## Tuning table

Approximate requirements for one cohort of about 350 samples with about 50 M read pairs each:

| Stage | Class | Jobs | Typical CPU-hours | Peak disk |
|---|---|---|---|---|
| Index build (`prepare_references`, once) | `process_high` + `process_high_memory` | 3 | 30–60 | 60 GB |
| STAR + Salmon per sample | `process_medium` | ~700 | 1,000–1,800 | ~2 TB BAMs (transient) |
| Junctions, featureCounts, editing pileups, IRFinder-S | `process_low`/`process_medium` | ~1,400 | 400–700 | 200 GB |
| Aggregation, normalization, MuSiC | `process_medium`/`process_high` | ~40 | 50–100 | 50 GB |
| HCP grid (per modality × ancestry) | `process_medium` | ~60 | ~60 | small |
| tensorQTL (per modality × ancestry) | `process_gpu_capable` | ~50 | 100–300 (CPU) | 100 GB |
| GxE tier 1 (per modality × exposure) | `process_gpu_capable` + `process_long` | ~20 | 100–300 (CPU) | small |
| G×C, GxE tier 2, LD, export | `process_medium` | ~100 | 20–50 | 100–150 GB export |

To shorten runs, raise `--queue_size` (if allowed), use GPUs for tensorQTL, and place the work directory on the fastest shared filesystem.

## HPC troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Every task fails at once with "command not found: singularity/apptainer" | Engine not on `PATH` on compute nodes | `process.beforeScript = 'module load apptainer'` in the site file |
| Jobs rejected at submission ("walltime exceeds queue limit", esub errors) | Queue rules do not match site limits | Fix `--queue_rules`; check `bqueues -l` / `scontrol show partition` |
| Jobs pending forever | Request larger than any node, or wrong queue/account | Lower `--max_memory`/`--max_cpus`; set `--cluster_options` |
| Multi-core LSF tasks killed far below their memory request | Per-slot memory interpretation | Keep `--lsf_per_job_mem true`; check `LSF_UNIT_FOR_LIMITS` |
| Exit 130/137/140 repeated on all attempts | Class too small for your data | Raise the class memory in the site file |
| "No such file" for an input that exists | Filesystem not bound in the container | Add it to `--bind_paths` |
| Each launch hangs about 2 minutes before starting | Nextflow trying to reach the internet | Use `run_xqtl.sh` (sets `NXF_OFFLINE=true`) |
| Head killed on the login node | Login-node process limits | Run the head as a batch job (see [above](#running-the-head-inside-a-batch-job)) and `-resume` |
| `-resume` reruns everything | Different launch directory or deleted work directory | Launch from the same directory with the same `-w` |
| Job arrays rejected | Arrays disabled or capped | `--job_array_size 0` (or a smaller size within the site cap) |
| GPU tasks fail with CUDA errors | Driver not visible or CPU image used | Use the `-cuda` SIF; check `nvidia-smi` on the GPU node; `--nv` is added automatically |
| Scheduler complains about submission rate | Too many submissions per second | `--submit_rate_limit '2 sec'` and lower `--queue_size` |
